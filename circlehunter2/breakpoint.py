from collections import defaultdict
from collections import deque

from natsort import natsorted
from scipy import stats
import networkx as nx
import numpy as np
import pysam

from circlehunter2.bamReader import fetch_breakpoints
from circlehunter2.bamReader import count_cross_reads


DISCORDANT = "DISCORDANT"
CONTINUOUS = "CONTINUOUS"

NEXT_EDGE_TYPE = {
    DISCORDANT: CONTINUOUS,
    CONTINUOUS: DISCORDANT
}


class BreakpointEstimator:
    def __init__(self, bam_file, error=5, extend=500,
                 flag_include=1, flag_exclude=1548,
                 mapq=10, mismatch_ratio=0.02,
                 barcode_tag="CB", buffer_size=256):
        self._bam = pysam.AlignmentFile(bam_file, "rb")
        self._error = error
        self._extend = extend
        self._fetch_args = dict(
            flag_include=flag_include, flag_exclude=flag_exclude,
            mapq=mapq, mismatch_ratio=mismatch_ratio,
            barcode_tag=barcode_tag
        )
        self._buffer_size = buffer_size

    @staticmethod
    def to_estimate_axis(start, end, reverse, data):
        if reverse:
            data = end - data + 1
        else:
            data -= start + 1
        return data

    @staticmethod
    def to_chromosome_axis(start, end, reverse, data):
        if reverse:
            data = end - data + 1
        else:
            data = start + data - 1
        return data

    def discordant_likelihood(self, hypo, obs):
        """
        Calculate the likelihood of each breakpoint hypo when observe a discordant reads end at obs.
        """
        # uniform distribution
        pdf = 1 / hypo
        # this is used for error, make the pdf shape symmetrical in obs
        pdf[hypo < obs] = 1 / (2 * obs - hypo[hypo < obs])
        # breakpoint can't be smaller than obs discordant read end
        pdf[hypo < obs - self._error] = 0
        # return log pdf and avoid overflow
        return np.log(pdf + np.exp(-10))

    def discordant_posterior(self, hypo, samples):
        """
        Calculate posterior for given hypothesis when observe samples.
        """
        if len(samples) == 0:  # no sample observed, set all posterior to minimum
            return np.full_like(hypo, -10, dtype=np.float64)
        # we will observe same end positions, this will help us cache result for same obs
        likelihood = {}
        # we have log the pdf, so add all log pdf of observe samples is the posterior
        # log pdf will avoid times operation which will overflow when the numbers all smaller than 1
        posterior = np.zeros_like(hypo, dtype=np.float64)
        for sample in samples:
            posterior += likelihood.get(sample,
                                        self.discordant_likelihood(hypo, sample))
        return posterior

    def clipped_likelihood(self, hypo, obs):
        """
        Calculate likelihood for given hypothesis when observe a clipped reads end at obs.
        """
        pdf = stats.norm.pdf(hypo, obs, self._error / 2)
        return np.log(pdf + np.exp(-10))

    def clipped_posterior(self, hypo, samples):
        """
        Calculate posterior for given hypothesis when observe samples.
        """
        if len(samples) == 0:  # no sample observed, set all posterior to minimum
            return np.full_like(hypo, -10, dtype=np.float64)
        # we will observe same end positions, this will help us cache result for same obs
        likelihood = {}
        posterior = np.zeros_like(hypo, dtype=np.float64)
        for sample in samples:
            # we have log the pdf, so add all log pdf of observe samples is the posterior
            # log pdf will avoid times operation which will overflow when the numbers all smaller than 1
            posterior += likelihood.get(sample,
                                        self.clipped_likelihood(hypo, sample))
        return posterior

    def estimate(self, chrom1, start1, end1, reverse1, chrom2, start2, end2, reverse2):
        # observe
        discordant, clipped = fetch_breakpoints(
            self._bam,
            chrom1, start1, end1, reverse1,
            chrom2, start2, end2, reverse2,
            **self._fetch_args, buffer_size=self._buffer_size
        )
        discordant_obs = self.to_estimate_axis(
            start1, end1, reverse1, discordant)
        clipped_obs = self.to_estimate_axis(start1, end1, reverse1, clipped)
        # hypo
        hypo = np.arange(0, end1 - start1) + 1
        # posterior
        discordant_posterior = self.discordant_posterior(hypo, discordant_obs)
        clipped_posterior = self.clipped_posterior(hypo, clipped_obs)
        log_posterior = discordant_posterior + clipped_posterior
        # shift back to positive range and exp
        posterior = np.exp(log_posterior - np.max(log_posterior))
        # grid search maximum likelihood estimation
        mle = hypo[np.argmax(posterior)]
        # calculate cdf and grid search 95% CI
        cdf = np.cumsum(posterior / np.sum(posterior))
        ci = hypo[np.sum(np.stack([cdf < 0.025, cdf < 0.975]), axis=1)]
        # shift back to chromosome axis
        mle = self.to_chromosome_axis(start1, end1, reverse1, mle)
        if reverse1:
            cr, cl = self.to_chromosome_axis(start1, end1, reverse1, ci)
        else:
            cl, cr = self.to_chromosome_axis(start1, end1, reverse1, ci)
        return cl, mle, cr

    def count_cross_reads(self, chrom, pos, reverse):
        return count_cross_reads(
            self._bam,
            chrom, pos, reverse,
            error=self._error, extend=self._extend,
            **self._fetch_args
        )


class PairedPeakGraph(nx.MultiGraph):
    def __init__(self, incoming_graph_data=None):
        super().__init__(incoming_graph_data=incoming_graph_data)
        self._precised = False

    def add_discordant_edge(self, chrom1, start1, end1, reverse1, chrom2, start2, end2, reverse2, count):
        """
        Add a discordant type edge into the breakpoint graph.
        """
        self.add_edge(
            (chrom1, start1, end1, reverse1),
            (chrom2, start2, end2, reverse2),
            key=DISCORDANT, count=count
        )

    def estimate_breakpoints(self, bam_file, error=5, extend=500,
                             flag_include=1, flag_exclude=1548,
                             mapq=10, mismatch_ratio=0.02,
                             barcode_tag="CB", buffer_size=256):
        """
        Estimate precise breakpoint for each breakpoint in this graph.
        Estimate result will add to the node attr.
        """
        estimator = BreakpointEstimator(
            bam_file, error, extend,
            flag_include, flag_exclude,
            mapq, mismatch_ratio,
            barcode_tag, buffer_size
        )
        for peak1, peak2, _ in self.edges:
            cl, mle, cr = estimator.estimate(*peak1, *peak2)
            cross, disjoint = estimator.count_cross_reads(
                peak1[0], mle, peak1[3]
            )
            nx.set_node_attributes(self, {peak1: dict(
                cl=cl, mle=mle, cr=cr, cross=cross, disjoint=disjoint
            )})

            cl, mle, cr = estimator.estimate(*peak2, *peak1)
            cross, disjoint = estimator.count_cross_reads(
                peak2[0], mle, peak2[3]
            )
            nx.set_node_attributes(self, {peak2: dict(
                cl=cl, mle=mle, cr=cr, cross=cross, disjoint=disjoint
            )})
        self._precised = True

    def to_breakpoint_graph(self):
        """
        Transform the graph to BreakpointGraph.
        Estimate breakpoint will act as node, paired peak info save to node attr.
        """
        assert self._precised

        graph = BreakpointGraph()

        # loop through edges
        for (chrom1, start1, end1, reverse1), (chrom2, start2, end2, reverse2), edge_type in self.edges:
            if edge_type != DISCORDANT:
                continue

            data1 = self.nodes[(chrom1, start1, end1, reverse1)]
            data2 = self.nodes[(chrom2, start2, end2, reverse2)]
            count = self.edges[(
                (chrom1, start1, end1, reverse1),
                (chrom2, start2, end2, reverse2),
                DISCORDANT
            )]["count"]

            graph.add_discordant_edge(
                chrom1, data1["mle"], reverse1, data1["cl"], data1["cr"],
                start1, end1, data1["cross"], data1["disjoint"],
                chrom2, data2["mle"], reverse2, data2["cl"], data2["cr"],
                start2, end2, data2["cross"], data2["disjoint"],
                count=count
            )

        return graph


class BreakpointGraph(nx.MultiGraph):
    def add_discordant_edge(self,
                            chrom1, pos1, reverse1, cl1, cr1,
                            peak_start1, peak_end1, cross1, disjoint1,
                            chrom2, pos2, reverse2, cl2, cr2,
                            peak_start2, peak_end2, cross2, disjoint2,
                            count):
        """
        Add a discordant type edge into the breakpoint graph.
        """
        self.add_node(
            (chrom1, pos1, reverse1),
            cl=cl1, cr=cr1,
            peak_start=peak_start1, peak_end=peak_end1,
            cross=cross1, disjoint=disjoint1
        )
        self.add_node(
            (chrom2, pos2, reverse2),
            cl=cl2, cr=cr2,
            peak_start=peak_start2, peak_end=peak_end2,
            cross=cross2, disjoint=disjoint2
        )
        self.add_edge(
            (chrom1, pos1, reverse1),
            (chrom2, pos2, reverse2),
            key=DISCORDANT, count=count
        )

    def add_continuous_edge(self, chrom, start, end, depth, coverage):
        """
        Add a continuous type edge into the breakpoint graph.
        """
        self.add_edge(
            (chrom, start, True), (chrom, end, False), key=CONTINUOUS,
            depth=depth, coverage=coverage
        )

    def all_possible_continuous_edges(self, cross_fraction=0.05):
        """
        Yield all possible continuous edges.
        """
        nodes = natsorted(list(self.nodes))
        for i in range(len(nodes)):
            if not nodes[i][2]:  # first node should be reverse
                continue
            for j in range(i + 1, len(nodes)):
                if nodes[j][2]:  # second node should not be reverse
                    continue
                if nodes[i][0] != nodes[j][0]:  # same chromosome
                    continue
                yield nodes[i][0], nodes[i][1], nodes[j][1]
                # is this breakpoint continuously?
                total = (
                    self.nodes[nodes[j]]["cross"] +
                    self.nodes[nodes[j]]["disjoint"]
                )
                if total==0:
                    # end if division by zero.
                    # break # report error
                    fraction = 0 # test
                else:
                    fraction = self.nodes[nodes[j]]["cross"] / total
                if fraction < cross_fraction:  # cross is too small
                    # end looking for right breakpoint(j) for this left breakpoint(i)
                    break
            # end j loop
        # end i loop

    def breakpoint_degree(self, node):
        """
        Calculate degree of two type edges for the given node.
        """
        degree = {
            DISCORDANT: 0, CONTINUOUS: 0
        }

        # loop through nodes
        for _, meta in self.adj[node].items():
            # loop through edge types
            for key in meta.keys():
                degree[key] += 1

        return degree

    def drop_invalid_nodes(self):
        """
        A valid node (breakpoint) should have at least one discordant type edge
        and one continuous type edge. Remove those node that is invalid.
        """
        dropped_nodes = set()
        nodes = list(self.nodes)

        # loop through nodes
        while nodes:
            node = nodes.pop()

            # skip dropped nodes
            if node in dropped_nodes:
                continue

            # calculate degree
            degree = self.breakpoint_degree(node)
            # is this a invalid node?
            if any(value < 1 for value in degree.values()):
                # note: the adjacency node may become invalid
                # we should them check again
                # push them back to the nodes queue
                for adj in self.adj[node]:
                    nodes.append(adj)
                # remove invalid node
                self.remove_node(node)
                dropped_nodes.add(node)

        # not useful
        return len(dropped_nodes)

    def sorted_continuous_edges(self):
        """
        Return sorted continuous edges by length.
        """
        edges = []
        for (node1, node2, key), _ in self.edges.items():
            if key == CONTINUOUS:
                node1, node2 = sorted((node1, node2))
                edges.append(
                    (node2[1] - node1[1], (node1, node2))
                )
        edges.sort()
        return (node for *_, node in edges)

    def find_discordant_children(self, father):
        """
        Find next node which is connected by a discordant type edge.
        """
        for adj, meta in self.adj[father].items():
            for key in meta.keys():
                if key == DISCORDANT:
                    yield adj

    def find_continuos_children(self, father):
        """
        Find next node which is connected by a continuous type edge.
        """
        for adj, meta in self.adj[father].items():
            for key in meta.keys():
                if key == CONTINUOUS:
                    yield adj

    def find_children(self, father, last_edge_type):
        """
        Find next node which can connected follow by last_edge_type
        """
        if last_edge_type == DISCORDANT:
            return self.find_continuos_children(father)
        elif last_edge_type == CONTINUOUS:
            return self.find_discordant_children(father)
        else:
            raise ValueError(f"invalid last_edge_type: {last_edge_type}")

    def find_circles(self):
        """
        Yield circles in this graph.
        """
        edges = self.sorted_continuous_edges()

        visited = set()

        while True:
            try:
                edge = next(edges)
            except StopIteration:
                break

            # skip visited edge as a source
            if edge in visited:
                continue
            else:
                visited.add(edge)

            circle_stack = list(edge)
            edge_stack = [CONTINUOUS, ]
            node_pointer = {node: i for i, node in enumerate(circle_stack)}

            stack = [
                self.find_discordant_children(circle_stack[-1])
            ]

            while stack:
                children = stack[-1]

                try:
                    child = next(children)
                except StopIteration:  # no more child for the last node
                    pop = circle_stack.pop()
                    node_pointer.pop(pop)
                    edge_stack.pop()
                    stack.pop()
                    continue

                # print(node_pointer, child)
                # print(edge_stack)

                if child not in node_pointer:  # new peace of circle, add to stack
                    node_pointer[child] = len(circle_stack)
                    circle_stack.append(child)
                    edge_stack.append(NEXT_EDGE_TYPE[edge_stack[-1]])
                    stack.append(self.find_children(child, edge_stack[-1]))
                    # remember visited continuous edge to skip in future
                    if edge_stack[-1] == CONTINUOUS:
                        edge = tuple(sorted(
                            (circle_stack[-1], circle_stack[-2])
                        ))
                        visited.add(edge)
                # this is not a valid ecDNA
                # Add filters if BFB or inversion ? 
                elif len(edge_stack) == node_pointer[child] or NEXT_EDGE_TYPE[edge_stack[-1]] == edge_stack[node_pointer[child]]:
                    continue
                # this is a valid ecDNA
                else:
                    yield circle_stack[node_pointer[child]:], edge_stack[node_pointer[child]:] + [NEXT_EDGE_TYPE[edge_stack[-1]]]

    @staticmethod
    def sorted_circle(nodes):
        """
        Sort the nodes in a circle to make it unique.
        """
        # the minimum node is the first node
        first = min(nodes)
        first_index = nodes.index(first)
        last_index = first_index - 1
        # reconstruct the cycle to make is fit the order but wont break the cycle
        next_index = first_index + 1 if first_index + 1 != len(nodes) else 0
        if nodes[last_index] > nodes[next_index]:
            nodes = nodes[first_index:] + nodes[:first_index]
        else:
            nodes = nodes[:next_index][::-1] + nodes[next_index:][::-1]
        return tuple(nodes)

    def find_uniq_circles(self):
        """
        Yield unique circles in this graph.
        """
        visited = set()
        for nodes, links in self.find_circles():
            # skip visited circle
            sorted_nodes = self.sorted_circle(nodes)
            if sorted_nodes in visited:
                continue
            visited.add(sorted_nodes)
            yield nodes, links

    def find_ecDNAs(self):
        """
        Yield all valid ecDNA in this graph.
        """
        # loop through circles
        for nodes, links in self.find_uniq_circles():
            segments = []
            # rolling through segments
            nodes, links = deque(nodes), deque(links)
            begin = nodes[0]
            while True:
                if links[0] == CONTINUOUS:  # this is a ecDNA segment
                    segments.append((nodes[0], nodes[1]))
                # rolling to next node
                nodes.rotate(-1)
                links.rotate(-1)
                if nodes[0] == begin:  # end of circle
                    break
            yield segments
