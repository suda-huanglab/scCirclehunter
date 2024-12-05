from argparse import ArgumentParser
import logging
import os
import re

from scipy.stats import poisson
from scipy.sparse import csr_matrix
import anndata as ad
import pandas as pd

from circlehunter2.bamReader import build_tag_track_multi_process
from circlehunter2.bamReader import build_tag_track_single_process
from circlehunter2.bamReader import collect_chromosome_sizes
from circlehunter2.bedTrack import BedTrack
from circlehunter2.breakpoint import PairedPeakGraph
from circlehunter2.breakpoint import DISCORDANT, CONTINUOUS


def run(bam_file, output,
        min_insert_size=5000,
        flag_include=1, flag_exclude=1548,
        mapq=10, mismatch_ratio=0.02,
        barcode_tag="CB",
        buffer_size=256,
        processes=2, window_size=64_000_000,
        no_progressbar=True,
        blacklist_file=None,
        cutoff=5, extend=500, error=5, cross=0.05):

    logger = logging.getLogger("circlehunter2")

    # check if processes number is valid
    processes = min(processes, os.cpu_count())

    # reading BAM input
    if processes > 1:
        logger.info(f"Reading BAM input with {processes} process ...")
        tag_collection, paired_tag_collection = build_tag_track_multi_process(
            bam_file,
            min_insert_size=min_insert_size,
            flag_include=flag_include, flag_exclude=flag_exclude,
            mapq=mapq, mismatch_ratio=mismatch_ratio,
            barcode_tag=barcode_tag, buffer_size=buffer_size,
            processes=processes, window_size=window_size,
            no_progressbar=no_progressbar
        )
    else:
        logger.info("Reading BAM input with single process ...")
        tag_collection, paired_tag_collection = build_tag_track_single_process(
            bam_file,
            min_insert_size=min_insert_size,
            flag_include=flag_include, flag_exclude=flag_exclude,
            mapq=mapq, mismatch_ratio=mismatch_ratio,
            barcode_tag=barcode_tag, buffer_size=buffer_size,
            no_progressbar=no_progressbar
        )

    logger.info(f"Filtering tags ...")

    # filter tags by chromosome sizes
    all_chrom_sizes = collect_chromosome_sizes(bam_file)
    chrom_sizes = {
        chrom: size
        for chrom, size in all_chrom_sizes.items()
        if re.match(r"chr[0-9XY]{1,2}", chrom)
    }

    tag_collection.set_chromosome_sizes(chrom_sizes)
    paired_tag_collection.set_chromosome_sizes(chrom_sizes)

    # pool all paired tags
    all_paired_tags = paired_tag_collection.pool_barcodes(
        list(paired_tag_collection.get_barcodes())
    )

    # filter tags by blacklist
    if blacklist_file is not None:
        blacklist = BedTrack.from_bed(
            all_chrom_sizes, blacklist_file, buffer_size=buffer_size
        )
        blacklist.set_chromosome_sizes(chrom_sizes)

        overlap = all_paired_tags.overlap(blacklist)
        keep = (overlap["overlaps1"] == 0) & (overlap["overlaps2"] == 0)
        all_paired_tags.filter(keep)

    logger.info(f"A total of {tag_collection.total} tag is passed filter")
    logger.info(
        f"A total of {paired_tag_collection.total} discordant tag is passed filter"
    )

    logger.info(f"Pileup tags ...")

    # pileup tags as depth
    all_pileup = tag_collection.pool_barcodes(
        list(tag_collection.get_barcodes())
    )
    
    # # determine cutoff
    # mean = all_pileup.mean()
    # cutoff_mean = poisson.isf(0.05, mean)
    #     # use inverse survival of mean depth of the whole genome if the value is larger than 5.

    # if cutoff is None:
    #     cutoff = cutoff_mean  # if None, mean coverage.
    # else:
    #     cutoff = max(cutoff, cutoff_mean)  # if larger than 5 (default), use mean.

    # determine cutoff
    if cutoff is None:
        # no custom cutoff, use inverse survival of mean depth of the whole genome
        mean = all_pileup.mean()
        cutoff = poisson.isf(0.05, mean)

    logger.info(f"Cutoff {cutoff} is used in high coverage judge")

    # find all paired peaks
    all_paired_peaks = all_paired_tags.call_peaks(
        cutoff=cutoff, distance=extend
    )
    all_paired_peaks.keep_bi_peaks()

    logger.info(f"Building breakpoint graph ...")

    graph = PairedPeakGraph(all_paired_peaks.to_graph_data())
    graph.estimate_breakpoints(
        bam_file=bam_file, error=error, extend=extend,
        flag_include=flag_include, flag_exclude=flag_exclude,
        mapq=mapq, mismatch_ratio=mismatch_ratio,
        barcode_tag=barcode_tag, buffer_size=buffer_size
    )
    graph = graph.to_breakpoint_graph()
    for edge in graph.all_possible_continuous_edges(cross_fraction=cross):
        depth, coverage = all_pileup.coverage(*edge, cutoff=cutoff)
        graph.add_continuous_edge(*edge, depth=depth, coverage=coverage)

    graph.drop_invalid_nodes()

    logger.info(f"Finding ecDNA by DFS in breakpoint graph ...")

    ecDNAs = list(graph.find_ecDNAs())

    logger.info(f"A total of {len(ecDNAs)} ecDNA is found")

    header = (
        "chrom", "start", "end","name", "score", "strand",
        "start_ci", "end_ci", "start_peak", "end_peak",
        "start_cross", "end_cross",
        "linked_reads", "depth_mean", "high_coverage","width"
    )

    with open(output, "w") as f:
        f.write("#{}\n".format("\t".join(header)))
        for n, ecDNA in enumerate(ecDNAs, start=1):
            total = len(ecDNA)
            for p, (node1, node2) in enumerate(ecDNA, start=1):
                # bed region info
                chrom = node1[0]
                start = node1[1] if node1[2] else node2[1]
                end = node2[1] if node1[2] else node1[1]
                width = end-start
                strand = "+" if node1[2] else "-"
                name = f"ecDNA_{n}_{total}_{p}"

                # breakpoint info
                start_ci = f"{graph.nodes[node1]['cl']}-{graph.nodes[node1]['cr']}"
                end_ci = f"{graph.nodes[node2]['cl']}-{graph.nodes[node2]['cr']}"
                start_peak = f"{graph.nodes[node1]['peak_start']}-{graph.nodes[node1]['peak_end']}"
                end_peak = f"{graph.nodes[node2]['peak_start']}-{graph.nodes[node2]['peak_end']}"
                start_cross = f"{graph.nodes[node1]['cross']}/{graph.nodes[node1]['disjoint']}"
                end_cross = f"{graph.nodes[node2]['cross']}/{graph.nodes[node2]['disjoint']}"
                if strand == "-":
                    start_ci, end_ci = end_ci, start_ci
                    start_peak, end_peak = end_peak, start_peak
                    start_cross, end_cross = end_cross, start_cross

                # discordant info
                current_node = p - 1
                if p == total:
                    next_node = 0
                else:
                    next_node = p

                linked_reads = graph.edges[
                    (ecDNA[current_node][1], ecDNA[next_node][0], DISCORDANT)
                ]["count"]

                # coverage info
                depth = graph.edges[(node1, node2, CONTINUOUS)]['depth']
                coverage = graph.edges[(node1, node2, CONTINUOUS)]['coverage']

                f.write(
                    f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\t"
                    f"{start_ci}\t{end_ci}\t{start_peak}\t{end_peak}\t"
                    f"{start_cross}\t{end_cross}\t{linked_reads}\t"
                    f"{depth:.6f}\t{coverage:.6f}\t{width}\n"
                )

    logger.info(f"Result saved to: {os.path.abspath(output)}")

    if (barcode_tag is not None) and (len(ecDNAs)>0) :
        logger.info(f"Counting reads for each barcode ...")
        ecDNAs_file = os.path.abspath(output)
        ecDNA = pd.read_table(ecDNAs_file, names=header, comment="#")
        matrix = csr_matrix([
            tag_collection.count(row["chrom"], row["start"], row["end"])
            for _, row in ecDNA.iterrows()
        ])
        adata = ad.AnnData(matrix.T, dtype=matrix.dtype)
        adata.obs.index = tag_collection.get_barcodes()
        adata.obs.index.name = "barcode"
        adata.var = ecDNA
        adata.var_names = ecDNA["name"]

        adata_file = os.path.splitext(ecDNAs_file)[0] + ".h5ad"
        adata.write(adata_file)
        logger.info(f"Read counts saved to: {adata_file}")

    logger.info(f"Finished")


def main():
    parser = ArgumentParser(
        description="Find ecDNA from bam file"
    )

    parser.add_argument("<bam>", help="input bam file")
    parser.add_argument("<bed>", help="output bed file")

    parser.add_argument(
        "-l", dest="min_insert_size", default=5000, type=int,
        help="minimum length of insert size"
    )

    parser.add_argument(
        "-f", dest="include", default=1, type=int,
        help="only include reads with all  of the FLAGs in INT present"
    )
    parser.add_argument(
        "-F", dest="exclude", default=1548, type=int,
        help="only include reads with none of the FLAGS in INT present"

    )
    parser.add_argument(
        "-q", dest="mapq", default=10, type=int,
        help="only include reads with mapping quality >= INT"

    )
    parser.add_argument(
        "-r", dest="ratio", default=0.02, type=float,
        help="only include reads with mismatch ratio <= FLOAT"

    )
    parser.add_argument(
        "-b", dest="barcode_tag", default=None, type=str,
        help=(
            "tag name of barcode, run as single cell mode if provided, "
            "output barcode counts for each ecDNA segments"
        )
    )

    parser.add_argument(
        "--buffer-size", dest="buffer_size", default=256, type=int,
        help="buffer size when reading BAM file"
    )

    parser.add_argument(
        "-p", dest="processes", default=1, type=int,
        help="processes number"
    )
    parser.add_argument(
        "--windows-size", dest="windows_size", default=64_000_000, type=int,
        help="windows size for each process to read"
    )

    parser.add_argument(
        "--no-progressbar", dest="no_progressbar", default=False, action="store_true",
        help="do not show progressbar"
    )

    parser.add_argument(
        "--blacklist", dest="blacklist", default=None, type=str,
        help="path of blacklist file"
    )

    parser.add_argument(
        "-c", dest="cutoff", default=None, type=int,
        help="cutoff of enrichment"
    )
    
    parser.add_argument(
        "--extend", dest="extend", default=500, type=int,
        help="extend size of enrichment region"
    )

    parser.add_argument(
        "--error", dest="error", default=5, type=int,
        help="alignment error"
    )
    parser.add_argument(
        "--cross", dest="cross", default=0.05, type=float,
        help="maximum cross fraction"
    )

    args = vars(parser.parse_args())

    logger = logging.getLogger("circlehunter2")
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(msg)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("Running circlehunter2 ...")

    for key, value in args.items():
        logger.info(f"Params {key} = {value}")

    run(
        bam_file=args["<bam>"], output=args["<bed>"],
        min_insert_size=args["min_insert_size"],
        flag_include=args["include"], flag_exclude=args["exclude"],
        mapq=args["mapq"], mismatch_ratio=args["ratio"],
        barcode_tag=args["barcode_tag"],
        buffer_size=args["buffer_size"],
        processes=args["processes"], window_size=args["windows_size"],
        no_progressbar=args["no_progressbar"],
        blacklist_file=args["blacklist"],
        cutoff=args["cutoff"], extend=args["extend"], error=args["error"],
        cross=args["cross"]
    )


if __name__ == "__main__":
    main()
