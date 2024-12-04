# scCirclehunter
![scCirclehunter overview](https://github.com/Maple-sakura/scCirclehunter/raw/main/img/F1.png)
## Step1 - Identification of ecDNA from scATAC-seq based on a pseudo-bulk algorithm (circlehunter2)

### Installation

```bash
micromamba install -f requirements.txt
micromamba install -f requirements-dev.txt
python setup.py install
```

### Usage

```bash
circlehunter2 -p 16 --blacklist /share/references/hg38/blacklist/hg38.blacklist.sorted.bed /mnt/2w/data2/andy/scATAC-Seq/rawdata/cellranger/GBM4349/outs/possorted_bam.bam /home/andy/Projects/circlehunter2/workspace/dev/data/test.bed
```

input BAM file have to be coordinate sorted, index of the BAM file should present with a `.bai` suffix, more options please check

```bash
circlehunter2 -h
```

### Output

cirlcehunter2 output an extended BED format file with columns:

1. `chrom`: contig name of which the segment come from
2. `start`: segment start
3. `end`: segment end
4. `name`: No. of this segment, format as `ecDNA_{ecDNA_no}_{total_segments}_{segment_no}`, you can use these info to form a complete ecDNA circle
5. `score`: no data in this columns, always a `.`
6. `strand`: orientation of this segment
7. `start_ci`: confidence interval of segment start
8. `end_ci`: confidence interval of segment end
9. `start_peak`: discordant enrich region of segment start
10. `end_peak`: discordant enrich region of segment end
11. `start_cross`: reads counts of which cross the segment start breakpoint and disjoint the breakpoint
12. `end_cross`: reads counts of which cross the segment end breakpoint and disjoint the breakpoint
13. `linked_reads`: reads counts of which connected this segment and the next one
14. `depth_mean`: mean of depth of this region, maybe false positive if this is smaller than `cutoff`[^cutoff]
15. `high_coverage`: fraction of region that depth is higher than `cutoff`, maybe false positive if this is pretty small

[cutoff]: `cutoff` is a user input params `-c` or determine by circlehunter2 automatically if not set by user, which is default to the inverse survival of poisson distribution use mean of the whole genome depth as $\lambda$. `cutoff` will be log in the stderr.

## Step2 - Assigning ecDNA to Cell Populations
Refer to the example in Refer to the example in **[demo/scCirclehunter_downstream_demo.html](https://github.com/Maple-sakura/scCirclehunter/blob/main/demo/scCirclehunter_downstream_demo.html)** for assigning ecDNA to cell populations. The functions utilized are available in the scripts/ directory.
 for assigning ecDNA to cell populations. The functions utilized are available in the **scripts/** directory.

