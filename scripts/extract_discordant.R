### usage.
# Rscript /home/jiangr/Scripts_by_jiangr/ecDNA_simulate/extract_disreads/extract_discordant.R \
#     --circlehunter2_bed /mnt/2w/data1/jiangr/scATAC/GBM/SRR10315837/circlehunter2/SRR10315837_circlehunter2_ecDNA.bed \
#     --cellranger_bam /mnt/2w/data1/jiangr/scATAC/GBM/SRR10315837/SRR10315837_CellRanger/outs/possorted_bam.bam \
#     --outpath /mnt/2w/data1/jiangr/scATAC/GBM/SRR10315837/circlehunter2/tmp \
#     --insert_size 2000

options(warn = -1)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(usage = "Usage: %prog [options]")
parser <- add_option(parser, "--circlehunter2_bed", dest = "circlehunter2_bed", type = "character",
    default = NULL, help = "result bed file of circlehunter2.")
parser <- add_option(parser, "--cellranger_bam", dest = "cellranger_bam", type = "character",
    default = NULL, help = "possorted_bam.bam output by cellranger-atac.")
parser <- add_option(parser, "--outpath", dest = "outpath", type = "character", default = NULL,
    help = "output path to save result file.")
parser <- add_option(parser, "--insert_size", dest = "insert_size", type = "numeric",
    default = 5000, help = "minimal insert size for discordant reads.")

opts <- parse_args(parser)

circlehunter2_bed <- opts$circlehunter2_bed
bam_file <- opts$cellranger_bam
outpath <- opts$outpath
insert_size <- opts$insert_size

setwd(outpath)

circlehunter2_bed <- data.table::fread(circlehunter2_bed, data.table = F)
colnames(circlehunter2_bed) <- gsub(colnames(circlehunter2_bed), pattern = "#", replacement = "")

for (i in c(1:nrow(circlehunter2_bed))) {
    tmp <- circlehunter2_bed[i, c("chrom", "start_peak", "end_peak")]
    name <- circlehunter2_bed[i, "name"]
    ec_chr <- circlehunter2_bed[i, "chrom"]
    out.tmp <- paste0(outpath, "/", name, ".tmp")

    out.bam <- paste0(outpath, "/", name, "_discordantReads.bam")
    out.txt <- paste0(outpath, "/", name, "_discordantcells.txt")
    write.table(x = tmp, file = out.tmp, row.names = F, col.names = FALSE, sep = "\t",
        quote = F)

    cmd1 <- sprintf("cat %s | while IFS= read -r line; do
        chr=$(echo $line | awk '{print $1}')
        start=$(echo $line | awk '{split($2, a, '-'); print a[1] '-' a[2]}')
        end=$(echo $line | awk '{split($3, a, '-'); print a[1] '-' a[2]}')
        samtools view -@ 10 -b %s $chr:$start $chr:$end > %s 
    done",
        out.tmp, bam_file, out.bam)

    cmd2 <- sprintf("samtools view %s | awk '($3 == %s) && (($9 > %d) || ($9 < -%d))' | awk '{for (i=12; i<=NF; i++) {if ($i ~ /^CB:/) {split($i, a, ':'); print a[3];}}}'|sort|uniq > %s",
        out.bam, ec_chr, insert_size, insert_size, out.txt)

    cmd1 <- gsub("\\$line", "\"$line\"", cmd1)
    cmd1 <- gsub("'-'", "\"-\"", cmd1)
    cmd1 <- gsub("\\$chr:\\$start", "\"$chr:$start\"", cmd1)
    cmd1 <- gsub("\\$chr:\\$end", "\"$chr:$end\"", cmd1)
    cmd2 <- gsub("(chr\\d+)", "\"\\1\"", cmd2)
    cmd2 <- gsub("':'", "\":\"", cmd2)

    system(cmd1)
    system(cmd2)
}

system("rm ./*tmp ")