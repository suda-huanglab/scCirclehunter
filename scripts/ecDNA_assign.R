# options(warn=-1)
# suppressPackageStartupMessages(library(ArchR))
# suppressPackageStartupMessages(library(GenomicRanges))
# suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
# suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(readr))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
# suppressPackageStartupMessages(library(biomaRt))
# suppressPackageStartupMessages(library(org.Hs.eg.db))
# suppressPackageStartupMessages(library(ggpubr))
# suppressPackageStartupMessages(library(diptest))
# suppressPackageStartupMessages(library(Rsamtools))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(GenomicAlignments))
# suppressPackageStartupMessages(library(rtracklayer))
# suppressPackageStartupMessages(library(mixtools))
# suppressPackageStartupMessages(library(mclust))
# suppressPackageStartupMessages(library(stringr))
# suppressPackageStartupMessages(library(patchwork))
# suppressPackageStartupMessages(library(LaplacesDemon)) #
# suppressPackageStartupMessages(library(mousetrap))
# suppressPackageStartupMessages(library(RColorBrewer))
# suppressPackageStartupMessages(library(ComplexHeatmap))
# suppressPackageStartupMessages(library(circlize))
# suppressPackageStartupMessages(library(ggvenn))
# suppressPackageStartupMessages(library(knitr))
# suppressPackageStartupMessages(library(janitor))
# suppressPackageStartupMessages(library(ggunchull))

# ensembl <- useMart('ensembl') mart <- useDataset('hsapiens_gene_ensembl',
# mart = ensembl)

anno_circlehunter2 <- function(ecDNA_bed_set = "") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    if (FALSE) {
        hg38_genes <- genes(txdb)
        hg38_genes_split <- split(hg38_genes, hg38_genes$gene_id)
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    hg38_genes <- genes(txdb, single.strand.genes.only = FALSE) %>%
        as.data.frame
    hg38_genes <- hg38_genes[!str_detect(hg38_genes$seqnames, "_"), ]
    hg38_genes <- hg38_genes %>%
        dplyr::select(group_name:strand) %>%
        dplyr::rename(gene_id = group_name)
    hg38_genes <- hg38_genes %>%
        filter(!gene_id %in% names(which((table(hg38_genes$gene_id) > 1))))
    rownames(hg38_genes) <- hg38_genes$gene_id
    hg38_genes <- GRanges(hg38_genes)
    hg38_genes_split <- split(hg38_genes, hg38_genes$gene_id)

    merge_bed <- c()

    for (i in c(1:length(ecDNA_bed_set))) {
        err <- try(read.table(ecDNA_bed_set[i]), silent = FALSE)
        if (!"try-error" %in% class(err)) {
            ecDNA_bed <- read.table(ecDNA_bed_set[i])
            colnames(ecDNA_bed) <- c("chrom", "start", "end", "name", "score", "strand",
                "start_ci", "end_ci", "start_peak", "end_peak", "start_cross", "end_cross",
                "linked_reads", "depth_mean", "high_coverage")
            ecDNA_bed$sample <- str_split(ecDNA_bed_set[i], pattern = "/", simplify = T)[,
                str_count(ecDNA_bed_set[i], pattern = "/") + 1]
            merge_bed <- rbind(merge_bed, ecDNA_bed)
        } else {
            print(paste0(ecDNA_bed_set[i], " has error output!"))
        }
    }

    merge_bed$label <- paste0(merge_bed$sample, "_", merge_bed$name)
    merge_bed$strand = "."
    merge_bed.raw <- merge_bed

    ecDNA_granges <- GRanges(merge_bed)

    In_ranges <- ecDNA_granges
    subjects <- hg38_genes_split

    Step1 <- findOverlaps(In_ranges, hg38_genes) %>%
        as_tibble
    Step2 <- as_tibble(In_ranges) %>%
        mutate(queryHits = as.integer(rownames(.)))
    Step1$ID <- Step2[match(Step1$queryHits, Step2$queryHits), "label"] %>%
        dplyr::pull(label)
    Step1$gene_id <- subjects[Step1$subjectHits, ] %>%
        as_tibble() %>%
        dplyr::pull(gene_id)
    OutPut_partial_geneid <- Step1 %>%
        dplyr::select(ID, gene_id)

    OutPut_partial_ngenes <- OutPut_partial_geneid %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(num_genes = length(unique(gene_id)))
    OutPut_partial_genesymbol <- OutPut_partial_geneid %>%
        as.data.frame()

    geneInfor <- getBM(attributes = c("entrezgene_id", "chromosome_name", "gene_biotype",
        "ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "band"),
        filters = "entrezgene_id", values = unique(as.character(OutPut_partial_geneid %>%
            pull(gene_id))), mart = mart)

    OutPut_partial_genesymbol$Symbol <- geneInfor[match(OutPut_partial_genesymbol$gene_id,
        geneInfor$entrezgene_id), "hgnc_symbol"]
    OutPut_partial_genesymbol <- na.omit(OutPut_partial_genesymbol)

    OutPut_partial_genesymbol <- OutPut_partial_genesymbol %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(num_genes = length(unique(gene_id)), genes = str_c(unique(Symbol),
            collapse = ",")) %>%
        as.data.frame()

    merge_bed$num_genes <- OutPut_partial_genesymbol[match(merge_bed$label, OutPut_partial_genesymbol$ID),
        "num_genes"]

    merge_bed$genes <- OutPut_partial_genesymbol[match(merge_bed$label, OutPut_partial_genesymbol$ID),
        "genes"]

    merge_bed$width <- merge_bed$end - merge_bed$start
    merge_bed$width_bp <- prettyNum(merge_bed$width, big.mark = ",", scientific = FALSE)

    merge_bed$ec_group <- paste0(merge_bed$sample, "|", str_extract(merge_bed$name,
        "ecDNA_(\\d)+_(\\d)+"))
    tmp <- merge_bed %>%
        group_by(ec_group) %>%
        dplyr::summarise(group_width = prettyNum(sum(width), big.mark = ",", scientific = FALSE),
            grp_wd = sum(width))
    merge_bed$group_width <- tmp[match(merge_bed$ec_group, tmp$ec_group), "group_width"] %>%
        dplyr::pull(group_width)
    merge_bed$grp_wd <- tmp[match(merge_bed$ec_group, tmp$ec_group), "grp_wd"] %>%
        dplyr::pull(grp_wd)

    return(merge_bed)
}

calculate_ecDNA_CPM <- function(fragment, qcfile, ecDNA_region, library = "scATAC") {
    # for scATAC-seq is__cell_barcode; for multi-scRNA+ATAC, is_cell indicated
    # if one cell have passed quality control.

    message(paste0(str_c(rep("#", 30), collapse = ""), " Perform calculate_ecDNA_CPM."))

    fragment <- data.table::fread(data.table = F, cmd = paste0("zcat ", fragment,
        " | grep -v '#' "))
    colnames(fragment) <- c("seqnames", "start", "end", "CB", "count")
    qcfile <- read.csv(qcfile)
    if (library == "scATAC") {
        pass_cells <- qcfile %>%
            dplyr::filter(is__cell_barcode == 1) %>%
            dplyr::pull(barcode) %>%
            sort
    } else if (library == "scMulti") {
        pass_cells <- qcfile %>%
            dplyr::filter(is_cell == 1) %>%
            dplyr::pull(barcode) %>%
            sort
    }
    fragment <- fragment[fragment$seqnames %in% paste0("chr", c(1:22)), ]
    fragment <- fragment[fragment$CB %in% pass_cells, ]
    AllReads <- fragment %>%
        dplyr::group_by(CB) %>%
        dplyr::summarise(AllReads = sum(count))
    RegionReads <- subsetByOverlaps(GRanges(fragment), ecDNA_region) %>%
        as.data.frame %>%
        dplyr::group_by(CB) %>%
        dplyr::summarise(RegionReads = sum(count))
    res <- dplyr::left_join(AllReads, y = RegionReads, by = "CB")
    res <- res %>%
        dplyr::mutate(Region_cpm = (RegionReads)/(AllReads) * 1e+06, Region_logcpm = log(Region_cpm +
            1))
    return(res)

}

run_GMM <- function(ecDNA_logCPM, keepCells = NA, cell_withsplit, runGMM = TRUE) {

    message(paste0(str_c(rep("#", 30), collapse = ""), " Perform run_GMM."))
    res <- ecDNA_logCPM

    # Remove cells without any reads within the region
    if (length(na.omit(keepCells)) > 0) {
        res.filtered <- res %>%
            dplyr::filter(CB %in% keepCells) %>%
            dplyr::filter(RegionReads > 0, !is.na(Region_logcpm)) %>%
            as.data.frame
        res.raw <- res %>%
            dplyr::filter(CB %in% keepCells) %>%
            as.data.frame
    } else {
        res.filtered <- res %>%
            dplyr::filter(RegionReads > 0, !is.na(Region_logcpm)) %>%
            as.data.frame
        res.raw <- res
    }

    # Run predict model
    GM_model <- normalmixEM(as.numeric(res.filtered[, "Region_logcpm"]), k = 2, epsilon = 1e-10,
        maxit = 1000, maxrestarts = 10)
    res.filtered <- cbind(res.filtered, GM_model$posterior)
    GM_model_mu <- GM_model$mu

    if (GM_model_mu[2] > GM_model_mu[1]) {
        res.filtered$predict <- ifelse(res.filtered$comp.2 > res.filtered$comp.1,
            "ecDNA+", "ecDNA-")
    } else {
        res.filtered$predict <- ifelse(res.filtered$comp.1 > res.filtered$comp.2,
            "ecDNA+", "ecDNA-")
    }
    # plot.mixEM(GM_model,whichplots=2,breaks = 100)

    res.merge <- left_join(res.raw, y = res.filtered[, c("CB", "comp.1", "comp.2",
        "predict")], by = "CB")
    res.merge$predict <- ifelse(is.na(res.merge$predict), "ecDNA-", res.merge$predict)
    rownames(res.merge) <- res.merge$CB

    # Plot result with ggplot2
    mu1 <- GM_model$mu[1]
    mu2 <- GM_model$mu[2]
    sigma1 <- GM_model$sigma[1]
    sigma2 <- GM_model$sigma[2]
    lambda1 <- GM_model$lambda[1]
    lambda2 <- GM_model$lambda[2]

    if (mu1 > mu2) {
        temp_mu <- mu1
        mu1 <- mu2
        mu2 <- temp_mu
        temp_sigma <- sigma1
        sigma1 <- sigma2
        sigma2 <- temp_sigma
        temp_lambda <- lambda1
        lambda1 <- lambda2
        lambda2 <- temp_lambda
    }

    res.merge$withsplit <- ifelse(res.merge$CB %in% cell_withsplit, "Discordant+",
        "Discordant-")

    ec.condition <- matrix(c(sum(res.merge$withsplit == "Discordant+" & res.merge$predict ==
        "ecDNA+"), sum(res.merge$withsplit == "Discordant+" & res.merge$predict !=
        "ecDNA+"), sum(res.merge$withsplit == "Discordant-" & res.merge$predict ==
        "ecDNA+"), sum(res.merge$withsplit == "Discordant-" & res.merge$predict !=
        "ecDNA+")), nrow = 2, dimnames = list(c("ecDNA+", "ecDNA-"), c("Discordant+",
        "Discordant-")))
    print(ec.condition)

    ec.condition %>%
        knitr::kable(booktabs = TRUE)
    fisher.test.p <- janitor::fisher.test(ec.condition)$p.value
    sprintf("P-value of fisher.test is %g.", round(fisher.test.p, 2)) %>%
        message()

    p1 <- res.merge %>%
        dplyr::filter(!is.na(Region_logcpm)) %>%
        ggplot(aes(x = Region_logcpm)) + geom_histogram(aes(y = ..density..), bins = 100,
        fill = "skyblue", color = "black") + stat_function(fun = function(x) dnorm(x,
        mean = mu1, sd = sigma1) * lambda1, color = "black") + stat_function(fun = function(x) dnorm(x,
        mean = mu2, sd = sigma2) * lambda2, color = "red") + theme_pubr() + scale_x_continuous(limits = c(min(res.merge$Region_logcpm),
        max(res.merge$Region_logcpm))) + theme(legend.position = "bottom", legend.title = element_blank(),
        axis.line = element_blank(), panel.border = element_rect(colour = "black",
            fill = NA, size = 1), axis.ticks = element_blank(), axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0), text = element_text(size = 12)) + labs(x = "",
        y = "")

    p2 <- ggplot() + theme_classic() + geom_jitter(data = res.merge[res.merge$withsplit ==
        "Discordant-", ], aes(x = Region_logcpm, y = -0.2, color = withsplit), alpha = 0.1,
        size = 2, height = 0.1) + geom_jitter(data = res.merge[res.merge$withsplit ==
        "Discordant+", ], aes(x = Region_logcpm, y = -0.5, color = withsplit), alpha = 0.6,
        size = 2, height = 0.1) + scale_y_continuous(breaks = c(-0.5, -0.1), labels = c("Discordant+",
        "Discordant-"), limits = c(-0.7, 0.1)) + scale_x_continuous(limits = c(min(res.merge$Region_logcpm),
        max(res.merge$Region_logcpm))) + scale_color_manual(values = c(`Discordant+` = "red",
        `Discordant-` = "black")) + theme(legend.position = "bottom", legend.title = element_blank(),
        text = element_text(size = 12), legend.text = element_text(size = 15), panel.border = element_rect(colour = "black",
            fill = NA, size = 0.5), axis.ticks.y = element_blank(), axis.text.y = element_text(color = c("red",
            "black"), size = 15), plot.margin = margin(0, 0, 0, 0)) + labs(y = "")

    dip_test.p <- round(dip.test(res.filtered$Region_logcpm)$p.value, 3)
    dip_test.text <- paste0("=", dip_test.p)

    if (runGMM) {
        p1 <- p1 + annotate("text", x = quantile(res.filtered$Region_logcpm, 0.5),
            y = 0.5, label = bquote(italic(P)[italic("dip.test")] ~ .(dip_test.text)),
            size = 6, color = "#f58231")
        # annotate('text',x=mu1,y=0.25,label = 'chrDNA',color='black',size=8)+
        # annotate('text',x=mu2,y=0.25,label = 'ecDNA',color='red',size=8)
    }

    ggHist <- p1/p2 + patchwork::plot_layout(heights = c(3, 1))

    mixtoolsClustering <- list(plot = ggHist, result = res.merge)
    return(mixtoolsClustering)
}

get_ecDNA_GeneScore <- function(ecDNA_region, genescoreMatrix) {

    message(paste0(str_c(rep("#", 30), collapse = ""), " Perform get_ecDNA_GeneScore."))

    ensembl <- biomaRt::useMart("ensembl", host = "dec2021.archive.ensembl.org")  # may not be connected due to network, should be solved after retrying several times or trying mirror.
    mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

    dt <- ecDNA_region %>%
        as.data.frame
    dt$name <- paste0("region", c(1:nrow(dt)))
    dt <- split(dt, dt$name)
    dt <- lapply(dt, function(x) {
        chr <- x[, "seqnames"] %>%
            gsub(pattern = "chr", replacement = "")
        start <- x[, "start"]
        end <- x[, "end"]
        biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position",
            "strand", "ensembl_gene_id", "entrezgene_id", "external_gene_name", "hgnc_symbol"),
            filters = c("chromosome_name", "start", "end"), values = list(chr, start,
                end), mart = mart)
    })
    dt <- do.call(rbind, dt)

    ec_genes <- dt %>%
        dplyr::filter(hgnc_symbol != "") %>%
        dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
        dplyr::pull(hgnc_symbol) %>%
        unique
    overlap_genes <- intersect(ec_genes, rownames(genescoreMatrix))

    message(sprintf("A total of %s genes carried by ecDNA overlap with GeneScoreMatrix.",
        length(overlap_genes)))

    ecGene.genescore <- genescoreMatrix[overlap_genes, ]
    ecGene.genescore <- log2(ecGene.genescore + 1)

    if (grepl(colnames(ecGene.genescore)[1], pattern = "#")) {
        colnames(ecGene.genescore) <- stringr::str_split(colnames(ecGene.genescore),
            pattern = "#", simplify = TRUE)[, 2]
    }

    return(ecGene.genescore)

}


run_kmean <- function(ecGene_genescore, cell_withsplit) {
    message(paste0(str_c(rep("#", 30), collapse = ""), " Perform run_kmean."))

    k2 <- kmeans(t(ecGene_genescore), 2)
    k2_df <- data.frame(k2_cluster = k2$cluster)
    k2_df$CB <- rownames(k2_df)
    k2_df$k2_cluster <- paste0("C", k2_df$k2_cluster)

    count_C1 <- sum(k2_df$k2_cluster == "C1")
    count_C2 <- sum(k2_df$k2_cluster == "C2")

    if (count_C1 > count_C2) {
        k2_df <- k2_df %>%
            mutate(k2_cluster = case_when(k2_cluster == "C1" ~ "C2", k2_cluster ==
                "C2" ~ "C1", TRUE ~ k2_cluster))
    }

    coordinate <- fpc::discrproj(t(ecGene_genescore), clvecd = k2$cluster, method = "dc",
        clnum, ignorepoints = FALSE, ignorenum = 0)$proj
    coordinate <- data.frame(coordinate)
    colnames(coordinate) <- paste0("dc", c(1:ncol(coordinate)))

    summarys <- coordinate
    summarys$k2_cluster <- k2_df[match(rownames(summarys), k2_df$CB), "k2_cluster"]
    summarys$withsplit <- ifelse(rownames(summarys) %in% cell_withsplit, "Discordant+",
        "Discordant-")

    cluster1.pct <- sum(summarys$withsplit == "Discordant+" & summarys$k2_cluster ==
        "C1")/sum(summarys$k2_cluster == "C1")
    cluster2.pct <- sum(summarys$withsplit == "Discordant+" & summarys$k2_cluster ==
        "C2")/sum(summarys$k2_cluster == "C2")
    ec.cluster <- ifelse(cluster2.pct > cluster1.pct, "C2", "C1")


    sprintf("%s (%d/%d) is ecDNA+ cluster; %d/%d carry discordant reads.", ec.cluster,
        sum(summarys$k2_cluster == ec.cluster), nrow(summarys), sum(summarys$withsplit ==
            "Discordant+" & summarys$k2_cluster == ec.cluster), sum(summarys$k2_cluster ==
            ec.cluster)) %>%
        message()

    ec.condition <- matrix(c(sum(summarys$withsplit == "Discordant+" & summarys$k2_cluster ==
        ec.cluster), sum(summarys$withsplit == "Discordant+" & summarys$k2_cluster !=
        ec.cluster), sum(summarys$withsplit == "Discordant-" & summarys$k2_cluster ==
        ec.cluster), sum(summarys$withsplit == "Discordant-" & summarys$k2_cluster !=
        ec.cluster)), nrow = 2, dimnames = list(c("ecDNA+", "ecDNA-"), c("Discordant+",
        "Discordant-")))
    print(ec.condition)

    ec.condition %>%
        knitr::kable(booktabs = TRUE)
    fisher.test.p <- janitor::fisher.test(ec.condition)$p.value
    sprintf("P-value of fisher.test is %g.", round(fisher.test.p, 2)) %>%
        message()

    summarys$aplhas <- ifelse(summarys$k2_cluster == ec.cluster, 1, 0.4)
    summarys$ecDNA <- ifelse(summarys$k2_cluster == ec.cluster, "ecDNA+", "ecDNA-")

    summarys$withsplit <- factor(summarys$withsplit, levels = c("Discordant+", "Discordant-"),
        ordered = TRUE)

    if (ec.cluster == "C1") {
        summarys$k2_cluster <- factor(summarys$k2_cluster, levels = c("C1", "C2"),
            ordered = TRUE)
    } else {
        summarys$k2_cluster <- factor(summarys$k2_cluster, levels = c("C2", "C1"),
            ordered = TRUE)
    }

    delta_r <- 0.02
    th_r <- 0.02
    delta <- diff(range(summarys$dc1)) * delta_r
    th <- diff(range(summarys$dc1)) * th_r

    ggplot(summarys, aes(x = dc1, y = dc2)) + geom_point(aes(alpha = aplhas, color = withsplit),
        size = 1.8) + theme_classic2() + scale_alpha_continuous(range = c(0.4, 1)) +
        scale_color_manual(values = c(`Discordant+` = "red", `Discordant-` = "black")) +
        guides(alpha = "none") + ggunchull::stat_unchull(aes(fill = k2_cluster),
        alpha = 0.2, th = th, delta = delta) + scale_fill_manual(values = c("red",
        "blue")) + facet_grid(ecDNA ~ .) + geom_point(data = summarys[summarys$withsplit ==
        "Discordant+", ], color = "red") + labs(x = "Coordinate_1", y = "Coordinate_2") +
        theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(size = 15),
            panel.border = element_rect(colour = "black", fill = NA, size = 1), legend.text = element_text(size = 20)) ->
            kmean.plot

    summarys$CB <- rownames(summarys)
    summarys$predict <- summarys$ecDNA
    summarys <- summarys %>%
        dplyr::select(CB, dc1, dc2, k2_cluster, withsplit, predict)

    kmeanClustering <- list(plot = kmean.plot, result = summarys)

    return(kmeanClustering)
}

run_All <- function(script_path, bam_path, circlehunter2_path, genescore_matrix = NA,
    fragment_path, qc.csv_path, out_path, tumorCells = NA, min_insertsize = 2000,
    ecDNA_id, sc_library = "scATAC", exportFigs = TRUE, useGMM = FALSE) {
    message(paste0(str_c(rep("#", 40), collapse = ""), " Run ecDNA identity assignment."))

    CPM_out <- paste0(out_path, "/result/", ecDNA_id, "_CPM.txt")
    ecDNA_out <- paste0(out_path, "/result/", ecDNA_id, "_determine.txt")

    if (!dir.exists(paste0(out_path, "/result/"))) {
        dir.create(paste0(out_path, "/result/"))
    }
    if (!dir.exists(paste0(out_path, "/result/png"))) {
        dir.create(paste0(out_path, "/result/png"))
    }

    circlehunter2.bed <- fread(circlehunter2_path, data.table = F)
    colnames(circlehunter2.bed)[which(grepl(colnames(circlehunter2.bed), pattern = "chrom"))] <- "seqnames"
    ecDNA.region <- circlehunter2.bed %>%
        dplyr::filter(grepl(name, pattern = ecDNA_id)) %>%
        GRanges

    cmd <- sprintf("mkdir -p %s && Rscript %s --circlehunter2_bed %s --cellranger_bam %s --outpath %s --insert_size %d",
        out_path, script_path, circlehunter2_path, bam_path, out_path, min_insertsize)
    cat(paste0(cmd, "\n"))
    system(cmd, intern = TRUE)

    runKmean <- ifelse(is.null(nrow(genescore_matrix)), FALSE, TRUE)

    if (!file.exists(CPM_out)) {
        ecDNA.logCPM <- calculate_ecDNA_CPM(fragment = fragment_path, library = sc_library,
            qcfile = qc.csv_path, ecDNA_region = ecDNA.region)

        write.table(ecDNA.logCPM, CPM_out, col.names = T, row.names = FALSE, quote = FALSE)
    } else {
        ecDNA.logCPM <- read.table(CPM_out, header = TRUE)
    }

    ecDNA.segments <- file.path(out_path, list.files(path = out_path, pattern = ecDNA_id))
    discordant.cell.txt <- ecDNA.segments[grepl(ecDNA.segments, pattern = "discordantcells.txt")]
    lapply(discordant.cell.txt, function(x) {
        readLines(x)
    }) %>%
        unlist %>%
        unique -> cell.withsplit

    if (length(na.omit(tumorCells)) == 0) {
        runGMM <- ifelse(dip.test(ecDNA.logCPM$Region_logcpm)$p.value < 0.05, TRUE,
            FALSE)
    } else {
        ecDNA.logCPM <- ecDNA.logCPM %>%
            dplyr::filter(CB %in% tumorCells)
        runGMM <- ifelse(dip.test(ecDNA.logCPM$Region_logcpm)$p.value < 0.05, TRUE,
            FALSE)
    }

    sprintf("Fitting Gaussian Mixture Model : %s.", runGMM)

    GMM.result <- run_GMM(ecDNA_logCPM = ecDNA.logCPM, keepCells = tumorCells, cell_withsplit = cell.withsplit,
        runGMM = runGMM)

    if (runKmean) {
        ec.genescore <- get_ecDNA_GeneScore(ecDNA_region = ecDNA.region, genescoreMatrix = genescore_matrix)
        kmean.result <- run_kmean(ecGene_genescore = ec.genescore, cell_withsplit = cell.withsplit)
        message(paste0(str_c(rep("#", 30), collapse = ""), " Plotting result."))
        grid.newpage()
        options(repr.plot.width = 16, repr.plot.height = 8)
        merge.plots <- GMM.result$plot | kmean.result$plot
        print(merge.plots)

    } else {
        message(paste0(str_c(rep("#", 30), collapse = ""), " Plotting result."))
        grid.newpage()
        options(repr.plot.width = 8, repr.plot.height = 8)
        print(GMM.result$plot)
        grid.newpage()

    }

    if (runKmean) {
        gmm.res <- GMM.result$result
        kmean.res <- kmean.result$result
        gmm.res <- gmm.res[intersect(rownames(kmean.res), rownames(gmm.res)), ]
        kmean.res <- kmean.res[intersect(rownames(kmean.res), rownames(gmm.res)),
            ]

        result.merge <- data.frame(row.names = gmm.res$CB, gmm.res[gmm.res$CB, c("CB",
            "AllReads", "RegionReads", "Region_cpm", "Region_logcpm", "comp.1", "comp.2")],
            gmm_predict = gmm.res[, "predict"], withsplit = gmm.res[, "withsplit"],
            k2_cluster = kmean.res[gmm.res$CB, "k2_cluster"], kmean_predict = kmean.res[gmm.res$CB,
                "predict"])

        if (useGMM | runGMM) {
            result.merge$final_decision <- result.merge$gmm_predict
        } else {
            result.merge$final_decision <- result.merge$kmean_predict
        }
        ecDNA.determine <- result.merge

        # Plot Venn plot.
        kmean.eccells <- ecDNA.determine %>%
            dplyr::filter(kmean_predict == "ecDNA+") %>%
            rownames
        gmm.eccells <- ecDNA.determine %>%
            dplyr::filter(gmm_predict == "ecDNA+") %>%
            rownames
        eccells.merge <- list(GMM = gmm.eccells, Kmean = kmean.eccells)
        venn.plot <- ggvenn::ggvenn(eccells.merge, show_elements = FALSE, show_percentage = FALSE,
            text_color = "red", stroke_color = "black", fill_color = c("#fc8d59",
                "#67a9cf"), text_size = 6, set_name_size = 8, fill_alpha = 0.4, stroke_alpha = 0.5)
        grid.newpage()
        options(repr.plot.width = 16, repr.plot.height = 6)
        print(venn.plot)

        if (TRUE) {
            ec.genescore <- ec.genescore[, intersect(colnames(ec.genescore), ecDNA.determine$CB)]
            ecDNA.determine <- dplyr::filter(ecDNA.determine, CB %in% intersect(colnames(ec.genescore),
                ecDNA.determine$CB))
        }

        mat <- ec.genescore %>%
            t() %>%
            as.matrix()
        mat <- mat[rownames(dplyr::arrange(ecDNA.determine, final_decision, withsplit)),
            ]

        withsplit.anno <- ecDNA.determine[rownames(mat), "withsplit"]
        names(withsplit.anno) <- ecDNA.determine$CB
        predict.anno <- ecDNA.determine[rownames(mat), "final_decision"]
        names(predict.anno) <- ecDNA.determine$CB

        col.h2 <- c("#00B050", "grey")
        names(col.h2) <- c("Discordant+", "Discordant-")
        col.h3 <- c("red", "grey")
        names(col.h3) <- c("ecDNA+", "ecDNA-")

        h2 <- Heatmap(rbind(Split_read = withsplit.anno), name = "discordant-reads",
            show_column_names = FALSE, col = col.h2, show_row_names = FALSE, heatmap_legend_param = list(title = "Discordant_read"))
        h3 <- Heatmap(rbind(ecDNA_predict = predict.anno), name = "ecDNA_predict",
            show_column_names = FALSE, heatmap_legend_param = list(title = "ecDNA_predict"),
            col = col.h3, show_row_names = FALSE)
        ha1 <- HeatmapAnnotation(ecDNA_logCPM = anno_points(ecDNA.determine[rownames(mat),
            "Region_logcpm"]), height = unit(2, "cm"))
        h1 <- Heatmap(mat = t(mat), show_column_names = FALSE, cluster_columns = FALSE,
            name = "h1", heatmap_legend_param = list(title = "logGeneScore)"), top_annotation = ha1)

        ht_list <- h3 %v% h2 %v% h1

        grid.newpage()
        options(repr.plot.width = 16, repr.plot.height = 2 + nrow(ec.genescore)/4)
        print(draw(ht_list))
        cut.point <- sum(ecDNA.determine$final_decision == "ecDNA+")/nrow(ecDNA.determine)
        decorate_annotation("ecDNA_logCPM", {
            grid.lines(c(cut.point, cut.point), c(0, 1), gp = gpar(lty = 2, lwd = 4,
                col = "red"))
        })
        grid.newpage()

        write.table(result.merge, ecDNA_out, col.names = T, row.names = FALSE, quote = FALSE)
        ecDNA.result <- list(plot = ht_list, result = result.merge)

    } else {
        if (runGMM) {
            ecDNA.result <- GMM.result
            ecDNA.result$result$final_decision <- ecDNA.result$result$predict
            write.table(ecDNA.result$result, ecDNA_out, col.names = T, row.names = FALSE,
                quote = FALSE)
        } else if (runKmean) {
            ecDNA.result <- kmean.result
            write.table(ecDNA.result$result, ecDNA_out, col.names = T, row.names = FALSE,
                quote = FALSE)
        } else if (useGMM) {
            GMM_out <- paste0(out_path, "/result/", ecDNA_id, "_poorGMM.txt")
            ecDNA.result <- GMM.result
            ecDNA.result$result$final_decision <- ecDNA.result$result$predict
            write.table(ecDNA.result$result, GMM_out, col.names = T, row.names = FALSE,
                quote = FALSE)
        } else {
            warning("GMM fitting performs poorly, set useGMM=TRUE to stick with fitting or provide genescore_matrix to run Kmean.")
        }
    }
    if (exportFigs == TRUE) {
        if (runGMM | useGMM) {
            while (!is.null(dev.list())) {
                dev.off()
            }
            png(filename = paste0(out_path, "/result/png/", ecDNA_id, "_gmm.png"),
                width = 2400, height = 1800, res = 300)
            print(GMM.result$plot)
            while (!is.null(dev.list())) {
                dev.off()
            }
        }
        if (runKmean) {
            while (!is.null(dev.list())) {
                dev.off()
            }
            png(filename = paste0(out_path, "/result/png/", ecDNA_id, "_kmean.png"),
                width = 2400, height = 1800, res = 300)
            print(kmean.result$plot)
            while (!is.null(dev.list())) {
                dev.off()
            }
            png(paste0(out_path, "/result/png/", ecDNA_id, "_heatmap.png"), width = 2400,
                height = 200 * (2 + nrow(ec.genescore)/4), res = 300)
            draw(ht_list)
            decorate_annotation("ecDNA_logCPM", {
                grid.lines(c(cut.point, cut.point), c(0, 1), gp = gpar(lty = 2, lwd = 4,
                  col = "red"))
            })
            while (!is.null(dev.list())) {
                dev.off()
            }
            grid.newpage()
        }
    }
    return(ecDNA.result)
}