########### This script is to visualize the results of fisher's test ###########
library(Gviz)
library(GenomicRanges)
library(rtracklayer)

# change this dir to point where your files locate
bw_dir <- "C:/Users/lenovo/Desktop/pausing_scripts/proseq"
# bw_dir <- "/media/yixin/Elements1/data/processed/comparaReg/proseq/PROseq-HUMAN-CD14-1-1"

# arguments for gene
gene_name <- "ENSG00000175538"
chrom <- 11
gene_start <- 74454841
gene_end <- 74467729
strand <- "-"
tss_expand <- 250

# arguments for plotting
expand_region <- 1000
plot_start <- gene_start - expand_region
plot_end <- gene_end + expand_region

options(ucscChromosomeNames=FALSE)

############ label tss
axisTrack <-
    GenomeAxisTrack(IRanges(start = plot_start, end = plot_end),
                    chromosome = chrom)

# plotTracks(axisTrack, showId=TRUE)
make_tss_track <- function(strand) {
    if (strand == "+") {
        tss_start = gene_start
        tss_end = gene_start + tss_expand
    } else if (strand == "-") {
        tss_start = gene_end - tss_expand
        tss_end = gene_end
    } else {
        print("strand is neither + nor -")
    }
    tss_track <- AnnotationTrack(
        start = tss_start,
        end = tss_end,
        name = "TSS",
        shape = "box",
        chromosome = chrom)
    return(tss_track)
}

tss_track <- make_tss_track(strand)

############ plot gene region
aTrack <- AnnotationTrack(start = gene_start, end = gene_end, chromosome = chrom,
                          strand = strand, name = gene_name)
# plotTracks(list(axisTrack, aTrack), from = plot_start, to = plot_end, showId=TRUE)

############  plot bw data
import_bw <- function(bw_name) {
    bw_rng <- import.bw(file.path(bw_dir, bw_name),
              which = GenomicRanges::GRanges(
                  seqnames = chrom,
                  ranges = IRanges::IRanges(start = plot_start,
                                            end = plot_end)))
    return(bw_rng)
}

hs_cd14_plus <- import_bw('PROseq-HUMAN-CD14-1-1_dedup_QC_end_plus.bw')
hs_cd14_minus <- import_bw('PROseq-HUMAN-CD14-1-1_dedup_QC_end_minus.bw')

dtrack_plus <- DataTrack(range = hs_cd14_plus,
                         name = "hs_cd14_plus", chromosome = chrom,
                         window = -1, windowSize = 250,
                         type = "h", col = 'blue', strand = "+")

dtrack_minus <- DataTrack(range = hs_cd14_minus,
                         name = "hs_cd14_plus", chromosome = chrom,
                         window = -1, windowSize = 250,
                         type = "h", col = 'red', strand = "-")

plotTracks(list(axisTrack, aTrack, tss_track, dtrack_plus, dtrack_minus),
           from = plot_start, to = plot_end, chromosome = chrom)
