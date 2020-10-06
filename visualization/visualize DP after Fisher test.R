###############################################################################
#### This script is to visualize fifferential pausing with Gviz ###############
###############################################################################
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(tuSelecter2)
library(dplyr)

fisher_path = 'C:/Users/lenovo/Desktop/pausing_scripts/output/fisher/'
testing_results = read.csv(paste0(fisher_path, 'human_fisher_results.tsv'), header = T, sep = '\t')
testing_0.01 = testing_results[which(testing_results$q_value < 0.01), ] #number 2461

#add chr, start, end, to df
human_mart = readRDS('C:/Users/lenovo/Desktop/pausing_scripts/human_mart.RDS')[,c(1:5)]
human_mart = human_mart[!duplicated(human_mart),]
names(human_mart) = c('chrom', 'start', 'end', 'strand', 'gene_id')
testing_0.01 = inner_join(testing_0.01, human_mart, by = 'gene_id')

#plot
cd4_minus_bwfile = 'C:/Users/lenovo/Desktop/pausing_scripts/proseq/PROseq-HUMAN-CD4-1-1_dedup_QC_end_minus.bw'
cd4_minus <- import(con = cd4_minus_bwfile, format = 'bigWig')
cd4_plus_bwfile = 'C:/Users/lenovo/Desktop/pausing_scripts/proseq/PROseq-HUMAN-CD4-1-1_dedup_QC_end_plus.bw'
cd4_plus <- import(con = cd4_plus_bwfile, format = 'bigWig')

cd14_minus_bwfile = 'C:/Users/lenovo/Desktop/pausing_scripts/proseq/PROseq-HUMAN-CD14-1-1_dedup_QC_end_minus.bw'
cd14_minus <- import(con = cd14_minus_bwfile, format = 'bigWig')
cd14_plus_bwfile = 'C:/Users/lenovo/Desktop/pausing_scripts/proseq/PROseq-HUMAN-CD14-1-1_dedup_QC_end_plus.bw'
cd14_plus <- import(con = cd14_plus_bwfile, format = 'bigWig')

options(ucscChromosomeNames=FALSE)
dTrack1 <- DataTrack(cd4_plus, name = "cd4_plus",chromosome = '4', type = 'h', ylim = )
dTrack2 <- DataTrack(cd14_plus, name = "cd14_plus",chromosome = '4', type = 'h')
plotTracks(list(dTrack1, dTrack2), from = 17577192, to = 17607972)

cd14_tq = readRDS('C:/Users/lenovo/Desktop/pausing_scripts/proseq/tq/PROseq-HUMAN-CD14-1-1.RDS')
cd4_tq = readRDS('C:/Users/lenovo/Desktop/pausing_scripts/proseq/tq/PROseq-HUMAN-CD4-1-1.RDS')
plot_model(cd14_tq,'ENSG00000002549')
plot_model(cd4_tq,'ENSG00000002549')


#add annotation track
human_tx = readRDS('C:/Users/lenovo/Desktop/pausing_scripts/hsapiens_transcript_grng.RDS')
sel_chr <- sort(str_subset(unique(seqnames(human_tx)), "^[0-9]+.?$|^X"))
human_tx <- human_tx[(seqnames(human_tx) %in% c('1'))]

grtrack <- GeneRegionTrack(human_tx, name = "Human Tx", start=16252, end=108262,
                           transcriptAnnotation = 'ensembl_transcript_id',
                           group = 'ensembl_gene_id',
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")
plotTracks(list(dTrack1, dTrack2, grtrack), from = 16252, to = 108262)
