########### This script is to visualize the results of fisher's test ###########

library(Gviz)
library(GenomicRanges)
library(rtracklayer)

############ label tss 
axisTrack <- GenomeAxisTrack(IRanges(start = 451197, end = 452197,names = 'TSS'),
                             chromosome = 'chr1')
plotTracks(axisTrack, from = 450000, to = 460000, showId=TRUE)
###############################################################

############ plot gene region 
aTrack <- AnnotationTrack(start = 450703, end = 451697, chromosome = 'chr1',
                          strand = '-', name = 'ENSG00000284733')
plotTracks(list(axisTrack, aTrack), from = 450000, to = 460000, showId=TRUE)
##########################################################

############  plot bw data 
hs_cd14_plus = import.bw('C:/Users/lenovo/Desktop/pausing_scripts/proseq/PROseq-HUMAN-CD14-1-1_dedup_QC_end_plus.bw')

hs_cd14_minus = import.bw('C:/Users/lenovo/Desktop/pausing_scripts/proseq/PROseq-HUMAN-CD14-1-1_dedup_QC_end_minus.bw')
options(ucscChromosomeNames=FALSE)
dtrack_plus = DataTrack(range = hs_cd14_plus, start = 450703, end = 451697, 
                         name = "hs_cd14_plus", chromosome = "1", type="h", col = 'red')

plotTracks(list(axisTrack, aTrack, dtrack_plus), from = 450000, to = 460000)
