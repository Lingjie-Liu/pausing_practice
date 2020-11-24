#! /usr/bin/env Rscript

####
## This script is to prepare the gene bodies bedfiles with windows
####

# library(argparse)
library(GenomicRanges)
library(rtracklayer)
# library(dplyr)

# parser <- ArgumentParser()
# parser$add_argument("-a","--annotation", dest="annotation", type = "character",
# 					                    required=TRUE, help="bed file containing of non overlapping genes")
# parser$add_argument("-w","--windowsize", dest="windowsize", type = "integer",
# 					                    required=TRUE, help="window size for to seperate the gene bodies")
# #parser$add_argument("-o","--output", dest="output", type="character",
# #					                    default = "stdout", help="File to store output in. Defaults to stdout.")
# args <- parser$parse_args()
#

# a couple of parameters
# window size
wd_size <- 100
# number of bins
bin_n <- 100
total_len <- wd_size * bin_n

# bed_file <- 'C:/Users/lenovo/Desktop/pausing_scripts/revised_sytenic/test_non_overlapping_genes.bed'
bed_file <- '/home/yixin/Desktop/github_repo/pausing_practice/data/test_non_overlapping_genes.bed'

# convert bed file into granges object
# Note your test file only has genes from positive strand, so I give you an example about positive strand
gene_gr <- import(bed_file)
gene_gr$score <- NULL

# filtered out genes with less than total_len
gene_gr <- gene_gr[width(gene_gr) > total_len]
# width(gene_gr)
# fixed length for every gene
gene_gr <- promoters(gene_gr, upstream = 0, downstream = total_len)
gene_gr_ls <- tile(gene_gr, n = bin_n)
names(gene_gr_ls) <- gene_gr$name

gene_gr_wd <- unlist(gene_gr_ls)
# add the gene name back to a column of the metadata if you prefer
gene_gr_wd$name <- names(gene_gr_wd)

gene_gr_wd

# # function for producing all windows for a gene
# produce_windows <- function(gr, ws){
#   if (as.character(strand(gr)) == '+'){
#     start = start(gr) + 1000
#     end = start(gr) + 10000
#     start_array = seq(start, end-ws, ws)
#     end_array = seq(start+ws-1, end, ws)
#   }
#   else if (as.character(strand(gr)) == '-'){
#     start = end(gr) - 10000
#     end = end(gr) - 1000
#     start_array = seq(start, end-ws, ws)
#     end_array = seq(start+ws-1, end, ws)
#   }
#   sub_gr = GRanges(seqnames = seqnames(gr),
#                    ranges = IRanges(start = start_array, end = end_array),
#                    strand = strand(gr),
#                    name = gr$name)
#   return(sub_gr)
# }
#
#
# result = produce_windows(gene_gr[1], 50)
# result


