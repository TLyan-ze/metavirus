library(readr)
library(dplyr)
library(tidyr)
library(stringr)
Args <- commandArgs(trailingOnly = TRUE) 
blastp2rdrp.out <- Args[1]
proteinId2Family <- Args[2]

blastp2rdrp.info <- read_delim(blastp2rdrp.out, delim = '\t', col_names = FALSE)

colnames(blastp2rdrp.info) <- c("qseqid", "qlen", "sseqid", "slen", 
                                "qstart", "qend", "sstart", "send",
                                "evalue", "bitscore", "length", "pident", "mismatch",
                                "gaps", "stitle", "qcovhsp", "scovhsp")
##对每个蛋白的比对选bitscore最大的那个比对结果
blastp2rdrp.bestMatch.info <- data.frame()
for (pid in unique(blastp2rdrp.info$qseqid)) {
    tmp <- blastp2rdrp.info %>% filter(qseqid == pid)
    tmp.use <- tmp[which.max(tmp$bitscore),]
    blastp2rdrp.bestMatch.info <- rbind(blastp2rdrp.bestMatch.info, tmp.use)
}
##获取contigId-proteinId-Family对应信息
blastp2rdrp.bestMatch.use <- blastp2rdrp.bestMatch.info %>%
    separate(stitle, into = c("accession", "rdrpLength", "Species", "Family"), sep = "\\|") %>%
    separate(qseqid,into = c("kmer","contigIndex","proteinIndex"), sep = "_") %>%
    mutate(
        contigId = str_c(kmer, contigIndex, sep = "_"),
        proteinId = str_c(kmer, contigIndex, proteinIndex, sep = "_"),
    ) %>%
    select(contigId, proteinId, Family, bitscore)
##若一个contigs上预测出多个蛋白,则保留bitscore最大的那个
contigs <- unique(blastp2rdrp.bestMatch.use$contigId)  
blastp2rdrp.bestMatch.out <- data.frame()
for (cid in contigs) {
    tmp <- blastp2rdrp.bestMatch.use[which(blastp2rdrp.bestMatch.use$contigId == cid),]
    tmp.use <- tmp[which.max(tmp$bitscore),]
    blastp2rdrp.bestMatch.out <- rbind(blastp2rdrp.bestMatch.out, tmp.use)
}
blastp2rdrp.bestMatch.out <- blastp2rdrp.bestMatch.out %>% select(-bitscore)
write_delim(blastp2rdrp.bestMatch.out, proteinId2Family, delim = '\t', quote_escape = FALSE)

