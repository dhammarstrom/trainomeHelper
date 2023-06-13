
library(tidyverse)
### Load transcripts

raw_dat <- read_csv("ome-data/vol_transcript_rsem.csv")

data("vol_samples")

metadata <- vol_samples


tiny_dat <-raw_dat[1:32,]

rownames(t(tiny_dat[1, 2:197]))


count_dfs <- split(tiny_dat[, 2:197], tiny_dat$transcript_id)

x <- count_dfs[[1]]

trans_countdfs <- function(x, metadata) {

  df <- merge(data.frame(seq_sample_id = rownames(t(x)), counts = t(x)),
        data.frame(metadata), by = "seq_sample_id")

  return(df)

}

counts_designs <- lapply(count_dfs, trans_countdfs, metadata = metadata)

counts_designs[[1]]
