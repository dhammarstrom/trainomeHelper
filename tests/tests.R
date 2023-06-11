### Tests ####


library(tidyverse);library(trainomeMetaData)
library(glmmTMB)

# download_ome(download = "vol_transcript_rsem")



volrnaseq <- read_csv("ome-data/vol_transcript_rsem.csv")

colnames(volrnaseq)

transcript <- volrnaseq[1,] %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "counts",
               cols = FP11w0L:FP9w2preR) %>%
  mutate(counts = as.integer(round(counts,0))) %>%
  print()



test_data <- vol_samples %>%
  inner_join(transcript)



model_settings <- list(formula = list(counts ~  time + time:condition + (1|participant),
                                      counts ~  time + time:condition + (1|participant)),
                       family = list(glmmTMB::nbinom2(), glmmTMB::genpois()))







test_results <- singlefit(model_settings = model_settings, data = test_data)

test_results



