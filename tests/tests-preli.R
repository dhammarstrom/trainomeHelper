## This script contains tests for the functions in seq_wrapper.R

library(trainomeHelper); library(tidyverse); library(glmmTMB); library(parallel)

# load sample data
data(rna_seq_sample); data(rna_seq_metadata)


# Change all counts to integer
rna_seq_sample_integer <- rna_seq_sample %>%
  mutate_at(2:197, ~ as.integer(round(., 0))) %>%
  print()
rna_seq_sample_integer[1, 2:197] <- NA
rna_seq_sample_integer[1,]

## Check the number of cores available
ncores <- parallel::detectCores()

## Save a list of arguments used by the fitting algorithm
args <- list(formula = y ~ time * condition + (1|participant),
             family = glmmTMB::nbinom2)





fits <- seq_wrapper(    fitting_fun = glmmTMB::glmmTMB,
                        arguments = args,
                        data = rna_seq_sample_integer,
                        metadata = rna_seq_metadata,
                        samplename = "seq_sample_id",
                        summary_fun = NULL,
                        eval_fun = NULL,
                        additional_vars = NULL,
                        exported = list(),
                        subset = 1:10,
                        cores = ncores)



rna_seq_sample_integer[1:10,] %>%
  pull(transcript_id)
