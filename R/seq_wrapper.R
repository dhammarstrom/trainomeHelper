

dat <- data.frame(time = rep(c(0,1), each = 10),
                   count = rnbinom(20, size = 20, mu = 10),
                  id = rep(seq(1:10)))


#' A flexible upper-level wrapper for iterative modelling using any available fitting algorithm.
#'
#' @param fiting_fun A model fitting function like stats::lm, glmmTMB::glmmTMB, lme4::lmer
#' @param arguments A list of argumnets to be passed to the fitting function, this should not contain data.
#' @param data A data frame with targets (i.e. genes, transcripts) as rows and sample names as colums.
#' The first column is assumed to contain target names/identification
#' @param metadata A data frame with sample names (corresponding to column names in the target matrix)
#' and design variables.
#' @param samplename A character value indicating the variable by which metadata can merge with the target data.
#' This defaults to "seq_sample_id" as this is used in the trainomeMetaData package.
#' @param summary_fun A custom (user-created) function for evaluating/summarizing models. If NULL, a list of fitted models are returned
#' @param eval_fun A custom (user-created) function for model diagnostics/evaluation. If NULL, no evaluation/diagnostics of models are returned
#' @param cores An integer indicating the number of cores to be used in parallel computations. If NULL, a sequential for loop is used. If "max", all availlable cores are used
#' @return A nested list with three upper levels slots: models, a list of fitted objects; summaries, a list of summaries created from the summary_fun function; diagnostics, a list of diagnostics created from eval_fun
#' @details
#' @export



seq_wrapper <- function(fitting_fun = glmmTMB::glmmTMB,
                        arguments,
                        data,
                        metadata,
                        cores) {

  ## Prepare data #####

  ## Count tables from the trainomeMetaData package comes with the first column
  # being trainscript id's. This splits the data into data frames based on number of
  # targets
  count_dfs <- split(data[,1:ncol(data)], data[,1])

  x <- count_dfs[[1]]

  # To merge meta data with meta data
  trans_countdfs <- function(x, metadata) {

    df <- merge(data.frame(seq_sample_id = rownames(t(x)), counts = t(x)),
                data.frame(metadata), by = samplename)

    return(df)

  }

  # This results in a list of data sets, each named by the target id.
  counts_designs <- lapply(count_dfs, trans_countdfs, metadata = metadata)




  ### Fitting models in parallel #######

  # Determine the number of cores
  if(is.null(cores)) num_cores <- parallel::detectCores()
  if(cores >= parallel::detectCores()) num_cores <- parallel::detectCores()
  if(cores < parallel::detectCores()) num_cores <- cores


  # Create a cluster using the number of cores specified
  cl <- makeCluster(num_cores)

  # Create the fitting function
  fit_fun <- function(fitting_fun, arguments) {

    # Fit the model using the selected machinery and available arguments
    model <- do.call(fitting_fun, arguments)
    # Return the model
    return(model)

  }


  arguments <- list(formula = counts~time + (1|participant))

  # Combine data to arguments
  args2 <- lapply(counts_designs, function(x) append(list(data = x), arguments))

  args2[[1]]

  # Parallel
  parLapply(cl, inputs, fun)


}






m1 <- seq_wrapper(data = dat, formula = count ~ time + (1|id))

summary(m1)




