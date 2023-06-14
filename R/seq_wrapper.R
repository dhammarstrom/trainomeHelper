

dat <- data.frame(time = rep(c(0,1), each = 10),
                   count = rnbinom(20, size = 20, mu = 10),




                                    id = rep(seq(1:10)))



#' A function to fit models with a chosen fitting algorithm. To be used in seq_wrapper.
#' @param FUN Name of a fitting function, like glmmTMB::glmmTMB
#' @param ARG A lsit of arguments that can be evaluated by the fitting function
fit_fun <- function(FUN, ARG) {

  # Fit the model using the selected machinery and available arguments
  model <- do.call(FUN, ARG)

  # Remove abundant info in the call
  # model$call <- ARG$formula

  # return the model
  return(model)

}






#' A flexible upper-level wrapper for iterative modelling using any available fitting algorithm.
#'
#' @param fiting_fun A model fitting function like stats::lm, glmmTMB::glmmTMB, lme4::lmer
#' @param arguments A list of arguments to be passed to the fitting function, this should not contain data.
#' @param data A data frame with targets (i.e. genes, transcripts) as rows and sample names as colums.
#' The first column is assumed to contain target names/identification
#' @param metadata A data frame with sample names (corresponding to column names in the target matrix)
#' and design variables.
#' @param samplename A character value indicating the variable by which metadata can merge with the target data.
#' This defaults to "seq_sample_id" as this is used in the trainomeMetaData package.
#' @param summary_fun A custom (user-created) function for evaluating/summarizing models. If NULL, a list of fitted models are returned
#' @param eval_fun A custom (user-created) function for model diagnostics/evaluation. If NULL, no evaluation/diagnostics of models are returned
#' @param exported A list of functions, values etc. to be passed to summary_fun and eval_fun. This list must contain any functions that should be used in model summarise or evaluations.
#' @param cores An integer indicating the number of cores to be used in parallel computations. If NULL, a sequential for loop is used. If "max", all available cores are used
#' @return A nested list with three upper levels slots: models, a list of fitted objects; summaries, a list of summaries created from the summary_fun function; diagnostics, a list of diagnostics created from eval_fun
#' @details
#' @export
#'
seq_wrapper <- function(fitting_fun = glmmTMB::glmmTMB,
                        arguments,
                        data,
                        metadata,
                        samplename = "seq_sample_id",
                        summary_fun = NULL,
                        eval_fun = NULL,
                        exported = list(),
                        cores) {



  ## Sanity checks

  # checks if arguments for the provided fitting function matches arguments
  stopifnot("arguments do not match named argumenst of the selected \nmodel fitting function (fitting_fun)." = names(arguments) %in% names(formals(glmmTMB::glmmTMB)) )

  # Check if the data has a character vector for first column (indicating transcript id).
  stopifnot("The first column of the data is not character or factor, check if this column indicate target identifications." = is.character(data[,1]))

  # Check if the sample name is present in the meta data
  stopifnot("The samplename does not exist in the metadata,\nno variable for matching metadata and data." = samplename %in% colnames(metadata))



  ## Prepare data #####

  ## Count tables from the trainomeMetaData package comes with the first column
  # being trainscript id's. This splits the data into data frames based on number of
  # targets
  count_dfs <- split(data[,1:ncol(data)], data[,1])



  # To merge meta data with meta data and fit the model, return model object
  transform_merge_fit <- function(x,
                                  samp_name = samplename,
                                  metdat = metadata,
                                  arg_list = arguments,
                                  ffun = fitting_fun) {

    transposed <- data.frame(seq_sample_id = rownames(t(x)), counts = as.numeric(t(x)))

    colnames(transposed)[1] <- samp_name

    df <- merge(transposed, data.frame(metdat), by = samp_name)

    ## Keep only data needed for fitting
    if("formula" %in% names(arg_list))  parsed <- all.vars(as.formula(arg_list$formula))
    if("model" %in% names(arg_list))  parsed <- all.vars(as.formula(arg_list$model))
    if(all(c("fixed", "random") %in% names(arg_list)))  parsed <- c(all.vars(as.formula(arg_list$fixed)), all.vars((arg_list$random[[1:length(arg_list$random)]])))


    df <- df[,parsed, drop = FALSE]


    ## Remove attributes from the list of arguments
    one_entry <- function(x) {
      for (i in length(x)) attr(x[[i]], "names") <- NULL
      environment(x) <- NULL
      return(x)
    }

    arguments_final <- append(arg_list, list(data = df))
    arguments_final_noattr <- lapply(arguments_final, one_entry)



    model <- fit_fun(ffun, arguments_final_noattr)

    return(model)

  }

  ### Fitting models in parallel #######

  # Determine the number of cores
  if(is.null(cores)) num_cores <- 1
  if(cores >= parallel::detectCores()) num_cores <- parallel::detectCores()
  if(cores < parallel::detectCores()) num_cores <- cores


  # Create a cluster using the number of cores specified
  cl <- parallel::makeCluster(num_cores)

  ## Export data to clusters
  parallel::clusterExport(cl, c("metadata",
                                "arguments",
                                "fit_fun",
                                "samplename",
                                "fitting_fun"))



  # Parallel execution of the fitting process
  cat("Transforming, merging and modelling data.\n")
  model_fits <- pbapply::pblapply(cl = cl, X = count_dfs, FUN = transform_merge_fit)


  parallel::stopCluster(cl)


  ## Store results in a list
  results <- list(model_fits = model_fits,
                  model_summarises = NULL,
                  model_evaluations = NULL)

  # Summary and evaluation functions ###############

  ## Summary functions
  if(!is.null(summary_fun)) {

    cat("Summarizing models\n")
    # Create a cluster using the number of cores specified
    cl2 <- parallel::makeCluster(num_cores)

    ## Export data to clusters
    parallel::clusterExport(cl2, c("exported"))

    # Using the user provided function

    results$model_summarises <- pbapply::pblapply(cl = cl2, X = model_fits, FUN = summary_fun)


    parallel::stopCluster(cl2)

  }

  if(!is.null(eval_fun)){

    cat("Evaluating models\n")

    # Create a cluster using the number of cores specified
    cl3 <- parallel::makeCluster(num_cores)

    ## Export data to clusters
    parallel::clusterExport(cl3, c("exported"))

    # Using user-provided function
    results$model_evaluations <- pbapply::pblapply(cl = cl3, X = model_fits, FUN = eval_fun)


    parallel::stopCluster(cl3)

  }


  return(results)



}



custom_sum <- function(x) {

  res <- data.frame(marginaleffects::avg_comparisons(x,
                                   variables = list(time = c("w0", "w2pre"))))



  return(res)

}








##
tdf <- tiny_dat %>%
  mutate_if(is.numeric, round) %>%
  mutate_if(is.numeric, as.integer) %>%
  print()


glmmTMB::glmmTMB(formula = counts ~ time, family = glmmTMB)



fits <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                        arguments = list(formula =counts ~ time + (1|participant), family = glmmTMB::nbinom2(link = "log")),
                        data = data.frame(tdf),
                        metadata = data.frame(metadata),
                        samplename = "seq_sample_id",
                    summary_fun = custom_sum,
                    exported = list(marginaleffects::avg_comparisons),
                        cores = 4)


bind_rows(fits$model_summarise) %>%
  mutate(transcript = names(fits$model_summarises)) %>%
  ggplot(aes(p.value)) + geom_histogram()



fits$model_summarises
