#' A function to fit models with a chosen fitting algorithm. To be used in seq_wrapper.
#' @param FUN Name of a fitting function, like glmmTMB::glmmTMB
#' @param ARG A lsit of arguments that can be evaluated by the fitting function
fit_fun <- function(FUN, ARG) {

  # Fit the model using the selected machinery and available arguments

  model <- tryCatch(do.call(FUN, ARG),
           error = function(e) {
             message("Error in fitting model: ", e)
             return(NULL)
           })

  # Remove abundant info in the call
  # model$call <- ARG$formula

  # return the model
  return(model)

}


#' Transform, merge and fit models. The function is used inside seq_wrapper to combine
#' metadata with target-level data and perform the model fitting. The function is used in a
#' call to pbapply::pblapply.
#' @param x A list of target-specific counts
#' @param samp_name Sample names from the upper level function
#' @param arg_list Arguments from the upper level function
#' @param ffun the fitting function from the upper level function
transform_merge_fit <- function(x,
                                samp_name = samplename,
                                metdat = metadata,
                                arg_list = arguments,
                                add_vars = additional_vars,
                                ffun = fitting_fun) {


  transposed <- data.frame(seq_sample_id = rownames(t(x[,-1])), y = as.numeric(t(x[,-1])))

  colnames(transposed)[1] <- samp_name

  df <- merge(transposed, data.frame(metdat), by = samp_name)

  ## Keep only data needed for fitting
  if("formula" %in% names(arg_list))  parsed <- all.vars(as.formula(arg_list$formula))
  if("model" %in% names(arg_list))  parsed <- all.vars(as.formula(arg_list$model))
  if(all(c("fixed", "random") %in% names(arg_list)))  parsed <- c(all.vars(as.formula(arg_list$fixed)), all.vars((arg_list$random[[1:length(arg_list$random)]])))

  ## Keep also additional variables that exists in the meta data data set
  if(!is.null(add_vars)) parsed <- c(parsed, add_vars)


  df <- df[,parsed, drop = FALSE]


  ## Remove attributes from the list of arguments (this solves an issue when using glmmTMB)


  arguments_final <- append(arguments, list(data = df))

for(i in 1:length(arguments_final)) {

  if(class(arguments_final[[i]]) == "formula") environment(arguments_final[[i]]) <- NULL

}


  model <- fit_fun(ffun, arguments_final)

  return(model)

}



#' A flexible upper-level wrapper for iterative modelling using any available fitting algorithm.
#'
#' @param fiting_fun A model fitting function like stats::lm, glmmTMB::glmmTMB, lme4::lmer
#' @param arguments A list of arguments to be passed to the fitting function, this should not contain data. Note that the formula must have y as the dependent variable.
#' @param data A data frame with targets (i.e. genes, transcripts) as rows and sample names as colums.
#' The first column is assumed to contain target names/identification
#' @param metadata A data frame with sample names (corresponding to column names in the target matrix)
#' and design variables.
#' @param samplename A character value indicating the variable by which metadata can merge with the target data.
#' This defaults to "seq_sample_id" as this is used in the trainomeMetaData package.
#' @param additional_vars A vector of additional variables that is contained in the metadata data set that is needed to fit the model. By default the metadata is reduced to variables contained in the slots formula/model/fixed/random in additional arguments. More variables may be needed for offsets, weights etc.
#' @param summary_fun A custom (user-created) function for evaluating/summarizing models. If NULL, a list of fitted models are returned
#' @param eval_fun A custom (user-created) function for model diagnostics/evaluation. If NULL, no evaluation/diagnostics of models are returned
#' @param exported A list of functions, values etc. to be passed to summary_fun and eval_fun. This list must contain any functions that should be used in model summarise or evaluations.
#' @param subset A sequence, random samples or integers to indicate which rows to keep in data. This is useful if you want to test the model in a subset of targets. If keft to the default (NULL), all rows will be used.
#' @param cores An integer indicating the number of cores to be used in parallel computations. If NULL, a sequential for loop is used. If "max", all available cores are used.
#' @return A nested list with three upper levels slots: models, a list of fitted objects; summaries, a list of summaries created from the summary_fun function; evaluations, a list of diagnostics created from eval_fun.
#' @details This function provides a flexible wrapper to fit, summarize and evaluate statistical models fitted to high dimensional omics-type data.
#' @export
seq_wrapper <- function(fitting_fun = glmmTMB::glmmTMB,
                        arguments,
                        data,
                        metadata,
                        samplename = "seq_sample_id",
                        additional_vars = NULL,
                        summary_fun = NULL,
                        eval_fun = NULL,
                        exported = list(),
                        subset = NULL,
                        cores) {


  ## Prepare data #####
  data <- data.frame(data)
  metadata <- data.frame(metadata)

  ## Subset
  if(!is.null(subset)) data <- data[subset,]

  ## Sanity checks

  # checks if arguments for the provided fitting function matches arguments
  stopifnot("arguments do not match named argumenst of the selected \nmodel fitting function (fitting_fun)." = names(arguments) %in% names(formals(glmmTMB::glmmTMB)) )

  # Check if the data has a character vector for first column (indicating transcript id).
  stopifnot("The first column of the data is not character or factor, check if this column indicate target identifications." = is.character(data[,1]))

  # Check if the sample name is present in the meta data
  stopifnot("The samplename does not exist in the metadata,\nno variable for matching metadata and data." = samplename %in% colnames(metadata))







  ## Count tables from the trainomeMetaData package comes with the first column
  # being trainscript id's. This splits the data into data frames based on number of
  # targets
  count_dfs <- split(data[,1:ncol(data)], data[,1])



  # To merge meta data with meta data and fit the model, return model object

  ### Fitting models in parallel #######

  # Determine the number of cores
  if(is.null(cores)) num_cores <- 1
  if(cores >= parallel::detectCores()) num_cores <- parallel::detectCores()
  if(cores <= parallel::detectCores()) num_cores <- cores


  # Create a cluster using the number of cores specified
  cl <- parallel::makeCluster(num_cores)

  ## Export data to clusters
  parallel::clusterExport(cl, varlist = c("metadata",
                                "arguments",
                                "fit_fun",
                                "samplename",
                                "additional_vars",
                                "fitting_fun"),
                          envir = environment())



  # Parallel execution of the fitting process
  cat("Transforming, merging and modelling data.\n")
  model_fits <- pbapply::pblapply(cl = cl, X = count_dfs, FUN = trainomeHelper:::transform_merge_fit)


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
    parallel::clusterExport(cl2, varlist = c("exported"), envir = environment())

    # Using the user provided function

    results$model_summarises <- pbapply::pblapply(cl = cl2, X = model_fits, FUN = summary_fun)


    parallel::stopCluster(cl2)

  }

  if(!is.null(eval_fun)){

    cat("Evaluating models\n")

    # Create a cluster using the number of cores specified
    cl3 <- parallel::makeCluster(num_cores)

    ## Export data to clusters
    parallel::clusterExport(cl3, varlist = c("exported"), envir = environment())

    # Using user-provided function
    results$model_evaluations <- pbapply::pblapply(cl = cl3, X = model_fits, FUN = eval_fun)


    parallel::stopCluster(cl3)

  }


  return(results)



}


