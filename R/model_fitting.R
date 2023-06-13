
library(tibble)

model_settings <- list(formula = list(counts ~  time + time:condition + (1|participant),
                                      counts ~  time + time:condition + (1|participant)),
                       family = list("nbinom2", "nbinom"),
                       ziformula = list(NULL, NULL),
                       dispformula = list(~1, ~1),
                       weights = list(NULL, NULL),
                       offset = list(NULL, NULL),
                       contrasts = list(NULL, NULL),
                       na.action = list(NULL, NULL),
                       se = list(TRUE, TRUE),
                       verbose = FALSE,
                       doFit = TRUE,
                       control = glmmTMBControl(),
                       REML = FALSE,
                       start = NULL,
                       map = NULL,
                       sparseX = NULL)


model_settings   <- tibble(model = c("m1", "m2"),
           formula = c(counts ~  time + time:condition + (1|participant),
           counts ~  time + time:condition + (1|participant)),
           family = c("nbinom2", "nbinom"))


nrow(model_settings)


#' Fits model(s) for a single target
#'
#' @param model_settings A list containing two names lists. 'formula' should contain model formula for each model to be fitted written for glmmTMB::glmmTMB. 'family' should contain family functions (or character description). 'formula' and 'family' should be equal lengths and submitted in the same order.
#' @param data A data frame containing relevant formatted variables for the provided formula/family.
#' @return
#' @export
singlefit <- function(model_settings = list(), summary_fun = NULL, data) {

  # Check if model_settings is a list
  if(!is.list(model_settings)) stop("model_settings must be a list containing at least one pair of model formula and family, see ?singlefit for details")

  # Check if model_settings contains sublists with correct names and matching number of components
  if(!(length(model_settings$formula) == length(model_settings$family))) stop("model_settings must be a list containing at least one pair of model formula and family, see ?singlefit for details")

  # Set the number of models based on model_settings
  nmodels <- length(model_settings$formula)


  # Create a storage for models
  fitted_models <- list()

  # Loop over models to
  for(i in 1:nmodels) {

    fitted_models[[i]] <- tryCatch(
      {
         glmmTMB::glmmTMB(formula = model_settings$formula[[i]],
                                                family = model_settings$family[[i]],
                                                data = data)

      },
      error = function(cond) {
        return(NA)
      }
    )
  }


  if(is.null(summary_fun)) {





  }




  return(fitted_models)



}





#' Iterative fitting over multiple targets in a RNA-seq data set
#'
#' @param model_settings A list containing two names lists. 'formula' should contain model formula for each model to be fitted written for glmmTMB::glmmTMB. 'family' should contain family functions (or character description). 'formula' and 'family' should be equal lengths and submitted in the same order.
#'




itfit <- function(model_settings = list(), meta_data, count_data) {

  # Prepare data





}



