library(flexmix)
library(ggplot2)

#' Choose the Best Trajectory Model
#'
#' This function fits different Group-Based Trajectory Models (GBTM)
#' to the provided data and ranks them based on Bayesian Information Criterion (BIC),
#' Entropy Information Criterion (EIC), Average Posterior Probability of Assignment (APPA),
#' and Odds of Correct Classification (OCC).
#'
#' @param data A data frame containing the variables needed for the model.
#' @param id_col The name of the column representing individual IDs.
#' @param treatment_col The name of the column representing the treatment outcome (binary).
#' @param time_col The name of the column representing time points.
#' @param J The number of groups (trajectories) to be fitted in the models.
#' @param additional_covariates (optional) A string representing additional covariates to be included in the model.
#' @return A list of model names ranked from best to worst based on the specified criteria.
#' @examples
#' # Assuming 'data' is your dataframe:
#' best_models <- choose_best_model(data, "id", "treatment", "time", 4, "oxygen")
#' @export

choose_best_model <- function(data, treatment_col, time_col, J, additional_covariates = NULL) {
  results <- data.frame(Model = character(),
                        BIC = numeric(),
                        EIC = numeric(),
                        APPA = numeric(),
                        OCC = numeric(),
                        Score = numeric(),
                        stringsAsFactors = FALSE)

  for (order in 1:3) {
    # Base formula part
    base_formula_part <- paste("cbind(", treatment_col, ",1-", treatment_col, ") ~ I(", time_col, ")")

    # Update formula for quadratic and cubic terms
    formula <- base_formula_part
    if (order >= 2) {
      formula <- paste(formula, "+ I(", time_col, "^2)")
    }
    if (order == 3) {
      formula <- paste(formula, "+ I(", time_col, "^3)")
    }

    # Add additional covariates if provided
    if (!is.null(additional_covariates)) {
      formula <- paste(formula, "+", additional_covariates)
    }

    model_formula <- as.formula(formula)

    # Fit the model
    lcgaMix <- try(flexmix(model_formula, data = data, k = J, model = FLXMRglm(family = "binomial"), control = list(iter.max = 1000, minprior = 0)), silent = TRUE)
    if (class(lcgaMix) == "try-error") {
      next
    }

    # Calculate metrics
    BIC_value <- BIC(lcgaMix)
    EIC_value <- EIC(lcgaMix)
    dpost <- data.frame(posterior(lcgaMix))
    dpost$class <- lcgaMix@cluster
    APPA_value <- mean(diag(as.matrix(aggregate(. ~ class, FUN = mean, data = dpost))[, -1]))
    prior <- lcgaMix@prior
    OCC_value <- mean((APPA_value / (1 - APPA_value)) / (prior / (1 - prior)))
    Score <- (BIC_value-1)/BIC_value + EIC_value + APPA_value + OCC_value

    # Add to results
    results <- rbind(results, data.frame(Model = paste("Model Polynomial Order", order), BIC = BIC_value, EIC = EIC_value, APPA = APPA_value, OCC = OCC_value, Score = Score))
  }

  # Rank models
  results$Rank <- rank(-results$Score) # Lower rank number is better
  return(results)
}

# Example usage
source("generate_trajdata")
dat.rshAL <- generate_group_data(n = 2000)
best_models <- choose_best_model(data = dat.rshAL, treatment_col = "A", time_col = "time", J = 4, additional_covariates = c("L","V"))
best_models
