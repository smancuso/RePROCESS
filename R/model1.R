#' PROCESS (v3.0) Model 1
#'
#' Implementation of the PROCESS (v3.0) Model 1 for moderation and conditional effects based on Hayes (2018).
#' Automatically generates interaction plots and Johnson-Neyman plots of areas of significance.
#' Predictor or moderating variables can be continuous, binary, or multicategorical. However, both cannot be multicategorical.
#'
#' @author Sam Mancuso
#' @param Y Dependent variable.
#' @param X Independent variable.
#' @param W Moderating variable.
#' @param CV Covariates. Defaults to \code{NULL}.
#' @param ci Confidence interval range. Defaults to \code{0.95}.
#' @param B Number of bootstrap resamples. Defaults to \code{5000}.
#' @param seed Seed value to replicate analyses. Defaults to \code{1234}.
#' @param boot.inf Obtains bootstrapped confidence intervals for regression. Defaults to \code{FALSE}.
#' @param digits Number of digits to output. Defaults to \code{4}.
#' @param mean.center Mean centre the predictors \code{X} and \code{W}. Defaults to \code{FALSE}.
#' @param se.HC Use HC consistent standard errors. Options are: \code{c("none", "HC0", "HC1", "HC2", "HC3", "HC4")}. Defaults to \code{"none"}.
#' @param p.probe Maximum \emph{p}-value of interaction term at which conditional effects should be probed. Defaults to \code{.10}.
#' @param cond.type Conditioning values for \code{W}. Options are \code{"perc"} for percentiles (\code{.16, .50, .84}) or \code{"mean"} for \code{M, M +/- SD}. Defaults to \code{"perc"}.
#' @param JN Calculates the Johnson-Neyman areas of statistical significance. Defaults to \code{FALSE}.
#' @param plot.JN Plots the Johnson-Neyman areas of statistical significance. Defaults to \code{FALSE}.
#' @param plot.int Plots the interaction plot. Defaults to \code{FALSE}.
#' @param x.multicat Logical value to identify if \code{X} is a multicategorical variable. Defaults to \code{FALSE}.
#' @param x.code Coding format for multicategorical \code{X}. Options are \code{c("indicator", "sequential", "helmert", "effect")}. Defaults to \code{"indicator"}.
#' @param w.multicat Logical value to identify if \code{W} is a multicategorical variable. Defaults to \code{FALSE}.
#' @param w.code Coding format for multicategorical \code{W}. Options are \code{c("indicator", "sequential", "helmert", "effect")}. Defaults to \code{"indicator"}.
#' @param data Name of data.frame or tibble.
#' @keywords PROCESS macro, Moderation, Conditional Effects, Model 1.
#' @import dplyr
#' @import rlang
#' @importFrom broom tidy
#' @importFrom car Anova
#' @importFrom magrittr %>%
#' @importFrom lmtest coeftest waldtest
#' @importFrom sandwich vcovHC
#' @export model1
#' @references Hayes, A. F. (2018). \emph{Introduction to Mediation, Moderation, and Conditional Process Analysis} (2e). The Guilford Press: New York.

model1 <- function(Y,
                   X,
                   W,
                   CV = NULL,
                   ci = .95,
                   B = 5000,
                   boot.inf = FALSE,
                   seed = 1234,
                   digits = 4,
                   probe.p = .10,
                   mean.center = FALSE,
                   se.HC = c("none", "HC0", "HC1", "HC2", "HC3", "HC4"),
                   x.multicat = FALSE,
                   x.code = c("indicator", "sequential", "helmert", "effect"),
                   w.multicat = FALSE,
                   w.code = c("indicator", "sequential", "helmert", "effect"),
                   cond.type = c("perc", "mean"),
                   JN = FALSE,
                   plot.JN = FALSE,
                   plot.int = FALSE,
                   data) {

  # Helper functions ----
  # Sequential coding
  contr.seq <- function(n) {
    # n = number of categories
    # Create matrix with all zero
    seq_contrast <- matrix(
      0,
      nrow = n,
      ncol = n - 1
      )

    # Apply sequential coding
    for (mat_col in 1:n-1) {
      for (mat_row in 2:n) {
        if(mat_row >= mat_col + 1) seq_contrast[mat_row, mat_col] <- 1
      }
    }

    rownames(seq_contrast) <- seq(1:n)

    seq_contrast

  }

  # Data Check ----
  # Ensure specified variables exist in data object
  var_spec <- c(Y, X, W, CV)
  var_names <- names(data)

  # Specified variables that don't exist in data object
  var_miss <- var_spec[which(!(var_spec %in% var_names))]

  # If variables don't exist, stop the function
  if (length(var_miss) > 0) {

    # List missing variables
    stop_msg <- paste0(
      "The following variables were not found: ",
      paste0(var_miss, collapse = ", ")
    )

    stop(stop_msg)

  }

  # Check whether required arguments have been specified
  if (is.null(Y)) stop("Y not specified")
  if (is.null(X)) stop("X not specified")
  if (is.null(W)) stop("W not specified")

  # Check if X is numeric or labelled
  X_class <- data %>%
    pull(X) %>%
    class

  if (!(X_class %in% c("numeric", "labelled"))) {
    stop("Class of X is not numeric or labelled")
  }

  # Check if W is numeric or labelled
  W_class <- data %>%
    pull(W) %>%
    class

  if (!(W_class %in% c("numeric", "labelled"))) {
    stop("Class of X is not numeric or labelled")
  }

  # Check if both X and W are multicategorical
  # Stop if both are multicategorical
  if (all(c(x.multicat, w.multicat))) {

    stop("Both X and W are multicategorical:
         \nFunction currently only supports one multicategorical variable.")

  }

  # Data Wrangling ----
  data_complete <- data %>%
    select(Y, X, W, CV) %>%
    filter(complete.cases(.))

  # Mean-center predictors if requested
  if (mean.center) {
    # Do not mean centre if X is a categorical predictor
    if (!x.multicat) {
      # Mean-centre X
      data_complete <- data_complete %>%
        mutate(!!sym(X) := !!sym(X) - mean(!!sym(X)))
    }
    # Do not mean centre if W is a multicategorical predictor
    if (!w.multicat) {
      # Mean-centre W
      data_complete <- data_complete %>%
        mutate(!!sym(W) := !!sym(W) - mean(!!sym(W)))
    }
  }

  ## X contrasts
  # Get number of categories
  Xk <- data_complete %>%
    pull(X) %>%
    unique %>%
    as.vector %>%
    length

  # Set contrasts only if X is multicategorical with 3 or more categories
  if(x.multicat & Xk > 2) {

    # Convert X to factor
    data_complete <- data_complete %>%
      mutate(!!sym(X) := factor(!!sym(X))) %>%
      # Convert to data.frame
      as.data.frame

    # Indicator (Dummy/Treatment)
    if (x.code == "indicator") contrasts(data_complete[, X]) <- contr.treatment(Xk)
    # Sequential (Custom function)
    if (x.code == "sequential") contrasts(data_complete[, X]) <- contr.seq(Xk)
    # Helmert does not need to be changed
    if (x.code == "helmert") contrasts(data_complete[, X]) <- contr.helmert(Xk)
    # Effect (Deviation)
    if (x.code == "effect") contrasts(data_complete[, X]) <- contr.sum(Xk)

  }

  ## W contrasts
  # Get number of categories
  Wk <- data_complete %>%
    pull(W) %>%
    unique %>%
    as.vector %>%
    length

  # Set contrasts only if W is multicategorical with 3 or more categories
  if(w.multicat & Wk > 2) {

    # Convert W to factor
    data_complete <- data_complete %>%
      mutate(!!sym(W) := factor(!!sym(W))) %>%
      # Convert to data.frame
      as.data.frame

    # Indicator (Dummy/Treatment)
    if (w.code == "indicator") contrasts(data_complete[, W]) <- contr.treatment(Wk)
    # Sequential (Custom function)
    if (w.code == "sequential") contrasts(data_complete[, W]) <- contr.seq(Wk)
    # Helmert does not need to be changed
    if (w.code == "helmert") contrasts(data_complete[, W]) <- contr.helmert(Wk)
    # Effect (Deviation)
    if (w.code == "effect") contrasts(data_complete[, W]) <- contr.sum(Wk)

  }

  # Regressions ----
  # NB Fit using do.call to prevent errors with waldtest()
  ## Model without the interaction term
  # Generate formula
  formula_base <- paste0(
    Y,
    "~",
    paste0(c(X, W, CV), collapse = "+")
  )

  # Fit regression
  fit_base <- do.call(
    "lm",
    list(
      formula = as.formula(formula_base),
      data = data_complete
      )
    )

  # Model with interaction term
  formula_int <- paste0(formula_base, "+", X, ":", W)

  # Fit regression
  fit_int <- do.call(
    "lm",
    list(
      formula = as.formula(formula_int),
      data = data_complete
    )
  )

  # Compare R^22 changes (NB: R^2 does not change if using robust SEs)
  r2_change = summary(fit_int)$r.squared - summary(fit_base)$r.squared

  # Use heteroskedasticity-consistent standard errors
  # Change "none" to "const" for the vcovHC function
  se.HC <- ifelse(se.HC == "none", "const", se.HC)

  # Base model
  fit_base_hc <- coeftest(
    fit_base,
    vcov. = vcovHC(fit_base, type = se.HC)
  )

  # Get HC-adjusted standard errors and t-tests
  fit_int_hc <- coeftest(
    fit_int,
    vcov. = function(x) {vcovHC(x, type = se.HC)}
  )

  # Get adjusted F-statistic
  fit_int_anova <- waldtest(
    fit_int,
    vcov = function(x) {vcovHC(x, type = se.HC)}
  )

  # Compare base and interaction models
  comp_fit <- waldtest(
    fit_base,
    fit_int,
    vcov = function(x) {vcovHC(x, type = se.HC)}
  )

  ## Tidy output
  # Extract model terms
  model_terms <- fit_int_hc %>%
    tidy %>%
    pull(term)

  # Add HC-adjusted confidence intervals
  fit_int_hc <- fit_int_hc %>%
    tidy %>%
    mutate(
      LLCI = estimate - qt(1 - ((1 - ci) / 2), fit_int$df.residual) * std.error,
      ULCI = estimate + qt(1 - ((1 - ci) / 2), fit_int$df.residual) * std.error
    ) %>%
    select(-term) %>%
    round(., digits) %>%
    mutate(term = model_terms) %>%
    select(term, estimate, std.error, statistic, p.value, LLCI, ULCI) %>%
    mutate(p.value = ifelse(p.value < .001, "<0.001", p.value)) %>%
    setNames(., c("", "coeff", "se", "t", "p", "LLCI", "ULCI"))

  # Model Summary
  model_summary <- data.frame(
    R = sqrt(summary(fit_int)$r.squared),
    R.sq = summary(fit_int)$r.squared,
    MSE = summary(fit_int)$sigma^2,
    F.value = fit_int_anova[2, 3],
    df1 = abs(fit_int_anova[2, 2]),
    df2 = fit_int_anova[1, 1],
    p = fit_int_anova[2, 4]
  ) %>%
    round(., digits) %>%
    mutate(p = ifelse(p < .001, "<0.001", p))

  # Interaction Summary
  int_summary <- data.frame(
    R2.change = r2_change,
    F.value = comp_fit[2, 3],
    df1 = comp_fit[2, 2],
    df2 = comp_fit[2, 1],
    p = comp_fit[2, 4]
  ) %>%
    round(digits) %>%
    mutate(p = ifelse(p < .001, "<0.001", p)) %>%
    mutate(term = paste0(X, " x ", W)) %>%
    select(term, R2.change, F.value, df1, df2, p) %>%
    setNames(., c("", "R2.change", "F.value", "df1", "df2", "p"))

  # Output Results ----
  cat(
    "\n", "---------------------------------------------------------------------",
    "\n", "PROCESS Procedure for R Version 0.1",
    "\n", "Based on SPSS Macro by Andrew F. Hayes  www.afhayes.com",
    "\n", "Coded by Sam Mancuso   sammancuso.com",
    "\n", "---------------------------------------------------------------------",
    "\n", "Model : 1",
    "\n", "    Y : ", Y,
    "\n", "    X : ", X,
    "\n", "    W : ", W,
    "\n",
    "\n", "Covariates:",
    "\n", " ", ifelse(is.null(CV),"None", paste0(CV, collapse = " ")),
    "\n",
    "\n", "Sample size: ", nrow(data_complete),
    sep = ""
    )
  if(w.multicat) {
    # Generate matrix
    cond_values <- contrasts(data_complete[, W])
    colnames(cond_values) <- paste0("W", c(1:(Wk - 1)))
    rownames(cond_values) <- data_complete %>%
      pull(W) %>%
      levels

    cat(
      "\n",
      "\n", "Coding of categorical W variable (", W, ") for analysis:",
      "\n",
      "\n",
      sep = ""
      )
    print(cond_values)
  }

  if(x.multicat) {
    # Generate matrix
    cond_values <- contrasts(data_complete[, X])
    colnames(cond_values) <- paste0("X", c(1:(Xk - 1)))
    rownames(cond_values) <- data_complete %>%
      pull(X) %>%
      levels

    cat(
      "\n",
      "\n", "Coding of categorical X variable (", X, ") for analysis:",
      "\n",
      "\n",
      sep = ""
    )
    print(cond_values)
  }

  cat(
    "\n", "---------------------------------------------------------------------",
    "\n",
    "\n", "OUTCOME VARIABLE:",
    "\n", " ", Y,
    "\n",
    "\n", "Model Summary",
    "\n",
    sep = ""
  )

  # Print Model Summary
  print(model_summary, row.names = FALSE)

  cat(
    "\n",
    "Model",
    "\n",
    sep = ""
  )

  print(fit_int_hc, row.names = FALSE)

  cat(
    "\n",
    "Test(s) of unconditional interaction(s):",
    "\n",
    sep = ""
  )

  # Print Interaction Test
  print(int_summary, row.names = FALSE)

  # Conditional Effects ----
  # Probe interaction if interaction p-value < probe.p (default = .10)
  if (comp_fit[2, 4] < probe.p) {

    ## Binary W and continuous X
    if (Wk == 2 & X_class %in% c("numeric", "labelled")) {
      # Conditional effects
      cond.effect(
        w.type = "bin",
        x.type = "cont",
        mod = fit_int,
        Y = Y,
        X = X,
        W = W,
        CV = CV,
        se.HC = se.HC,
        ci = ci,
        digits = digits,
        JN.req = JN,
        plot.JN = plot.JN,
        plot.int = plot.int
      )
    }

    ## Multicategorical W and continuous X
    if (Wk > 2 & w.multicat & X_class %in% c("numeric", "labelled")) {
      # Conditional effects
      cond.effect(
        w.type = "multicat",
        x.type = "cont",
        mod = fit_int,
        Y = Y,
        X = X,
        W = W,
        CV = CV,
        se.HC = se.HC,
        ci = ci,
        digits = digits,
        JN.req = JN,
        plot.JN = plot.JN,
        plot.int = plot.int
      )
    }

    ## Continuous X and W
    if (
      all(
        W_class == "numeric",
        X_class == "numeric",
        x.multicat == FALSE,
        w.multicat == FALSE,
        Xk > 2,
        Wk > 2
        )
      ) {
      # Conditional effects
      cond.effect(
        w.type = "cont",
        x.type = "cont",
        cond.type = cond.type,
        mod = fit_int,
        Y = Y,
        X = X,
        W = W,
        CV = CV,
        se.HC = se.HC,
        ci = ci,
        digits = digits,
        JN.req = JN,
        plot.JN = plot.JN,
        plot.int = plot.int
      )
    }

    ## Binary X and Continuous W
    if (W_class %in% c("numeric", "labelled") & Xk == 2) {
      # Conditional effects
      cond.effect(
        w.type = "cont",
        x.type = "bin",
        cond.type = cond.type,
        mod = fit_int,
        Y = Y,
        X = X,
        W = W,
        CV = CV,
        se.HC = se.HC,
        ci = ci,
        digits = digits,
        JN.req = JN,
        plot.JN = plot.JN,
        plot.int = plot.int
      )
    }

    ## Multicategorical X and continuous W
    if (
      all(
        Xk > 2,
        x.multicat,
        w.multicat == FALSE,
        W_class %in% c("numeric", "labelled")
        )
      ) {
      # Conditional effects
      cond.effect(
        x.type = "multicat",
        w.type = "cont",
        mod = fit_int,
        Y = Y,
        X = X,
        W = W,
        CV = CV,
        se.HC = se.HC,
        ci = ci,
        digits = digits,
        JN.req = JN,
        plot.JN = plot.JN,
        plot.int = plot.int
        )
    }
  }

  # Message for W values in conditional tables
  W_message <- ""

  if (W_class %in% c("numeric")) {

    if (cond.type == "perc") {

      W_message <- "\nW values in conditional tables are the 16th, 50th, and 84th percentiles.\n"

    }

    if (cond.type == "mean") {

      W_temp <- data_complete %>%
        pull(W)

      if (mean(W_temp) - sd(W_temp) < min(W_temp)) {

        W_message <- "W values in conditional tables are the minimum, mean and 1 SD above the mean.

Note: One SD below the mean is below the minimum observed in the data for W,
      so the minimum measurement on W is used for conditioning instead.\n"

      } else{

        W_message <- "\nW values in conditional tables are the mean and +/- SD from the mean.\n"

      }
    }

  }

  # Bootstrapping if requested ----
  if (boot.inf) {

    # Set seed
    set.seed(seed)

    # Get regression formula
    boot_formula <- formula(fit_int)

    # Get model data
    data_model <- model.frame(fit_int)

    # Get number of observations
    n_obs <- nrow(data_model)

    # Bootstrap function
    boot.lm <- function(formula, data, n) {
      lm(
        formula,
        data = data[sample.int(n, replace = TRUE), ]
        )$coef
    }

    # Bootstrap coefficients
    mat_boot <- t(
      replicate(
        B,
        boot.lm(formula = boot_formula, data = data_model, n = n_obs)
        )
      )

    # Bootstrap output table
    boot_out <- data.frame(
      Coeff = coef(fit_int),
      BootMean = apply(mat_boot, 2, mean),
      BootSE = apply(mat_boot, 2, sd),
      BootLLCI = apply(mat_boot, 2, function(x) quantile(x, (1 - ci) / 2)),
      BootULCI = apply(mat_boot, 2, function(x) quantile(x, 1 - ((1 - ci) / 2)))
      ) %>%
      round(., digits)

    cat(
      "\n", "--------- BOOTSTRAP RESULTS FOR REGRESSION MODEL PARAMETERS ---------",
      "\n",
      "\n", "OUTCOME VARIABLE:",
      "\n", " ", Y,
      "\n",
      sep = ""
      )

    print(boot_out)

  }

  cat(
    "\n", "--------------------- ANALYSIS NOTES AND ERRORS ---------------------",
    "\n",
    "\n", "Level of confidence for all confidence intervals in output:",
    "\n  ", format(round(100 * ci, digits), nsmall = digits),
    "\n",
    ifelse(
      boot.inf,
      paste0("\nNumber of bootstrap samples for percentile confidence intervals:\n  ", B),
      ""
      ),
    "\n", W_message,
    ifelse(
      se.HC != "const",
      paste0(
        "\nNOTE: The ", toupper(se.HC), " standard error and covariance matrix estimator was used."
        ),
      ""
      ),
    "\n",
    sep = ""
    )
  # List of mean-centred variables
  cent_var <- NULL

  if (mean.center) {
    # Do not mean centre if X is a categorical predictor
    if (!x.multicat) {
      # Mean-centre X
      cent_var <- c(cent_var, X)
    }
    # Do not mean centre if W is a multicategorical predictor
    if (!w.multicat) {
      # Mean-centre Y
      cent_var <- c(cent_var, W)
    }
  }

  if (!is.null(cent_var)) {

    cat(
      "\n", "NOTE: The following variables were mean centered prior to analysis: ",
      "\n", "\t", paste0(cent_var, collapse = " "),
      "\n",
      sep = ""
      )

  }

}
