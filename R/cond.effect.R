#' PROCESS (v3.0) Conditional Effects and Plots for Model 1
#'
#' Implementation of the PROCESS (v3.0) Model 1 conditional effects based on Hayes (2018).
#' Automatically generates interaction plots and Johnson-Neyman plots of areas of significance.
#' Predictor or moderating variables can be continuous, binary, or multicategorical. However, both cannot be multicategorical.
#'
#' @author Sam Mancuso
#' @param mod Linear model of class \code{"lm"}.
#' @param Y Dependent variable.
#' @param X Independent variable.
#' @param W Moderating variable.
#' @param CV Covariates. Defaults to \code{NULL}.
#' @param digits Number of digits to output. Defaults to \code{4}.
#' @param se.HC Use HC consistent standard errors. Options are: \code{c("none", "HC0", "HC1", "HC2", "HC3", "HC4")}. Defaults to \code{"none"}.
#' @param p.probe Maximum \emph{p}-value of interaction term at which conditional effects should be probed. Defaults to \code{.10}.
#' @param cond.type Conditioning values for \code{W}. Options are \code{"perc"} for percentiles (\code{.16, .50, .84}) or \code{"mean"} for \code{M, M +/- SD}. Defaults to \code{"perc"}.
#' @param JN.req Calculates the Johnson-Neyman areas of statistical significance. Defaults to \code{FALSE}.
#' @param plot.JN Plots the Johnson-Neyman areas of statistical significance. Defaults to \code{FALSE}.
#' @param plot.int Plots the interaction plot. Defaults to \code{FALSE}.
#' @param x.type Format of \code{X} variable: binary, multicategorical, or continuous. Corresponding options are \code{c("bin", "multicat", "cont")}.
#' @param w.type Format of \code{W} variable: binary, multicategorical, or continuous. Corresponding options are \code{c("bin", "multicat", "cont")}.#' @param data Name of data.frame or tibble.
#' @keywords PROCESS macro, Moderation, Conditional Effects, Model 1.
#' @import dplyr
#' @import rlang
#' @importFrom broom tidy
#' @importFrom car Anova
#' @importFrom magrittr %>%
#' @importFrom lmtest coeftest
#' @export cond.effect
#' @references Hayes, A. F. (2018). \emph{Introduction to Mediation, Moderation, and Conditional Process Analysis} (2e). The Guilford Press: New York.

cond.effect <- function(mod,
                        Y, X, W, CV,
                        data, se.HC,
                        ci, digits,
                        x.type = c("bin", "multicat", "cont"),
                        w.type = c("bin", "multicat", "cont"),
                        cond.type = "perc",
                        plot.int = FALSE,
                        plot.JN = FALSE,
                        JN.req = FALSE) {

  # Extract data from model
  data_model <- model.frame(mod)

  # Extract terms of interest (X and X:W)
  mod_sub <- mod %>%
    coeftest(., vcov. = vcovHC(., type = se.HC)) %>%
    tidy %>%
    filter(
      term %in%
        grep(
          paste0(c(X, ":"), collapse = "|"),
          term,
          value = TRUE
        )
    ) %>%
    select(term, estimate, std.error)

  # Name of interaction term(s)
  int_term <- mod_sub %>%
    filter(term %in% grep(":", term, value = TRUE)) %>%
    pull(term)

  # Get X term(s)
  x_term <- mod_sub %>%
    filter(
      term %in% intersect(
        grep(X, term, value = TRUE),
        grep(":", term, value = TRUE, invert = TRUE)
      )
    ) %>%
    pull(term)

  # Continuous X ----
  if (x.type == "cont") {

    # Binary W ----
    if (w.type == "bin") {

      # Get W-values
      cond_values <- data_model %>%
        pull(W) %>%
        unique

      cond_values <- cond_values %>%
        as.vector

      # Get beta values
      b1 <- coef(mod)[X]
      b3 <- coef(mod)[int_term]


      # Get conditional effects
      cond_effects <- b1 + cond_values * b3

      # Calculate Standard Errors for conditional effects
      # Get variance-covariance matrix
      mod_cov <- vcovHC(mod, type = se.HC)

      Vb1 <- mod_cov[X, X]
      Vb3 <- mod_cov[int_term, int_term]
      Vb1b3 <- mod_cov[X, int_term]

      # Create blank vector
      se_cond <- sqrt(Vb1 + (Vb3 * cond_values ^ 2) + (2 * Vb1b3 * cond_values))
    }

    # Multicategorical W ----
    if (w.type == "multicat") {
      # Get beta coefficients
      b1 <- coef(mod)[X]
      bw <- coef(mod)[int_term]

      # Get W levels
      W_levels <- data_model %>%
        pull(W) %>%
        levels %>%
        as.numeric

      # Create contrast matrix
      cond_values <- contrasts(data_model[, W])

      # Use matrix multiplication to get conditional effects
      cond_effects <- b1 + cond_values %*% as.matrix(bw) %>%
        as.vector

      # Calculate Standard Errors for conditional effects
      # Get number of categories
      k <- length(cond_effects)

      # Get variance-covariance matrix
      mod_cov <- vcovHC(mod, type = se.HC)

      Vb1 <- mod_cov[X, X]
      Vbw <-vector(mode = "numeric", length = k  - 1)
      Vb1bw <-vector(mode = "numeric", length = k  - 1)

      for(i in 1:length(int_term)) {

        Vbw[i] <- mod_cov[int_term[i], int_term[i]]
        Vb1bw[i] <- mod_cov[X, int_term[i]]

      }

      # SE of conditional values
      se_cond <- sqrt(
        Vb1 +
          rowSums(t(t(cond_values ^ 2) * Vbw)) +
          2 * rowSums(t(t(cond_values) * Vb1bw))
      )

    }

  }

  # Continuous W ----
  if (w.type == "cont") {

    # Get W-values
    cond_values <- data_model %>%
      pull(W)

    # Conditioning values
    if (cond.type == "perc") {
      # Use percentiles. NB Does not correspond to SPSS percentiles.
      cond_values <- quantile(cond_values, c(.16, .50, .84))
    }

    if (cond.type == "mean") {

      # Check if M - 1SD is below minimum value of W
      if (mean(cond_values) - sd(cond_values) < min(cond_values)) {

        # Use minimum measurement of W instead
        cond_values <- c(
          min(cond_values),
          mean(cond_values),
          mean(cond_values) + sd(cond_values)
        )

      } else if (mean(cond_values) + sd(cond_values) > max(cond_values)) {

        # Use minimum measurement of W instead
        cond_values <- c(
          mean(cond_values) - sd(cond_values),
          mean(cond_values),
          max(cond_values)
        )

      } else {

        cond_values <- c(
          mean(cond_values) - sd(cond_values),
          mean(cond_values),
          mean(cond_values) + sd(cond_values)
        )

      }
    }

    cond_values <- cond_values %>%
      unname

    # Set W-levels as conditioning values
    W_levels <- cond_values

    # Continuous or Binary X ----
    if (x.type %in% c("cont", "bin")) {

      # Get beta values
      # Get beta coefficients
      b1 <- coef(mod)[X]
      b3 <- coef(mod)[int_term]

      # Get conditional effects
      cond_effects <- b1 + cond_values * b3

      # Calculate Standard Errors for conditional effects
      # Get variance-covariance matrix
      mod_cov <- vcovHC(mod, type = se.HC)

      Vb1 <- mod_cov[X, X]
      Vb3 <- mod_cov[int_term, int_term]
      Vb1b3 <- mod_cov[X, int_term]

      # Create blank vector
      se_cond <- sqrt(Vb1 + (Vb3 * cond_values ^ 2) + (2 * Vb1b3 * cond_values))

      # # Johnson-Neyman significance regions ----
      if (JN.req) {
        # Get critical t-value
        tcrit <- qt(1 - (1 - ci) / 2, mod$df.residual)

        a <- tcrit ^ 2 * Vb3 - b3 ^ 2
        b <- 2 * (tcrit ^ 2 * Vb1b3 - b1 * b3)
        c <- tcrit ^ 2 * Vb1 - b1 ^ 2

        jn <- c(
          (-b + sqrt(b ^ 2 - 4 * a * c))/ (2 * a),
          (-b - sqrt(b ^ 2 - 4 * a * c))/ (2 * a)
        ) %>%
          unname %>%
          sort

        perc.below <- c(
          100 * sum(data_model %>% pull(W) < jn[1]) / nrow(data_model),
          100 * sum(data_model %>% pull(W) < jn[2]) / nrow(data_model)
        )

        perc.above <- c(
          100 * sum(data_model %>% pull(W) > jn[1]) / nrow(data_model),
          100 * sum(data_model %>% pull(W) > jn[2]) / nrow(data_model))

        jn_out <- data.frame(
          Value = jn,
          perc.below,
          perc.above
        ) %>%
          round(., digits) %>%
          setNames(., c("Value", "% below", "% above"))

        # Values of W
        W_data <- pull(data_model, W)

        W_data <- seq(
          min(W_data),
          max(W_data),
          (max(W_data) - min(W_data) + 1) / 21
        )

        # Insert jn values into W_data
        if (all(perc.above == 100) | all(perc.below == 100)) {
          W_data <- sort(W_data)
        } else {
          W_data <- sort(c(jn, W_data))
        }

        # Create table to output JN results
        jn_cond_out <- tibble(
          !!sym(W) := W_data,
          Effect = b1 + !!sym(W) * b3,
          se = sqrt(Vb1 + (Vb3 * (!!sym(W)) ^ 2) + (2 * Vb1b3 * !!sym(W))),
          t = Effect / se,
          p = 2 * pt(-abs(t), df = mod$df.residual),
          LLCI = Effect - se * qt(1 - ((1 - ci) / 2), mod$df.residual),
          ULCI = Effect + se * qt(1 - ((1 - ci) / 2), mod$df.residual)
        ) %>%
          round(., digits) %>%
          as.data.frame
      }
    }
  }

  # Output Results ----------------------------------------------------------
  # Multicategorical X ----
  if (x.type == "multicat") {

    # Print multicategorical X
    # Get X values
    X_values <- data_model %>%
      pull(X) %>%
      levels

    # Create blank list to store conditional effects
    cond_effects <- vector(mode = "numeric", length = length(x_term))

    # Loop through conditioning values
    for (i in 1:length(cond_values)) {

      # Centre W at conditioning values
      data_temp <- data_model %>%
        mutate(!!sym(W) := !!sym(W) - cond_values[i])

      # Fit full model
      fit_full <- do.call(
        "lm",
        list(
          formula = formula(mod),
          data = data_temp
        )
      )

      # Get HC-adjusted F-statistic
      # waldtest will not work, so a work-around with Anova from the 'car'
      # package is used
      comp_fit <- suppressMessages(
        car::Anova(fit_full,
                   type = 3,
                   vcov. = vcovHC(fit_full, type = se.HC)
        )
      )

      # Get residual df
      comp_fit_df_resid <- comp_fit %>%
        tidy %>%
        filter(term == "Residuals") %>%
        pull(df)

      # Get df, F, and p for X
      comp_fit_x <- comp_fit %>%
        tidy %>%
        filter(term == X)

      comp_fit_out <- data.frame(
        "F.value" = comp_fit_x$statistic,
        df1 = comp_fit_x$df,
        df2 = comp_fit_df_resid,
        p = comp_fit_x$p.value
      ) %>%
        round(., digits) %>%
        mutate(p = ifelse(p < .001, "<0.001", p))

      # Conditional effects
      cond_effects <- fit_full %>%
        coeftest(., vcov. = vcovHC(., type = se.HC)) %>%
        tidy %>%
        filter(term %in% x_term) %>%
        select(-term) %>%
        rename(
          Effect = estimate,
          se = std.error,
          t = statistic,
          p = p.value
        ) %>%
        mutate(
          LLCI = Effect - se * qt(1 - ((1 - ci) / 2), fit_full$df.residual),
          ULCI = Effect + se * qt(1 - ((1 - ci) / 2), fit_full$df.residual)
        ) %>%
        round(., digits) %>%
        mutate(p = ifelse(p < .001, "<0.001", p)) %>%
        mutate(term = x_term) %>%
        select(term, Effect, se, t, p, LLCI, ULCI) %>%
        setNames(., c("", "Effect", "se", "t", "p", "LLCI", "ULCI"))

      # Esitmated conditional means
      ecm_data <- tibble(
        !!sym(X) := X_values,
        !!sym(W) := cond_values[i]
      ) %>%
        as.data.frame

      ecm <- predict(mod, newdata = ecm_data)

      ecm_out <- tibble(
        !!sym(X) := as.numeric(X_values),
        !!sym(Y) := ecm
      ) %>%
        as.data.frame %>%
        round(., digits)

      # Print results to screen
      cat(
        "\n",
        "\n", "Moderator value(s):",
        "\n", W, "\t", format(round(cond_values[i], digits), nsmall = digits),
        "\n",
        "\n",
        sep = ""
      )

      print(cond_effects, row.names = FALSE)

      cat(
        "\n",
        "\n", "Test of conditional means",
        "\n",
        sep = ""
      )

      print(comp_fit_out, row.names = FALSE)

      cat(
        "\n",
        "\n", "Estimated conditional means being compared:",
        "\n",
        sep = ""
      )

      print(ecm_out, row.names = FALSE)

      cat(
        "----------------",
        "\n",
        sep = ""
      )
    }
  } else {

    cond_out <- tibble(
      !!sym(W) := W_levels,
      Effect = cond_effects,
      se = se_cond
    ) %>%
      mutate(
        t = Effect/se,
        p = 2 * pt(-abs(t), df = mod$df.residual)
      ) %>%
      as.data.frame

    # Output conditional effects
    cond_out <- cond_out %>%
      mutate(
        LLCI = Effect - qt(1 - ((1 - ci) / 2), mod$df.residual) * se,
        ULCI = Effect + qt(1 - ((1 - ci) / 2), mod$df.residual) * se
      ) %>%
      round(., digits) %>%
      as.data.frame

    cat(
      "\n", "--------",
      "\n", "    Focal predictor: ", paste0(X, " (X)"),
      "\n", "          Moderator: ", paste0(W, " (W)"),
      "\n",
      "\n", "Conditional effects of the focal predictor at values of the moderator(s):",
      "\n",
      "\n",
      sep = ""
    )

    print(cond_out, row.names = FALSE)

  }

  # Print Johnson-Neyman Results (W must be continuous)
  if (JN.req & w.type == "cont") {

  if (exists("jn_out")) {

    if (all(perc.above == 100) | all(perc.below == 100)) {
      cat("\n", "There are no statistical significance transition points within
          the observed range of the moderator found using the Johnson-Neyman method.\n",
          sep = "")

    } else {

      cat(
        "\n", "Moderator value(s) defining Johnson-Neyman significance region(s):",
        "\n",
        sep = "")

      print(jn_out, row.names = FALSE)
      }
    }

    if (exists("jn_cond_out")) {

      cat(
        "\n", "Conditional effect of focal predictor at values of the moderator:",
        "\n",
        sep = ""
      )

      print(jn_cond_out, row.names = FALSE)

    }

  }

# Plots -------------------------------------------------------------------
  # Interaction plots ----
  if (plot.int) {

    # Extract data to plot
    Yplot <- pull(data_model, Y)
    Xplot <- pull(data_model, X)
    Wplot <- pull(data_model, W)

    # Continuous W ----
    if (w.type == "cont") {

      # Binary or multicategorical X ----
      if(x.type %in% c("bin", "multicat")) {

        plot(
          range(Wplot),
          range(predict(mod)),
          type = "n",
          xlab = paste0(W, " (W)"),
          ylab = paste0(Y, " (Y)")
        )

        # Check if X is a factor
        if (class(Xplot) %in% "factor") {

          X_unique <- Xplot %>%
            levels %>%
            sort

        } else {

          X_unique <- sort(unique(Xplot))
        }

        # Vectorise legend text
        legend_text <- vector(mode = "character", length = length(X_unique))

        # Generate new data for plotting
        for (i in seq_along(X_unique)) {

          newdata <- tibble(
            !!sym(X) := X_unique[i],
            !!sym(W) := range(Wplot[Xplot == X_unique[i]])
          ) %>%
            as.data.frame

          # Set variable classe for X
          if (class(Xplot) == "factor") {
            newdata[, X] <- factor(newdata[, X])
          }

          if (!is.null(CV)) {

            for (j in 1:length(CV)) {
              # Extract covariate data
              CVplot <- pull(data_model, CV[j])

              # Create covariate
              newdata <- newdata %>%
                mutate(!!sym(CV[j]) := CVplot[Xplot == X_unique[i]]) %>%
                as.data.frame
            }

          }

          # Plot line
          lines(pull(newdata, W), predict(mod, newdata), lty = i)

          # Add legend text
          legend_text[i] <- paste0(X, " = ", X_unique[i], " (X", i, ")")

        }

        # Add legend to plot
        legend(
          "topleft",
          lty = seq_along(X_unique),
          legend = legend_text,
          cex = .60
        )

      }

      if (x.type == "cont") {

        plot(
          range(Xplot),
          range(predict(mod)),
          type = "n",
          xlab = paste0(X, " (X)"),
          ylab = paste0(Y, " (Y)")
        )

        # Generate new data for plotting
        for (i in seq_along(unique(cond_values))) {

          newdata <- cbind(
            unique(Wplot)[i],
            range(Xplot[Wplot == cond_values[i]])
          )

          if (!is.null(CV)) {

            for (j in 1:length(CV)) {

              # Temporary CV
              CVplot <- pull(data_model, CV[j])

              newdata <- cbind(
                newdata,
                range(CVplot[Wplot == cond_values[i]])
              )
            }

          }

          newdata <- newdata %>%
            as.data.frame %>%
            setNames(., c(W, X, CV))

          lines(pull(newdata, X), predict(mod, newdata), lty = i)

        }

        # Legend Text
        if(cond.type == "perc") {

          legend_text <- c(
            paste0(W, " = ", round(cond_values[1], digits), " (W = 16th perc)"),
            paste0(W, " = ", round(cond_values[2], digits), " (W = 50th perc)"),
            paste0(W, " = ", round(cond_values[3], digits), " (W = 84th perc)")
          )

        } else {

          legend_text <- c(
            ifelse(
              mean(Wplot) - sd(Wplot) < min(Wplot),
              paste0(W, " = ", round(min(Wplot), digits), " (W = Min)"),
              paste0(W, " = ", round(cond_values[1], digits), " (W = M - 1SD)")
            ),
            paste0(W, " = ", round(cond_values[2], digits), " (W = M)"),
            paste0(W, " = ", round(cond_values[3], digits), " (W = M + 1SD)")
          )
        }

        legend(
          "topleft",
          lty = seq_along(unique(Wplot)),
          legend = legend_text,
          cex = .60
        )

      }
    }

    # Continuous X and binary or multicategorical W
    if (w.type %in% c("bin", "multicat") & x.type == "cont") {

      plot(
        range(Xplot),
        range(predict(mod)),
        type = "n",
        xlab = paste0(X, " (X)"),
        ylab = paste0(Y, " (Y)")
      )

      # Check if X is a factor
      if (class(Wplot) %in% "factor") {

        W_unique <- Wplot %>%
          levels %>%
          sort

      } else {

        W_unique <- sort(unique(Wplot))
      }

      # Vectorise legend text
      legend_text <- vector(mode = "character", length = length(W_unique))

      # Generate new data for plotting
      for (i in seq_along(W_unique)) {

        newdata <- tibble(
          !!sym(W) := W_unique[i],
          !!sym(X) := range(Xplot[Wplot == W_unique[i]])
        ) %>%
          as.data.frame

        # Set variable classe for X
        if (class(Wplot) == "factor") {
          newdata[, W] <- factor(newdata[, W])
        }

        if (!is.null(CV)) {

          for (j in 1:length(CV)) {
            # Extract covariate data
            CVplot <- pull(data_model, CV[j])

            # Create covariate
            newdata <- newdata %>%
              mutate(!!sym(CV[j]) := CVplot[Wplot == W_unique[i]]) %>%
              as.data.frame
          }

        }

        # Plot line
        lines(pull(newdata, X), predict(mod, newdata), lty = i)

        # Add legend text
        legend_text[i] <- paste0(W, " = ", W_unique[i], " (W", i, ")")

      }

      # Add legend to plot
      legend(
        "topleft",
        lty = seq_along(W_unique),
        legend = legend_text,
        cex = .60
      )

    }

  }

  # Johnson-Neyman plot ----
  if (all(plot.JN, JN.req)) {

    # Generate JN plot
    plot(
      x = W_data, y = jn_cond_out$Effect,
      type = "l", pch = 19,lwd = 1,
      xlab = paste0(W, " (W)"),
      ylab = paste0("Conditional effect of ", X, " (X) on ", Y, " (Y)"),
      xlim = c(min(W_data), max(W_data)),
      ylim = c(min(jn_cond_out$LLCI), max(jn_cond_out$ULCI))

    )
    points(W_data, jn_cond_out$LLCI, lwd = 2, lty = 2, type = "l")
    points(W_data, jn_cond_out$ULCI, lwd = 2, lty = 2, type = "l")
    abline(h = 0, untf = FALSE, lty = 3, lwd = 1)
    abline(v = jn[1], untf = FALSE, lty = 3,lwd = 1)
    abline(v = jn[2], untf = FALSE, lty = 3, lwd = 1)

  }

}
