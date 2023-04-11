#' LoadRData
#'
#' @return loaded R object
#'
#' @noRd
loadRData <- function( ### Loads an RData file, and returns it
                      ######################################################
                      fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
} ######################################################


#' create_conditions
#'
#' aim of this file is to keep track of new wrapper/workflow functions
#'
#' @noRd
create_conditions <- function(treatment = NULL,
                              strat = NULL,
                              survival = NULL) {
  conditions <- expand.grid(list(
    treatment = treatment,
    strat = strat,
    survival = survival
  ))
  return(conditions)
}


#' Creation of enrichment given CMS and condition
#'
#' performs hypergeometric tests
#'
#' @param tab table of subgroup 1 and subgroup 2
#' @param col_label name of subgroup 1
#' @param row_label name of subgroup 2
#' @param digits round of p-value to digit x
#'
#'
#' @return Returns enrichments in given treatment arms and their p-values\
#'
#' @noRd
enrichment <- function( ######################################################
                       tab, # table of counts for two covariates
                       col_label, # name of columns
                       row_label, # name of rows
                       digits # round of p-value to digit x
) {
  tab_test <- matrix(ncol = ncol(tab), nrow = nrow(tab))
  for (i in seq_len(nrow(tab_test))) {
    for (j in seq_len(ncol(tab_test))) {
      p_under <- phyper(
        q = tab[i, j],
        m = sum(tab[i, ]),
        n = sum(tab[-i, ]),
        k = sum(tab[, j]),
        lower.tail = TRUE # Under-representation
      )
      p_over <- phyper(
        q = tab[i, j] - 1,
        m = sum(tab[i, ]),
        n = sum(tab[-i, ]),
        k = sum(tab[, j]),
        lower.tail = FALSE # Over-representation
      )
      p <- which.min(c(p_under, p_over))
      if (p == 1) {
        p <- -p_under
      }
      if (p == 2) {
        p <- p_over
      }
      tab_test[i, j] <- p
    }
    row.names(tab_test) <- row.names(tab)
    colnames(tab_test) <- colnames(tab)
  }
  tab_test <- reshape2::melt(tab_test)
  colnames(tab_test) <- c(row_label, col_label, "value")
  tab_test$label <- abs(round(tab_test$value, digits = digits))
  tab_test$label <- unlist(lapply(tab_test$label, function(x) {
    if (x < 0.05) {
      paste0("p=", as.character(format.pval(x)), "*")
    } else {
      paste0("p=", as.character(format.pval(x)), "")
    }
  }))

  tab_test$Enrichment <- sign(tab_test$value) %>%
    as.character() %>%
    dplyr::recode("1" = "Over-enrichment", "-1" = "Under-enrichment")

  return(tab_test)
} ######################################################



#' Enrichment genomics
#'
#' Calculate enrichment of mutations in pre-defined subtypes
#'
#' @param bem A binary event matrix
#' @param clin A Binary Mutation Data Frame.
#' Row: Sample X1-Xn, Column: Gene of interest
#' @param sample_column Specify sample name
#' (First column name in Binary Mutation Data Frame)
#' @param min_mutants Specifies minimum amount of mutants per alteration
#' @param clin_alterations Clinical factor of interest
#' @param subtype Specify which (tumour) subtype to test
#'
#' @return A data frame of enriched somatic mutations in specified
#' subtypes and their adjusted p-values
#'
#' @details With a binary event matrix,
#' investigate enrichments of somatic mutations in subtypes
#'
#'
#' @noRd
enrichment_genomics <- function( ######################################################
                                bem, # binary event matrix
                                clin, # clinical data
                                sample_column = "sample", # which column identifies samples
                                min_mutants = 5, # how many mutants must be present at least per alteration
                                clin_alterations = NULL, # what clinical factor to add in bem
                                subtype = "primary.site" # which tumour subtype to test
) {
  # Check data types
  bem <- as.data.frame(bem) %>% mutate_all(as.character)
  clin <- as.data.frame(clin)
  intersecting_samples <- intersect(bem[, sample_column], clin[, sample_column])
  bem <- bem[bem[, sample_column] %in% intersecting_samples, ] %>%
    dplyr::select(-sample_column)
  clin <- clin[clin[, sample_column] %in% intersecting_samples, ] %>%
    dplyr::select(-sample_column)


  bem_cox_fig4 <- bem %>% mutate_all(as.numeric)

  sav_cox <- clin
  ### get only variants that more abundant than min_mutants
  bem_cox_fig4 <- bem_cox_fig4[, apply(
    bem_cox_fig4,
    2, function(x) sum(x)
  ) >= min_mutants]

  if (!any(is.null(clin_alterations))) {
    if (!all(apply(
      clin[, clin_alterations, drop = FALSE],
      2, function(x) length(levels(factor(x))) == 2
    ))) {
      stop("clin_alteration factor columns do not have exactly two levels.
           Please only choose two-level factor columns.")
    }
  }
  bem_cox_fig4 <- cbind(
    bem_cox_fig4,
    2 - sav_cox[, colnames(sav_cox) %in% clin_alterations, drop = FALSE] %>%
      mutate_all(factor) %>%
      mutate_all(as.numeric)
  )

  if (!subtype %in% colnames(sav_cox)) {
    stop("Subtype column is not in clinical data.")
  }

  sav_cox <- mutate_at(sav_cox, subtype, as.factor)

  tab_site <- list()
  for (x in colnames(bem_cox_fig4)) {
    tab2 <- table(bem_cox_fig4[, x] %>%
      dplyr::recode("1" = x, "0" = "wild type"), sav_cox[, subtype])

    itself <- paste0(subtype, levels(sav_cox[, subtype]))
    if (make.names(x) %in% make.names(itself)) {
      tab_site[[x]] <- NULL
    } else {
      tab_site[[x]] <- enrichment(tab2,
        col_label = "alteration", row_label = "status", 5
      )
    }
  }
  tab_site <- do.call(rbind, tab_site)
  tab_site <- tab_site[tab_site$status != "wild type", ]
  tab_site$Enrichment <- sign(tab_site$value) %>%
    as.character() %>%
    dplyr::recode("1" = "Over-enrichment", "-1" = "Under-enrichment")
  tab_site$abundance <- lapply(
    tab_site$status,
    function(x) sum(na.omit(bem_cox_fig4[, as.character(x)]))
  ) %>% unlist()

  #
  tab_site$fdr[(tab_site$Enrichment == "Over-enrichment")] <- p.adjust(abs(
    tab_site$value[(tab_site$Enrichment == "Over-enrichment")]
  ),
  method = "BH"
  )
  tab_site$fdr[(tab_site$Enrichment == "Under-enrichment")] <- p.adjust(abs(
    tab_site$value[(tab_site$Enrichment == "Under-enrichment")]
  ),
  method = "BH"
  )
  return(tab_site)
} #################################################



#' Make Mutex DF
#'
#' prepare data-frame for running the mutex algorithm
#'
#' @param bem binary event matrix, _AMP and _LOSS must be included in
#' column names of alterations
#'
#' @return dataframe including somatic mutations and their type for mutex input
#'
#' @noRd
make_mutex_df <- function( ######################################################
                          bem = NULL # binary event matrix, _AMP and _LOSS must be included in column names of alterations
) {
  names <- lapply(colnames(bem), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()
  types <- lapply(colnames(bem), function(x) strsplit(x, "_")[[1]][2]) %>% unlist()
  doubled <- names[which(unlist(lapply(
    names,
    function(x) length(which(x == names))
  )) == 2)]
  for (i in seq_len(ncol(bem))) {
    bem_tmp <- bem[, i]
    if (types[i] == "AMP") {
      bem_tmp[bem_tmp == "1"] <- "2"
    }
    if (types[i] == "LOSS") {
      bem_tmp[bem_tmp == "1"] <- "3"
    }
    bem[, i] <- bem_tmp
  }
  bem <- as.data.frame(t(as.matrix(bem)))
  bem$X <- names

  # Summarize
  # -1 because of bem$X is adding one above
  res <- as.data.frame(matrix(ncol = ncol(bem) - 1, nrow = 0))
  colnames(res) <- colnames(bem)[seq_len((ncol(bem) - 1))]
  for (i in seq_len(length(names))) {
    name <- names[i]
    bem_tmp <- bem[bem$X %in% name, ]
    bem_tmp <- bem_tmp[seq_len((length(bem_tmp) - 1))]
    # types <-  lapply(row.names(bem_tmp),
    # function(x) strsplit(x, "_")[[1]][2])%>% unlist

    new <- rep(NA, ncol(bem_tmp))
    new[apply(bem_tmp, 2, function(x) all(x == "0"))] <- 0
    new[apply(bem_tmp, 2, function(x) ("1" %in% x))] <- 1
    new[apply(bem_tmp, 2, function(x) ("2" %in% x))] <- 2
    new[apply(bem_tmp, 2, function(x) ("3" %in% x))] <- 3
    new[apply(bem_tmp, 2, function(x) ("3" %in% x) & ("1" %in% x))] <- 5
    new[apply(bem_tmp, 2, function(x) ("2" %in% x) & ("1" %in% x))] <- 4

    # message(i)

    res[i, ] <- new
  }

  res$gene <- names
  res <- res[!duplicated(res$gene), ]
  row.names(res) <- res$gene
  res <- res[, seq_len((ncol(res) - 1))]

  return(res)
} #################################################



#' Plot Cox
#'
#' Fit cox/logistic regression models for different readouts
#'
#' @param features feature dataframe
#' @param response readout dataframe with potential confounders
#' @param feature_name Feature_name like KRAS, a molecular marker
#' @param strat What feature to filter for in response
#' @param treatment (Double) treatment e.g. "FOLFIRI+Bevacizumab",
#' "FOLFIRI+Cetuximab"
#' @param survival Given in plot example, e.g. "OS", "PFS"
#' @param covariates Array of covariates
#' @param N Minimal number of cell lines per condition
#' @param rna_asso Either True or False
#' @param compare_treatments True if one wants to compare treatments
#' between each other
#'
#' @return list including model object and other statistics
#'
#' @noRd
plot_cox <- function( ######################################################
                     treatment, # treatment name
                     features = NULL, # feature vectors, $bem_cox
                     response = NULL, # response dataframe with survival and covariates, $sav_cox
                     feature_name = NULL, # a molecular marker
                     strat = NULL, # what feature to filter for in response
                     survival = "OS", # also PFS
                     covariates = c("1"), # array of covariates
                     N = N, # minimal number of cell lines per condition
                     rna_asso = NULL,
                     compare_treatments = FALSE) {
  response <- response[order(row.names(response)), ]
  features <- features[order(row.names(features)), ]
  if (all(row.names(features) == row.names(response)) == FALSE) {
    stop("Rownames of features and response do not match !")
  }
  # Filter for treatment
  if (compare_treatments == FALSE) {
    features <- features[response$treatment %in% treatment, ]
    response <- response[response$treatment %in% treatment, ]
  }

  # Standardize treatment names
  tmp <- levels(response$treatment)
  response$treatment <- as.character(response$treatment)
  response$treatment[response$treatment == tmp[1]] <- "t1"
  response$treatment[response$treatment == tmp[2]] <- "t2"
  response$treatment <- factor(response$treatment)

  # Filter for strat variable
  strats <- response[, strat]
  strats_levels <- levels(factor(strats))
  howmany_strats <- length(strats_levels)
  if (nrow(response) == 0) {
    stop("Stratification for stat variable failed, not present in dataset ?")
  }
  List <- list()

  if (compare_treatments) {
    if (survival == "OS") {
      survival <- "OSint"
    }
    if (survival == "PFS") {
      survival <- "PFSint"
    }
    if (survival == "ORR") {
      survival <- "ORRint"
    }
    if (survival == "DFS") {
      survival <- "DFSint"
    }
    if (survival == "RFS") {
      survival <- "RFSint"
    }
  }

  covariates <- paste("+", paste(covariates, collapse = "+"), sep = "")
  for (level in strats_levels) {
    logical <- strats == level
    # remove NA observations from strat observations
    logical[is.na(logical)] <- FALSE
    temp_response <- response[logical, ]
    temp_features <- features[logical, ]
    if (all(row.names(temp_features) == row.names(temp_response)) == FALSE) {
      stop("Rownames do not match in between strat levels!")
    }

    if (feature_name %in% names(temp_features) == FALSE) {
      stop("feature_name does not match features!")
    }
    temp_features <- temp_features[, feature_name]
    temp_response$biomarker <- temp_features
    if (survival %in% c("OS", "PFS", "DFS", "RFS")) {
      n1 <- length(which(temp_features == "1"))
      n2 <- length(which(temp_features == "0"))
    }
    if (survival == "ORR") {
      n1 <- length(which(temp_features[!is.na(temp_response[, "ORR"])] == "1"))
      n2 <- length(which(temp_features[!is.na(temp_response[, "ORR"])] == "0"))
    }

    # if interaction analysis
    if (compare_treatments) {
      temp_features <- factor(temp_features)
      levels(temp_features) <- c("1", "0")
      if (survival %in% c("OSint", "PFSint", "DFSint", "RFSint")) {
        tab <- table(temp_features, temp_response$treatment)
        n1 <- min(tab["1", ])
        n2 <- min(tab["0", ])
      }
      if (survival == "ORRint") {
        tab <- table(
          temp_features[!is.na(temp_response[, "ORR"])],
          temp_response$treatment[!is.na(temp_response[, "ORR"])]
        )
        n1 <- min(tab["1", ])
        n2 <- min(tab["0", ])
      }
    }
    temp_response$covariates <- rep(1, length(temp_response$biomarker))

    # if not minimal amount of mutants are present in one of the investigated populations
    if (!(n1 < N | n2 < N) | rna_asso) {
      if (survival == "OS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(OS.months, OS.event) ~ biomarker", covariates,
          sep = ""
        )), data = temp_response) # INT
        km_fit <- survival::survfit(as.formula(paste(
          "survival::Surv(OS.months, OS.event) ~ biomarker", "",
          sep = ""
        )),
        data = temp_response
        )
      }
      if (survival == "OSint") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(OS.months, OS.event) ~ biomarker + biomarker*treatment",
          covariates,
          sep = ""
        )), data = temp_response) # INT
        km_fit <- NULL
      }
      if (survival == "PFS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(PFS.months, PFS.event) ~ biomarker",
          covariates,
          sep = ""
        )), data = temp_response)
        km_fit <- survival::survfit(as.formula(paste(
          "survival::Surv(PFS.months, PFS.event) ~ biomarker", "",
          sep = ""
        )),
        data = temp_response
        )
      }
      if (survival == "PFSint") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(PFS.months, PFS.event) ~ biomarker + biomarker*treatment",
          covariates,
          sep = ""
        )), data = temp_response)
        km_fit <- NULL
      }
      if (survival == "ORR") {
        cox <- glm(as.formula(paste(
          "ORR ~ biomarker + biomarker", covariates,
          sep = ""
        )),
        data = temp_response, family = binomial(link = "logit")
        )
        km_fit <- list(as.data.frame(
          table(marker = temp_response$biomarker, resp = temp_response[, "ORR"])
        ))
      }
      if (survival == "ORRint") {
        cox <- glm(as.formula(paste(
          "ORR ~ biomarker + biomarker*treatment", covariates,
          sep = ""
        )),
        data = temp_response, family = binomial(link = "logit")
        )
        km_fit <- NULL
      }


      if (survival == "DFS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(DFS.months, DFS.event) ~ biomarker",
          covariates,
          sep = ""
        )), data = temp_response)
        km_fit <- survival::survfit(as.formula(paste(
          "survival::Surv(DFS.months, DFS.event) ~ biomarker", "",
          sep = ""
        )), data = temp_response)
      }
      if (survival == "DFSint") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(DFS.months, DFS.event) ~ biomarker + biomarker*treatment",
          covariates,
          sep = ""
        )), data = temp_response)
        # km_fit <- survival::survfit(as.formula(paste("survival::Surv(PFS.months, PFS.event) ~ biomarker","",sep="")), data = temp_response)
        km_fit <- NULL
      }
      if (survival == "RFS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(RFS.months, RFS.event) ~ biomarker",
          covariates,
          sep = ""
        )), data = temp_response)
        km_fit <- survival::survfit(as.formula(paste(
          "survival::Surv(RFS.months, RFS.event) ~ biomarker", "",
          sep = ""
        )),
        data = temp_response
        )
      }
      if (survival == "RFSint") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(RFS.months, RFS.event) ~ biomarker + biomarker*treatment",
          covariates,
          sep = ""
        )), data = temp_response)
        # km_fit <- survival::survfit(as.formula(paste("survival::Surv(PFS.months, PFS.event) ~ biomarker","",sep="")), data = temp_response)
        km_fit <- NULL
      }

      if (!survival %in% c(
        "OS", "OSint", "PFS", "PFSint", "ORR",
        "ORRint", "DFS", "DFSint", "RFS", "RFSint"
      )) {
        cox <- NULL
        km_fit <- NULL
      }
    } else {
      cox <- NULL
      km_fit <- NULL
    }

    List[[level]] <- list(
      cox = cox,
      fit = km_fit,
      strats_levels = level,
      n1 = n1,
      n2 = n2
    )
  }

  return(List)
} ######################################################



#' Calculate Association Cox
#'
#' wrapper for plot_cox to calculate cox/logistic regression
#' with different putative biomarkers in different subtypes
#'
#' @param plot_volcano logical set to false
#' @param features feature vectors, $bem_cox
#' @param response response dataframe with survival and covariates, $sav_cox
#' @param strat what feature to filter for in response
#' @param treatment Specify treatment, currently: "FOLFIRI+Bevacizumab"
#' and "FOLFIRI+Cetuximab"
#' @param survival readout
#' @param covariates confounding factors
#' @param N minimal number of tumors per condition
#' @param rna_asso (deprecated)
#' @param compare_treatments logical if treatments will be compared or not
#' @param min_redistribution how many patients get different status when using
#' MEs
#'
#' @return results of association tests for specific treatments
#'
#' @noRd
calc_asso_cox <- function( ######################################################
                          treatment, # treatment
                          plot_volcano = FALSE,
                          features = NULL, # feature vectors, $bem_cox
                          response = NULL, # response dataframe with survival and covariates, $sav_cox
                          strat = NULL, # what feature to filter for in response
                          survival = "OS", # also PFS
                          covariates = c("1"), # array of covariates
                          N = 10, # minimal number of tumors per condition
                          rna_asso = FALSE, # should be FALSE, lost its function
                          compare_treatments = FALSE, #
                          min_redistribution = 4 # how many patients get different status when using MEs
) {
  isFunctionCalled <- FALSE
  features[] <- mutate_all(.tbl = features, .funs = as.character)
  int <- intersect(features$sample, response$sample)
  features <- features[features$sample %in% int, ] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample")
  response <- response[response$sample %in% int, ] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample")
  if (survival %in% c("OS", "PFS", "ORR", "DFS", "RFS")) {
    prob <- "Pr(>|z|)"
  } else {
    prob <- "Pr(>|t|)"
  }

  for (i in seq_along(colnames(features))) {
    # message(paste("Start run #",as.character(i),sep=""))
    NAME <- colnames(features)[i]

    # FIT
    model <- plot_cox(
      features = features, response = response,
      feature_name = NAME, strat = as.character(strat),
      treatment = treatment, survival = survival,
      covariates = covariates, N = N, rna_asso = rna_asso,
      compare_treatments = compare_treatments
    )
    #####

    if (isFunctionCalled == FALSE & compare_treatments == FALSE) {
      asso_cox <- lapply(
        seq_along(model),
        function(x) data.frame(matrix(nrow = 0, ncol = 8))
      )
      for (x in seq_along(model)) {
        colnames(asso_cox[[x]]) <- c(
          "effectsize", "P.Value", "cfe",
          "strat", "npos", "nneg", "model", "medians"
        )
      }

      isFunctionCalled <- TRUE
    }
    if (isFunctionCalled == FALSE & compare_treatments == TRUE) {
      asso_cox <- lapply(
        seq_along(model),
        function(x) data.frame(matrix(nrow = 0, ncol = 10))
      )
      for (x in seq_along(model)) {
        colnames(asso_cox[[x]]) <- c(
          "effectsize", "P.Value", "cfe",
          "strat", "npos", "nneg", "pInt",
          "effInt", "model", "medians"
        )
      }

      isFunctionCalled <- TRUE
    }

    for (x in seq_along(model)) {
      if (!(model[[x]]$n1 < N | model[[x]]$n2 < N) | rna_asso) {
        if (rna_asso) {
          if (compare_treatments == FALSE) {
            asso_cox[[x]][i, ] <- list(
              model[[x]]$cox$coefficients["biomarker"],
              (summary(model[[x]]$cox))$coefficients["biomarker", prob],
              NAME,
              model[[x]]$strats_levels,
              model[[x]]$n1,
              model[[x]]$n2,
              list(model[[x]]),
              list(survminer::surv_median(model[[x]]$fit))
            )
          } else {
            asso_cox[[x]][i, ] <- list(model[[x]]$cox$coefficients["biomarker1"],
              (summary(model[[x]]$cox))$coefficients["biomarker1", prob],
              NAME,
              model[[x]]$strats_levels,
              model[[x]]$n1,
              model[[x]]$n2,
              pInt = NA,
              effInt = NA,
              list(model[[x]]),
              list(survminer::surv_median(model[[x]]$fit))
            )
          }
        } else {
          if (compare_treatments == FALSE) {
            asso_cox[[x]][i, ] <- list(
              model[[x]]$cox$coefficients["biomarker1"],
              (summary(model[[x]]$cox))$coefficients["biomarker1", prob],
              NAME,
              model[[x]]$strats_levels,
              model[[x]]$n1,
              model[[x]]$n2,
              list(model[[x]]),
              tryCatch(list(survminer::surv_median(model[[x]]$fit)),
                error = function(e) list(model[[x]]$fit)
              )
            )
          } else {
            asso_cox[[x]][i, ] <- list(model[[x]]$cox$coefficients["biomarker1"],
              (summary(model[[x]]$cox))$coefficients["biomarker1", prob],
              NAME,
              model[[x]]$strats_levels,
              model[[x]]$n1,
              model[[x]]$n2,
              pInt = (summary(model[[x]]$cox))$coefficients[
                "biomarker1:treatmentt2",
                "Pr(>|z|)"
              ], # INT
              effInt = model[[x]]$cox$coefficients["biomarker1:treatmentt2"], # INT
              list(model[[x]]),
              tryCatch(list(survminer::surv_median(model[[x]]$fit)),
                error = function(e) list(model[[x]]$fit)
              )
            )
          }
        }
      } else {
        if (compare_treatments == FALSE) {
          asso_cox[[x]][i, ] <- c(
            NA, NA, NAME, model[[x]]$strats_levels,
            model[[x]]$n1, model[[x]]$n2, NA, NA
          )
        } else {
          asso_cox[[x]][i, ] <- c(NA, NA, NAME, model[[x]]$strats_levels,
            model[[x]]$n1, model[[x]]$n2,
            pInt = NA,
            effInt = NA, model = NA, medians = NA
          )
        }
      }
    }
  }


  for (x in seq_along(model)) {
    asso_cox[[x]]$fdr <- p.adjust(asso_cox[[x]]$P.Value, method = "BH")
  }


  for (x in seq_len(length(asso_cox))) {
    asso_cox[[x]]$gene <- colnames(features)
    asso_cox[[x]]$effectsize <- as.numeric(asso_cox[[x]]$effectsize)
    asso_cox[[x]]$P.Value <- as.numeric(asso_cox[[x]]$P.Value)
    asso_cox[[x]]$fdr <- as.numeric(asso_cox[[x]]$fdr)
    asso_cox[[x]]$npos <- as.numeric(asso_cox[[x]]$npos)
    asso_cox[[x]]$nneg <- as.numeric(asso_cox[[x]]$nneg)
    # confidence score is the minimum of mutants or wild types
    asso_cox[[x]]$confidence <- unlist(lapply(
      seq_len(nrow(asso_cox[[x]])),
      function(y) min(asso_cox[[x]]$npos[y], asso_cox[[x]]$nneg[y])
    ))
    asso_cox[[x]]$survival <- rep(survival, nrow(asso_cox[[x]]))
    asso_cox[[x]]$treatment <- rep(treatment, nrow(asso_cox[[x]]))

    ### Added function for identifying if ME module, FILTERING
    min_redistributed <- lapply(seq_len(nrow(asso_cox[[x]])), function(y) {
      redistributed(
        alteration = asso_cox[[x]][y, ]$gene, # column name of mutations_df
        # data_mutations dataframe
        mutations = features %>% tibble::rownames_to_column("sample"),
        clinical = response %>% tibble::rownames_to_column("sample"),
        subtype_column = strat,
        subtype = asso_cox[[x]][y, ]$strat, # which subpopulation
        treatment = treatment # what subtype
      )
    })
    names(min_redistributed) <- paste(
      asso_cox[[x]]$gene,
      names(min_redistributed)
    )

    asso_cox[[x]] <- (asso_cox[[x]])[unlist(min_redistributed)
    >= min_redistribution, ]
    ###############################################

    asso_cox[[x]]$fdr <- p.adjust(asso_cox[[x]]$P.Value, method = "BH")
  }

  if (compare_treatments) {
    for (x in seq_len(length(asso_cox))) {
      asso_cox[[x]]$pInt <- as.numeric(asso_cox[[x]]$pInt)
      asso_cox[[x]]$effInt <- as.numeric(asso_cox[[x]]$effInt)
    }
  }

  asso_cox <- do.call(rbind, asso_cox)

  # filter non-tested
  asso_cox <- asso_cox[!is.na(asso_cox$P.Value), ]
  asso_cox$fdr_new <- p.adjust(asso_cox$P.Value, method = "BH")
  if (compare_treatments) {
    asso_cox$fdr_new_int <- p.adjust(asso_cox$pInt, method = "BH")
  }

  # change ORR sign
  asso_cox$effectsize[asso_cox$survival == "ORR"] <-
    -asso_cox$effectsize[asso_cox$survival == "ORR"]

  # apply new FDR for not comparing treatments
  if (!compare_treatments) {
    asso_cox <- asso_cox %>%
      dplyr::rename(fdr = fdr_new, fdr_old = fdr) %>%
      dplyr::select(-fdr_old)
  }

  # apply new FDR for not comparing treatments
  if (compare_treatments) {
    asso_cox <- asso_cox %>%
      dplyr::rename(fdr_int = fdr_new_int) %>%
      dplyr::select(-c(effectsize, P.Value, fdr, fdr_new))
  }

  return(asso_cox)
} ######################################################



#' Add me
#'
#' helper function adding mutually exclusive modules to binary event matrix
#' of somatic mutations
#'
#' @param data -
#' @param names -
#' @param name -
#' @param logic "or","and" argument for grouping mutations with 'or' or 'and'
#' logic
#'
#' @return -
#'
#' @noRd
add_me <- function( ######################################################
                   data,
                   names,
                   name,
                   logic = "or") {
  names <- names[names != ""]
  names1 <- paste0(names, "_SV")
  names2 <- paste0(names, "_AMP")
  names3 <- paste0(names, "_DEL")

  names <- c(names1, names2, names3)
  names <- names[names %in% colnames(data)]
  message(names, collapse = "_or_")
  alteration <- data[, names, drop = FALSE]
  if (logic == "or") {
    alteration <- as.character(as.numeric(apply(
      alteration, 1,
      function(x) any(as.logical(as.numeric(x)))
    )))
  }
  if (logic == "not") {
    alteration <- as.character(as.numeric(apply(
      alteration, 1,
      function(x) {
        as.logical(as.numeric(x[1])) &
          !as.logical(as.numeric(x[2]))
      }
    )))
  }
  data <- cbind(data, alteration)
  colnames(data)[ncol(data)] <- paste(names, collapse = "_or_")
  return(data)
}



#' Make Forest
#'
#' get confidence intervals and effect sizes
#'
#' @param df dataframe containing association test results
#'
#' @return dataframe with added columns containing additional test statistics
#'
#' @noRd
make_forest <- function(df) {
  testFunction <- function(x, y) {
    family <- tryCatch(x$family$family, error = function(e) NULL)
    if (is.null(family)) {
      family <- "none"
    }
    if (y == 1) {
      tmp <- tryCatch(x$conf.int["biomarker1", "exp(coef)"], error = function(e) NULL)
      tmp <- tryCatch(x$conf.int[1, "exp(coef)"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients[2, "Estimate"]), error = function(e) NULL)
      }
    }
    if (y == 2) {
      tmp <- tryCatch(x$conf.int["biomarker1", "upper .95"], error = function(e) NULL)
      tmp <- tryCatch(x$conf.int[1, "upper .95"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients[2, "Estimate"] -
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual)),
        error = function(e) NULL
        )
      }
    }
    if (y == 3) {
      tmp <- tryCatch(x$conf.int["biomarker1", "lower .95"], error = function(e) NULL)
      tmp <- tryCatch(x$conf.int[1, "lower .95"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients[2, "Estimate"] +
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual)),
        error = function(e) NULL
        )
      }
    }
    if (is.null(tmp)) {
      tmp <- NA
    }
    return(tmp)
  }
  df$hazard <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 1)
  ) %>% unlist()
  df$lower <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 2)
  ) %>% unlist()
  df$upper <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 3)
  ) %>% unlist()
  return(df)
}



#' Make Forest Interaction
#'
#' get confidence intervals and effect sizes
#'
#' @param df dataframe containing association test results (interactions)
#'
#' @return dataframe with added columns containing additional test statistics
#' (interactions)
#'
#' @noRd
make_forest_interaction <- function(df) {
  testFunction <- function(x, y) {
    family <- tryCatch(x$family$family, error = function(e) NULL)
    if (is.null(family)) {
      family <- "none"
    }
    if (y == 1) {
      tmp <- tryCatch(x$conf.int["biomarker1:treatmentt2", "exp(coef)"],
        error = function(e) NULL
      )
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients["biomarker1:treatmentt2", "Estimate"]),
          error = function(e) NULL
        )
      } # for ORR
    }
    if (y == 2) {
      tmp <- tryCatch(x$conf.int["biomarker1:treatmentt2", "upper .95"],
        error = function(e) NULL
      )
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients["biomarker1:treatmentt2", "Estimate"] -
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual)),
        error = function(e) NULL
        )
      }
    }
    if (y == 3) {
      tmp <- tryCatch(x$conf.int["biomarker1:treatmentt2", "lower .95"],
        error = function(e) NULL
      )
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients["biomarker1:treatmentt2", "Estimate"] +
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual)),
        error = function(e) NULL
        )
      }
    }
    if (is.null(tmp)) {
      tmp <- NA
    }
    return(tmp)
  }
  df$hazard <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 1)
  ) %>% unlist()
  df$lower <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 2)
  ) %>% unlist()
  df$upper <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 3)
  ) %>% unlist()
  return(df)
}


make_forest_comp <- function(df) {
  testFunction <- function(x, y) {
    family <- tryCatch(x$family$family, error = function(e) NULL)
    if (is.null(family)) {
      family <- "none"
    }
    if (y == 1) {
      tmp <- tryCatch(x$conf.int["biomarker1", "exp(coef)"],
        error = function(e) NULL
      )
      tmp <- tryCatch(x$conf.int[1, "exp(coef)"],
        error = function(e) NULL
      )
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients[2, "Estimate"]),
          error = function(e) NULL
        )
      }
    }
    if (y == 2) {
      tmp <- tryCatch(x$conf.int["biomarker1", "upper .95"], error = function(e) NULL)
      tmp <- tryCatch(x$conf.int[1, "upper .95"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients[2, "Estimate"] +
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual)),
        error = function(e) NULL
        )
      }
    }
    if (y == 3) {
      tmp <- tryCatch(x$conf.int["biomarker1", "lower .95"],
        error = function(e) NULL
      )
      tmp <- tryCatch(x$conf.int[1, "lower .95"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(exp(-x$coefficients[2, "Estimate"] -
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual)),
        error = function(e) NULL
        )
      }
    }
    if (is.null(tmp)) {
      tmp <- NA
    }
    return(tmp)
  }
  df$hazard <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]])), y = 1)
  ) %>% unlist()
  df$lower <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]])), y = 2)
  ) %>% unlist()
  df$upper <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]])), y = 3)
  ) %>% unlist()
  return(df)
}



#' Impute Median
#'
#' helper function to median-imput features
#'
#' @param x numeric vector
#'
#' @return numeric vector
#'
#' @noRd
impute_median <- function(x) {
  ind_na <- is.na(x)
  x[ind_na] <- median(x[!ind_na])
  as.numeric(x)
}



#' Fit Model via cross validation
#'
#' fits multiple regression model
#'
#' @param covariates feature matrix
#' @param fold how many folds for cross-validation
#' @param resamples how many times is the cross-validation repeated
#' @param alpha elastic net hyperparameter
#' @param method (cox/logistic) regression
#' @param continous.response logical if response is continuous
#'
#' @return list of model statistics and feature weights
#'
#' @noRd
fit_model <- function(covariates,
                      fold = 10,
                      resamples = 10,
                      alpha = 0,
                      method = "cox",
                      continous.response = FALSE) {
  response <- covariates$survival
  survival_var <- covariates$survival_var
  covariates <- covariates$covariates

  res <- list(
    models_null = NULL,
    models_null_performance = NULL,
    models = NULL,
    models_performance = NULL,
    feature_importance = NULL,
    survival_var = survival_var
  )

  models_list <- list()
  for (i in seq_len(resamples)) {
    meth <- if (survival_var %in% c("OS", "PFS", "DFS", "RFS")) {
      method
    } else {
      "glm_bin"
    }
    if (continous.response != FALSE) {
      meth <- "glm"
    }
    models_list[[i]] <- run_pipeline_benchmark(
      feature_path = covariates,
      response_path = response,
      kfold = fold,
      method = meth,
      hyperparam = if (method == "cox") {
        c("alpha" = alpha)
      } else {
        list(c(5), c(100))
      },
      cvglm = TRUE,
      returnFit = TRUE,
      cvseed = i,
      null.model = TRUE
    )
  }
  lapply(models_list, function(x) x[["score"]])
  score_null <- lapply(models_list, function(x) x[["score"]][, 1]) %>%
    unlist()
  score_null[is.na(score_null)] <- 0.5
  res$models_null_performance <- list(mean = mean(score_null), sd = sd(score_null))
  res$models_null <- models_list

  models_list <- list()
  for (i in seq_len(resamples)) {
    models_list[[i]] <- run_pipeline_benchmark(
      feature_path = covariates,
      response_path = response,
      kfold = fold,
      method = meth,
      hyperparam = if (method == "cox") {
        c("alpha" = alpha)
      } else {
        list(c(5), c(100))
      },
      cvglm = TRUE,
      returnFit = TRUE,
      cvseed = i
    )
  }

  score <- lapply(models_list, function(x) x[["score"]][, 1]) %>%
    unlist()
  score[is.na(score)] <- 0.5

  # bayes factor
  bayes <- sum(score > score_null) / length(score_null) /
    (1 - (sum(score > score_null) / length(score_null)))
  res$models_performance <- list(mean = mean(score), sd = sd(score), bayes = bayes)
  res$models <- models_list

  for (i in seq_len(resamples)) {
    for (j in seq_len(fold)) {
      fit <- models_list[[i]]$param[j, 1][[1]][[1]]
      beta <- coef(fit, s = "lambda.min")
      beta <- beta[-1, ]
      if (i == 1 & j == 1) {
        beta_df <- as.data.frame(beta)
      }
      beta_df <- cbind(beta_df, beta)
    }
  }
  beta_df$mean <- apply(beta_df, 1, mean)
  beta_df$sd <- apply(beta_df, 1, function(x) sd(x))
  beta_df_red <- beta_df[abs(beta_df$mean) - beta_df$sd > 0, ]

  res$feature_importance <- list(all = beta_df, red = beta_df_red)
  return(res)
}


#################################################
#' Compare treatment arms
#'
#' helper function for comparing treatment arms
#'
#' @param asso_cox_bev association tests for treatment arm 1
#' @param asso_cox_cet association tests for treatment arm 2
#' @param strat_variable tumour subtype variable
#' @param name name of parent frame (name)
#' @param cms name of parent frame (cms)
#' @param sav_cox_fig5 clinical dataframe
#' @param bem_cox_fig4_me mutational dataframe
#' @param add_row add row (binary)
#' @param readout readout (survival/binary)
#' @param fdr_threshold (deprecated)
#' @param verbose verbose
#' @param sample_column identifier for individuals
#' @param covariates confounding factors
#'
#' @return dataframe containing statistics for treatment comparisons
#'
#' @noRd
compare_treatments <- function(asso_cox_bev = parent.frame()$asso_cox_bev,
                               asso_cox_cet = parent.frame()$asso_cox_cet,
                               strat_variable = parent.frame()$strat_variable,
                               name = parent.frame()$name,
                               cms = parent.frame()$cms,
                               sav_cox_fig5 = parent.frame()$sav_cox_fig5,
                               bem_cox_fig4_me = parent.frame()$bem_cox_fig4_me,
                               add_row = FALSE,
                               readout = "OS",
                               fdr_threshold = 1, # deprecated
                               verbose = FALSE,
                               sample_column = "sample",
                               covariates = NULL) {
  # Check data types
  bem_cox_fig4_me <- as.data.frame(bem_cox_fig4_me)
  sav_cox_fig5 <- as.data.frame(sav_cox_fig5)

  intersecting_samples <- intersect(
    bem_cox_fig4_me[, sample_column],
    sav_cox_fig5[, sample_column]
  )
  bem_cox_fig4_me <- bem_cox_fig4_me[bem_cox_fig4_me[, sample_column] %in% intersecting_samples, ]
  sav_cox_fig5 <- sav_cox_fig5[sav_cox_fig5[, sample_column] %in% intersecting_samples, ]



  asso_cox <- bind_rows(
    asso_cox_bev %>% mutate(treatment = rep("t2_me", nrow(asso_cox_bev))),
    asso_cox_cet %>% mutate(treatment = rep("t1_me", nrow(asso_cox_cet)))
  )

  sign <- (as.numeric(asso_cox$fdr) < !is.na(asso_cox$fdr))
  message("number of significants: ", length(which(sign)))
  if (length(which(sign)) == 0) {
    return(NULL)
  }
  names(sign) <- rep("", length(sign))
  names(sign)[sign(asso_cox$effectsize) < 0 & sign] <- "sensitive"
  names(sign)[sign(asso_cox$effectsize) > 0 & sign] <- "resistant"
  asso_cox <- asso_cox[sign, ] %>% mutate(direction = names(sign)[sign])
  asso_cox$cfe <- lapply(
    asso_cox$cfe,
    function(x) {
      gsub(
        "primary.siteright-sided",
        "not-primary.siteleft-sided", x
      )
    }
  ) %>% unlist()

  # add artificial row
  if (!is.logical(add_row)) {
    asso_cox <- rbind(
      asso_cox,
      c(
        10, 0, add_row, name, 15, 15, NA, 0,
        add_row, 15, "t1_me", "sensitive"
      )
    )
  }

  resistant.cohort <- FALSE

  res <- list()
  for (i in seq_len(nrow(asso_cox))) {
    if (verbose) {
      message(i)
    }
    biomarker <- asso_cox$cfe[i]
    strat <- cms
    if (asso_cox$direction[i] == "sensitive") {
      if (resistant.cohort) {
        direction <- 0
      } else {
        direction <- 1
      }
    }
    if (asso_cox$direction[i] == "resistant") {
      if (resistant.cohort) {
        direction <- 1
      } else {
        direction <- 0
      }
    }
    if (asso_cox$treatment[i] %in% c("t2_me", "t1_me")) {
      sav_cox_strat <- sav_cox_fig5[(bem_cox_fig4_me[, biomarker] == direction) &
        (sav_cox_fig5[, strat_variable] == strat), ]
      bem_cox_strat <- bem_cox_fig4_me[(bem_cox_fig4_me[, biomarker] == direction) &
        (sav_cox_fig5[, strat_variable] == strat), ]
    }
    if (asso_cox$treatment[i] %in% c("t2_clin", "t1_clin")) {
      sav_cox_strat <- clin[clin[, biomarker] == direction, ]
      bem_cox_strat <- clin[clin[, biomarker] == direction, ]
    }
    sav_cox_strat$treatment_strat <- sav_cox_strat$treatment

    # Standardize treatment names
    tmp <- levels(sav_cox_strat$treatment_strat)
    sav_cox_strat$treatment_strat <- as.character(sav_cox_strat$treatment_strat)
    sav_cox_strat$treatment_strat[sav_cox_strat$treatment_strat == tmp[1]] <- "t1"
    sav_cox_strat$treatment_strat[sav_cox_strat$treatment_strat == tmp[2]] <- "t2"
    sav_cox_strat$treatment_strat <- factor(sav_cox_strat$treatment_strat)

    if (nrow(sav_cox_strat) > 0) {
      sav_cox_strat$treatment <- "yes"
    }
    m1 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = ""
    ), error = function(e) NULL)

    m2 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "PFS"
    ), error = function(e) NULL)

    m3 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "ORR"
    ), error = function(e) NA)


    m4 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "DFS"
    ), error = function(e) NULL)


    m5 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "RFS"
    ), error = function(e) NULL)

    smodel <- NULL
    if (readout == "OS") {
      smodel <- list(m1)
    }
    if (readout == "PFS") {
      smodel <- list(m2)
    }
    if (readout == "ORR") {
      smodel <- list(m3)
    }
    if (readout == "DFS") {
      smodel <- list(m4)
      m1 <- m4
    }
    if (readout == "RFS") {
      smodel <- list(m5)
      m1 <- m5
    }


    res[[i]] <- tryCatch(list(
      # biomarker p
      p = (summary(m1))$coefficients[paste0("treatment_stratt2", ""), "Pr(>|z|)"],
      # biomarker eff
      eff = (summary(m1))$coefficients["treatment_stratt2", "coef"],
      confidence = if (readout == "ORR") {
        min(table(sav_cox_strat$treatment_strat[!is.na(sav_cox_strat$ORR)]))
      } else {
        min(table(sav_cox_strat$treatment_strat))
      },
      treatment_int_p = NA,
      treatmentp = (summary(m1))$coefficients["treatment_stratt2", "Pr(>|z|)"],
      better = NA,
      alteration = biomarker,
      matrix = list(as.data.frame(sav_cox_strat)),
      treatment_info = asso_cox$treatment[i],
      direction = asso_cox$direction[i],
      strat_variable = name,
      model = smodel,
      mutantorwildtype = direction,
      p_pfs = tryCatch((summary(m2))$coefficients[paste0(
        "treatment_stratt2",
        ""
      ), "Pr(>|z|)"], error = function(e) NA), # biomarker p
      eff_pfs = tryCatch((summary(m2))$coefficients[
        "treatment_stratt2",
        "coef"
      ], error = function(e) NA), # biomarker eff
      p_orrs = tryCatch(coef(summary(m3))[2, 4], error = function(e) NA),
      eff_orr = tryCatch(coef(summary(m3))[2, 1], error = function(e) NA)
    ), error = function(e) NULL)
  }

  if (length(res) == 1) {
    res <- do.call(bind_rows, list(res)) %>% mutate_at(c(
      "p", "eff",
      "confidence", "p_pfs", "eff_pfs", "p_orrs", "eff_orr"
    ), as.numeric)
  } else {
    res <- do.call(bind_rows, res) %>% mutate_at(c(
      "p", "eff",
      "confidence", "p_pfs", "eff_pfs", "p_orrs", "eff_orr"
    ), as.numeric)
  }
  res$treatment_info <- unlist(lapply(
    res$treatment_info,
    function(x) paste0(x, "_sens")
  ))
  res2 <- res # save object

  resistant.cohort <- TRUE
  res <- list()
  for (i in seq_len(nrow(asso_cox))) {
    if (verbose) {
      message(i)
    }
    biomarker <- asso_cox$cfe[i]
    strat <- cms
    if (asso_cox$direction[i] == "sensitive") {
      if (resistant.cohort) {
        direction <- 0
      } else {
        direction <- 1
      }
    }
    if (asso_cox$direction[i] == "resistant") {
      if (resistant.cohort) {
        direction <- 1
      } else {
        direction <- 0
      }
    }
    if (asso_cox$treatment[i] %in% c("t2_me", "t1_me")) {
      sav_cox_strat <- sav_cox_fig5[(bem_cox_fig4_me[, biomarker] == direction) &
        (sav_cox_fig5[, strat_variable] == strat), ]
      bem_cox_strat <- bem_cox_fig4_me[(bem_cox_fig4_me[, biomarker] == direction) &
        (sav_cox_fig5[, strat_variable] == strat), ]
    }
    if (asso_cox$treatment[i] %in% c("t2_clin", "t1_clin")) {
      sav_cox_strat <- clin[clin[, biomarker] == direction, ]
      bem_cox_strat <- clin[clin[, biomarker] == direction, ]
    }
    sav_cox_strat$treatment_strat <- sav_cox_strat$treatment

    # Standardize treatment names
    tmp <- levels(sav_cox_strat$treatment_strat)
    sav_cox_strat$treatment_strat <- as.character(sav_cox_strat$treatment_strat)
    sav_cox_strat$treatment_strat[sav_cox_strat$treatment_strat == tmp[1]] <- "t1"
    sav_cox_strat$treatment_strat[sav_cox_strat$treatment_strat == tmp[2]] <- "t2"
    sav_cox_strat$treatment_strat <- factor(sav_cox_strat$treatment_strat)

    if (nrow(sav_cox_strat) > 0) {
      sav_cox_strat$treatment <- "yes"
    }
    m1 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")), # GENERALIZATION
      further = ""
    ), error = function(e) NULL)

    m2 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "PFS"
    ), error = function(e) NULL)

    m3 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "ORR"
    ), error = function(e) NA)

    m4 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "DFS"
    ), error = function(e) NULL)


    m5 <- tryCatch(survival_compare_treatment(
      df = sav_cox_strat,
      interaction = "treatment_strat",
      additional = paste0("+", paste(covariates, collapse = "+")),
      further = "", readout = "RFS"
    ), error = function(e) NULL)


    smodel <- NULL
    if (readout == "OS") {
      smodel <- list(m1)
    }
    if (readout == "PFS") {
      smodel <- list(m2)
    }
    if (readout == "ORR") {
      smodel <- list(m3)
    }
    if (readout == "DFS") {
      smodel <- list(m4)
      m1 <- m4
    }
    if (readout == "RFS") {
      smodel <- list(m5)
      m1 <- m5
    }


    res[[i]] <- tryCatch(list(
      # biomarker p
      p = (summary(m1))$coefficients[paste0("treatment_stratt2", ""), "Pr(>|z|)"],
      # biomarker eff
      eff = (summary(m1))$coefficients["treatment_stratt2", "coef"],
      confidence = if (readout == "ORR") {
        min(table(sav_cox_strat$treatment_strat[!is.na(sav_cox_strat$ORR)]))
      } else {
        min(table(sav_cox_strat$treatment_strat))
      },
      treatment_int_p = NA,
      treatmentp = (summary(m1))$coefficients["treatment_stratt2", "Pr(>|z|)"],
      better = NA,
      alteration = biomarker,
      matrix = list(as.data.frame(sav_cox_strat)),
      treatment_info = asso_cox$treatment[i],
      direction = asso_cox$direction[i],
      strat_variable = name,
      model = smodel,
      mutantorwildtype = direction,
      p_pfs = tryCatch((summary(m2))$coefficients[
        paste0("treatment_stratt2", ""),
        "Pr(>|z|)"
      ], error = function(e) NA), # biomarker p
      eff_pfs = tryCatch((summary(m2))$coefficients["treatment_stratt2", "coef"],
        error = function(e) NA
      ), # biomarker eff
      p_orrs = tryCatch(coef(summary(m3))[2, 4], error = function(e) NA),
      eff_orr = tryCatch(coef(summary(m3))[2, 1], error = function(e) NA)
    ), error = function(e) NULL)
  }

  if (length(res) == 1) {
    res <- do.call(bind_rows, list(res)) %>% mutate_at(c(
      "p", "eff",
      "confidence", "p_pfs", "eff_pfs", "p_orrs", "eff_orr"
    ), as.numeric)
  } else {
    res <- do.call(bind_rows, res) %>% mutate_at(c(
      "p", "eff",
      "confidence", "p_pfs", "eff_pfs", "p_orrs", "eff_orr"
    ), as.numeric)
  }

  res$treatment_info <- unlist(lapply(
    res$treatment_info,
    function(x) paste0(x, "_resi")
  ))

  # Put everthing together
  res <- rbind(res, res2)
  if (readout == "OS") {
    pvalue_here <- "P.Value"
    res$P.Value <- res$p
    res$effectsize <- res$eff
    res$p <- res$p
    res$eff <- res$eff
  }
  if (readout == "PFS") {
    pvalue_here <- "p_pfs"
    res$P.Value <- res$p_pfs
    res$effectsize <- res$eff_pfs
    res$p <- res$p_pfs
    res$eff <- res$eff_pfs
  }
  if (readout == "ORR") {
    pvalue_here <- "p_orrs"
    res$P.Value <- res$p_orrs
    res$effectsize <- res$eff_orr
    res$p <- res$p_orrs
    res$eff <- res$eff_orr
  }

  #
  if (readout == "DFS") {
    pvalue_here <- "P.Value"
    res$P.Value <- res$p
    res$effectsize <- res$eff
    res$p <- res$p
    res$eff <- res$eff
  }
  if (readout == "RFS") {
    pvalue_here <- "P.Value"
    res$P.Value <- res$p
    res$effectsize <- res$eff
    res$p <- res$p
    res$eff <- res$eff
  }


  res$fdr <- p.adjust(res[, pvalue_here, drop = TRUE], "BH")
  res$cfe <- res$alteration %>% make.unique()
  res$strat <- "1"
  res$gene <- res$alteration
  res$strat <- res$treatment_info

  # Change the variables for visualization,
  # resistant-markers and resistance-cohorts
  tmp <- unlist(lapply(
    seq_len(nrow(res)),
    function(i) {
      gsub(
        "sens$", res$direction[i],
        gsub("resi$", res$direction[i], res$strat[i])
      )
    }
  ))
  tmp <- unlist(lapply(
    seq_len(length(tmp)),
    function(i) {
      gsub(
        "sensitive", "sens",
        gsub("resistant", "resi", res$strat[i])
      )
    }
  ))
  res$direction <- unlist(lapply(
    res$strat,
    function(x) strsplit(x, "_")[[1]][3]
  )) %>%
    dplyr::recode("sens" = "sensitive", "resi" = "resistant")
  res$strat <- tmp

  # filter unneeded columns
  res <- res %>% dplyr::select(-c("p_pfs", "eff_pfs", "p_orrs", "eff_orr"))
  return(res)
}
#################################################



#################################################
#' Compare Survival Treatment
#'
#' helper function for comparison surival data
#'
#' @param df dataframe of intrest
#' @param interaction glm interaction test
#' @param additional additional covariates
#' @param further more covariates
#' @param interaction.with default: treatment_strat
#' @param readout default "OS"
#'
#' @noRd

survival_compare_treatment <- function(df,
                                       interaction,
                                       additional,
                                       further,
                                       interaction.with = "treatment_strat",
                                       readout = "OS") {
  if (readout %in% c("ORR")) {
    model <- glm(as.formula(paste0(
      "df$", readout, " ~
              ", interaction, "", additional, further, " + treatment_strat  +",
      interaction.with, "*", interaction, ""
    )),
    data = df,
    family = binomial(link = "logit")
    )
  }

  if (readout %in% c("OS", "PFS", "DFS", "RFS")) {
    model <- survival::coxph(as.formula(paste0(
      "survival::Surv(time = df$",
      readout, ".months, event = df$", readout, ".event) ~
    ", interaction, "", additional, further, " + treatment_strat  +",
      interaction.with, "*", interaction, ""
    )), data = df)
  }
  return(model)
}
#################################################


#' My Stringsplit
#'
#' helper function for string splitting
#'
#' @param x string
#' @param splits what splits
#' @param ...
#'
#' @return string
#'
#' @noRd
my_strsplit <- function(x, splits, ...) {
  for (split in splits)
  {
    x <- unlist(strsplit(x, split, ...))
  }
  return(x[!x == ""]) # Remove empty values
}



#################################################
#' Get Anno Combinations
#'
#' helper function
#'
#' @param asso
#'
#' @noRd
get_anno_combinations <- function(asso) {
  anno <- as.data.frame(stringi::stri_list2matrix(lapply(
    asso$gene,
    function(x) {
      gsub(
        "\\)", "",
        gsub("\\(", "", my_strsplit(x, c("_and_", "_or_")))
      )
    }
  ), byrow = TRUE))
  columns <- paste0("V", seq_len(ncol(anno)))
  anno <- anno %>%
    mutate(.id = seq_len(length(asso$gene))) %>%
    tidyr::pivot_longer(cols = columns) %>%
    mutate(values = 1) %>%
    dplyr::select(-name) %>%
    tidyr::pivot_wider(
      values_from = values, names_from = value,
      values_fn = list(values = function(x) unique(x))
    ) %>%
    replace(is.na(.), 0) %>%
    tibble::column_to_rownames(".id")
  row.names(anno) <- make.unique(asso$gene)
  return(anno)
}
#################################################


#################################################
#' Add Medians
#'
#' helper function for adding medians to summary_table
#'
#' @param summary make_table dataframe
#'
#' @return make_table
#'
#' @noRd
add_medians <- function(summary) {
  medians <- lapply(seq_len(nrow(summary)), function(i) {
    data <- summary[i, "matrix"][[1]][[1]]
    data <- data[names(summary[i, "model"][[1]][[1]]$y), ]
    if (summary$readout[i] %in% c("OS", "PFS", "DFS", "RFS")) {
      tmp1 <- tryCatch((survminer::surv_median(survminer::surv_fit(
        formula = summary[i, "model"][[1]][[1]]$y ~ treatment_strat,
        data = data
      )) %>%
        tibble::column_to_rownames("strata"))[, "median", drop = FALSE] %>% t(),
      error = function(e) c(NA, NA)
      )
      tmp <- tryCatch(as.data.frame(table(data$treatment_strat)),
        error = function(e) c(NA, NA)
      )
      tmp2 <- tmp %>%
        group_by(Var1) %>%
        summarise(sum = sum(Freq)) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("Var1") %>%
        t()
      if (all(is.na(tmp1))) {
        tmp1 <- data.frame(a = NA, b = NA)
        colnames(tmp1) <- paste0("treatment_strat=", colnames(tmp2))
      }
      tmp <- cbind(tmp1, tmp2)
    }
    if (summary$readout[i] %in% c("ORR")) {
      tmp <- tryCatch(as.data.frame(table(data$treatment_strat, data[, "ORR"])),
        error = function(e) c(NA, NA)
      )
      if (all(is.na(tmp))) {
        tmp <- data.frame(a = NA, b = NA)
        colnames(tmp) <- paste0("treatment_strat=", colnames(tmp))
      }
      if (any(!is.na(tmp))) {
        tmp1 <- tmp %>%
          group_by(Var1) %>%
          summarise(diff = Freq[2] / (sum(Freq))) %>%
          tibble::remove_rownames() %>%
          tibble::column_to_rownames("Var1") %>%
          t()
        tmp2 <- tmp %>%
          group_by(Var1) %>%
          summarise(sum = sum(Freq)) %>%
          tibble::remove_rownames() %>%
          tibble::column_to_rownames("Var1") %>%
          t()
        colnames(tmp1) <- paste0("treatment_strat=", colnames(tmp2))
        tmp <- cbind(tmp1, tmp2)
      }
    }
    return(tmp)
  })
  medians <- as.data.frame(do.call(rbind, medians))
  colnames(medians) <- c("mean1", "mean2", "n1", "n2")

  medians <- medians %>% mutate_if(is.numeric, round, 2)
  summary <- cbind(summary, medians)
  return(summary)
}
#################################################



#################################################
## Change medians
#'
#' Helper functions changing median for binary readouts;
#' change median survival for step 3 summary table;
#' please harmonize wtih change_median
#'
#' @param summary
#'
#' @noRd
change_medians <- function(summary) {
  medians <- lapply(seq_len(nrow(summary)), function(i) {
    data <- summary[i, "medians"]
    if (summary$readout[i] %in% c("OS", "PFS", "DFS", "RFS")) {
      tmp <- data[[1]][, c("strata", "median")] %>%
        tibble::remove_rownames %>%
        tibble::column_to_rownames("strata") %>%
        base::t()
      colnames(tmp) <- unlist(lapply(
        colnames(tmp),
        function(x) strsplit(x, "=")[[1]][2]
      ))
    }
    if (summary$readout[i] %in% c("ORR")) {
      tmp <- data[[1]][[1]] %>%
        group_by(marker) %>%
        summarise(diff = Freq[2] / (sum(Freq))) %>%
        tibble::remove_rownames %>%
        tibble::column_to_rownames("marker") %>%
        base::t()
    }
    return(tmp)
  })
  medians <- as.data.frame(do.call(rbind, medians))
  colnames(medians) <- c("mean1", "mean2")
  medians <- medians %>% mutate_if(is.numeric, round, 2)
  summary <- cbind(summary, medians)
  return(summary)
}
#################################################

#' Redistributed
#'
#' first helper function that checks if mutually exclusive modules in subtypes contain enough reclassified mutations compared to baseline
#'
#' @param alteration column name of mutations_df
#' @param mutations data_mutations dataframe
#' @param clinical data_clinical dataframe
#' @param subtype_column strat in calc_asso_cox
#' @param subtype which subpopulation
#' @param treatment what treatment
#'
#' @noRd
redistributed <- function(alteration, # column name of mutations_df
                          mutations, # data_mutations dataframe
                          clinical, # data_clinical dataframe
                          subtype_column, # strat in calc_asso_cox
                          subtype, # which subpopulation
                          treatment # what treatment
) {
  names <- na.omit(clinical$sample[(clinical[, as.character(subtype_column)] == subtype) &
    (clinical$treatment == as.character(treatment))])
  local_muts <- mutations[mutations$sample %in% names, ]
  local_mut <- local_muts[, alteration]

  containing_alterations <- unlist(strsplit(alteration, "_or_"))
  local_mut_derivatives <- local_muts[, containing_alterations, drop = FALSE]

  redistributed <- lapply(seq_len(ncol(local_mut_derivatives)), function(x) {
    tmp <- table(local_mut_derivatives[, x], local_mut)
    diag(tmp) <- NA
    tmp <- sum(na.omit(as.numeric(tmp)))
    # returns the redistributed elements of ME modules compared to its baseline
    return(tmp)
  })
  names(redistributed) <- colnames(local_mut_derivatives)
  if (length(containing_alterations) == 1) {
    redistributed[redistributed == 0] <- Inf
  }
  # how many were distributed by MEs compared to the most similar one gene baseline
  return(redistributed[which.min(redistributed)])
}



#' Get statistics
#'
#' Get confidence interval and effect sizes for association tests
#'
#' @param df
#'
#' @noRd
get_forest_new <- function(df) {
  testFunction <- function(x, y) {
    family <- tryCatch(x$family$family, error = function(e) NULL)
    if (is.null(family)) {
      family <- "none"
    }
    if (y == 1) {
      tmp <- tryCatch(x$conf.int["biomarker1", "exp(coef)"],
        error = function(e) NULL
      )
      tmp <- tryCatch(x$conf.int[1, "exp(coef)"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(-x$coefficients[2, "Estimate"], error = function(e) NULL)
      }
    }
    if (y == 2) {
      tmp <- tryCatch(x$conf.int["biomarker1", "upper .95"],
        error = function(e) NULL
      )
      tmp <- tryCatch(x$conf.int[1, "upper .95"], error = function(e) NULL)
      if (family == "binomial") {
        tmp <- tryCatch(-x$coefficients[2, "Estimate"] +
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual),
        error = function(e) NULL
        )
      }
    }
    if (y == 3) {
      tmp <- tryCatch(x$conf.int["biomarker1", "lower .95"],
        error = function(e) NULL
      )
      tmp <- tryCatch(x$conf.int[1, "lower .95"], error = function(e) NULL)
      # x$conf.int[1,"upper .95"]
      if (family == "binomial") {
        tmp <- tryCatch(-x$coefficients[2, "Estimate"] -
          x$coefficients[2, "Std. Error"] * qt(p = 0.975, df = x$df.residual),
        error = function(e) NULL
        )
      }
    }
    if (is.null(tmp)) {
      tmp <- NA
    }
    return(tmp)
  }
  df$hazard <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 1)
  ) %>% unlist()
  df$lower <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 2)
  ) %>% unlist()
  df$upper <- lapply(
    seq_len(nrow(df)),
    function(x) testFunction((summary(df$model[x][[1]][[1]])), y = 3)
  ) %>% unlist()
  return(df)
}

# Oncoprint config
alter_fun <- list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = grid::gpar(fill = "white", col = NA)
    )
  },
  LOSS = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = grid::gpar(fill = "darkblue", col = NA)
    )
  },
  AMP = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = grid::gpar(fill = "darkred", col = NA)
    )
  },
  MUT = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33,
      gp = grid::gpar(fill = "grey30", col = NA)
    )
  },
  MSI = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33,
      gp = grid::gpar(fill = "orange", col = NA)
    )
  },
  `NA` = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33,
      gp = grid::gpar(fill = "white", col = NA)
    )
  },
  FEAT = function(x, y, w, h) {
    grid::grid.rect(x, y, w - unit(0.5, "mm"), h * 0.2,
      gp = grid::gpar(fill = "darkorange", col = NA)
    )
  }
)


#' Make matrices
#'
#' Prepare matrices for multiple regression
#'
#' @param clinical clinical dataframe
#' @param mutations mutational dataframe
#' @param treatment_variable treatment column
#' @param treatment treatment factor string
#' @param include_2_order_interactions logical if including 2nd order
#' interactions
#' @param include_3_order_interactions logical if including 3rd order
#' interactions
#' @param subtype_column subtype column in clinical data
#' @param survival_var OS, PFS or ORR, for survival please name for example
#' OS.months and OS.event
#' @param excluded_clinical vector of excluded clinical covariates
#' @param min_abundance minimum number of tumours in mutations for
#' adding to the model
#'
#' @return list containing mutational and clinical dataframes
#'
#' @noRd
prepare_matrices <- function( ### Takes mutations and clinical data, imputes, builds and saves model
                             ### matrices
                             ######################################################
                             clinical,
                             mutations,
                             treatment_variable,
                             treatment,
                             include_2_order_interactions = FALSE,
                             include_3_order_interactions = FALSE,
                             subtype_column = NULL,
                             survival_var,
                             excluded_clinical = NULL,
                             min_abundance = 10) {
  treatment <- as.character(treatment)
  if (treatment != "yes") {
    clinical <- clinical[clinical[, treatment_variable] == treatment, ]
  }
  if (survival_var == "OS") {
    survival <- data.frame(
      vitalStatus = clinical$OS.event,
      overallSurvival = clinical$OS.months
    )
    row.names(survival) <- clinical[, "sample"]
  }
  if (survival_var == "PFS") {
    survival <- data.frame(
      vitalStatus = clinical$PFS.event,
      overallSurvival = clinical$PFS.months
    )
    row.names(survival) <- clinical[, "sample"]
  }
  if (survival_var == "ORR") {
    survival <- data.frame(overallSurvival = clinical$ORR, drop = NA)
    row.names(survival) <- clinical[, "sample"]
    survival <- survival[!is.na(survival$overallSurvival), , drop = FALSE]
  }
  if (survival_var == "none") {
    survival <- clinical
    row.names(survival) <- clinical[, "sample"]
  }
  if (survival_var == "DFS") {
    survival <- data.frame(
      vitalStatus = clinical$DFS.event,
      overallSurvival = clinical$DFS.months
    )
    row.names(survival) <- clinical[, "sample"]
  }
  if (survival_var == "RFS") {
    survival <- data.frame(
      vitalStatus = clinical$RFS.event,
      overallSurvival = clinical$RFS.months
    )
    row.names(survival) <- clinical[, "sample"]
  }
  rowname <- row.names(survival)
  if (survival_var != "none") {
    mutations <- ((mutations %>% tibble::remove_rownames() %>%
      tibble::column_to_rownames("sample"))[row.names(survival), ]) %>%
      mutate_all(as.numeric)
  } else {
    mutations <- ((mutations %>% tibble::remove_rownames() %>%
      tibble::column_to_rownames("sample"))) %>%
      mutate_all(as.numeric)
  }

  mutations <- mutations %>% mutate_all(impute_median)

  mutations <- mutations[, apply(
    mutations, 
    2,
    function(x) sum(x)) >= min_abundance]
  row.names(mutations) <- rowname

  if (treatment == "yes") {
    clinical <- clinical[, !names(clinical) %in% excluded_clinical] %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("sample")
  } else {
    clinical <- clinical[clinical[, treatment_variable] == treatment, !names(clinical) %in% excluded_clinical] %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("sample")
  }

  # make model matrix
  covariates <- merge(mutations, clinical, by = 0) %>%
    tibble::column_to_rownames("Row.names")

  # delete constant rows
  covariates <- covariates[, !(sapply(covariates,function(x) 
    (is.character(x) | is.factor(x)) & length(unique(x)) < 2))]

  if (treatment != "yes") {
    if (include_2_order_interactions) {
      covariates <- covariates[, !colnames(covariates) %in% treatment_variable]
      message(as.formula(paste0("~ .*(",paste(subtype_column, collapse = "+"),")")))

      covariates <- model.matrix(as.formula(paste0("~ .*(",
        paste(subtype_column, collapse = "+"), ")")),
      data = covariates,
      contrasts.arg = lapply(covariates[, sapply(covariates, is.factor),
                              drop = FALSE], 
                              contrasts, 
                              contrasts = FALSE))
    } else {
      covariates <- covariates[, !colnames(covariates) %in% treatment_variable]
      covariates <- model.matrix(~.,
        data = covariates,
        contrasts.arg = lapply(covariates[, sapply(covariates, is.factor),
          drop = FALSE], 
          contrasts, 
          contrasts = FALSE)
      )
    }
  } else {
    if (include_3_order_interactions) {
      message(as.formula(paste0("~ .*treatment*(",paste(subtype_column, collapse = "+"),")")))
      covariates <- model.matrix(as.formula(paste0(
        "~ .*treatment*(",
        paste(subtype_column, collapse = "+"), ")"
      )),
      data = covariates,
      contrasts.arg = lapply(covariates[, sapply(covariates, is.factor),
            drop = FALSE], 
            contrasts, 
            contrasts = FALSE)
      )
    } else {
      covariates <- model.matrix(~ . * treatment,
        data = covariates,
        contrasts.arg = lapply(covariates[, sapply(covariates, is.factor),
          drop = FALSE], 
          contrasts, 
          contrasts = FALSE)
      )
    }
  }

  if (survival_var != "none") {
    int <- intersect(row.names(survival), row.names(covariates))
  } else {
    int <- row.names(covariates)
  }

  if (survival_var %in% c("OS", "PFS", "DFS", "RFS")) {
    survival <- survival[int, c("vitalStatus", "overallSurvival"), drop = FALSE]
  }
  if (survival_var %in% c("ORR")) {
    survival <- cbind(
      survival[int, c("overallSurvival"), drop = FALSE],
      survival[int, c("overallSurvival"), drop = FALSE]
    )
  }
  covariates <- covariates[int, apply(covariates, 2,function(x) sum(x)) > 
                             min_abundance]

  survival <- if (is.null(survival)) {
    survival <- survival
  } else {
    survival <- as.matrix(survival)
  }
  return(list(
    survival = survival,
    covariates = as.matrix(covariates),
    survival_var = survival_var
  ))
} ######################################################



#' Calculate Mutual Exclusivity
#'
#' run the mutex algorithm, deprecated
#'
#' @param data_mut A binary event matrix, here mutation: yes, no
#' @param mutex_output_exists If object from java mutex already exists,
#' you can specify TRUE and use the file
#' in metadata/ranked-groups.txt instead.
#' @param min_variants mimimum number of variants for searching
#' mutually exclusive modules
#' @param prior_modules vector of one prior modules which is known and in use
#' @param save saving location for found modules
#'
#' @return List of 2: First element: Mutually exclusive binary event matrix
#'          Second element: table of mutually exclusive genes ranked by their
#'          mutual exclusivity score
#'
#' @noRd

calculate_mututal_exclusivity <- function(data_mut = data_mutations,
                                          mutex_output_exists,
                                          min_variants = 10,
                                          prior_modules = NULL,
                                          save = "/home/metadata") {
  bem_mutex <- data_mut %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample")
  ### get only variants that more abundant than 10
  bem_mutex <- bem_mutex[, as.logical(apply(
    bem_mutex, 2,
    function(x) sum(as.numeric(x))
  ) >= min_variants)]
  ngs <- make_mutex_df(bem_mutex)
  ngs <- ngs %>% tibble::rownames_to_column("Symbol")
  ngs[is.na(ngs)] <- 0
  write.table(
    x = ngs,
    file = "metadata/mutex_ngs_alterations_min10.csv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  if (mutex_output_exists == FALSE) {
    message("Calculating mutual exclusivity")
    system(paste("java -jar /mutex/target/mutex.jar", save), intern = TRUE)
    ngs_exclusive <- read_ranked_groups(save = save)
  }

  if (mutex_output_exists == TRUE) {
    message("Took existing mutual exclusivity table")
    ngs_exclusive <- read_ranked_groups(save = save)
  }
  colnames(ngs_exclusive) <- ngs_exclusive[1, ] %>%
    unlist()
  ngs_exclusive <- ngs_exclusive[2:nrow(ngs_exclusive), ]
  colnames(ngs_exclusive) <- as.character(seq_len(ncol(ngs_exclusive)))
  ngs_exclusive <- ngs_exclusive %>% mutate_all(as.character)


  ### add previous knowledge
  if (length(prior_modules) > 0) {
    know <- prior_modules
    ngs_exclusive[, (ncol(ngs_exclusive) + 1):(ncol(ngs_exclusive) + length(know))] <- ""
    j <- nrow(ngs_exclusive) + 1
    for (i in seq_len(nrow(ngs_exclusive))) {
      alts <- unlist(lapply(
        ((ngs_exclusive[i, 2:(ncol(ngs_exclusive))])),
        as.character
      ))
      if (any(know %in% alts)) {
        newalts <- c(alts[alts != ""], know[!know %in% alts])
        newalts <- c(newalts, rep("", ncol(ngs_exclusive) - 1 - length(newalts)))
        ngs_exclusive[j, ] <- c(1, newalts)
        j <- j + 1
      }
    }

    newalts <- c(alts[alts != ""], know[!know %in% alts])
    newalts <- c(newalts, rep("", ncol(ngs_exclusive) - 1 - length(newalts)))
    ngs_exclusive[j, ] <- c(0, know, rep("", ncol(ngs_exclusive) - 1 - length(know)))
  }
  ngs_exclusive$name <- unlist(lapply(
    seq_len(nrow(ngs_exclusive)),
    function(x) {
      tmp <- as.character(unlist(ngs_exclusive[x, 2:ncol(ngs_exclusive)]))
      tmp <- tmp[tmp != "" && tmp != "NA"]
      return(paste(tmp, collapse = ";"))
    }
  ))

  # remove duplicate modules
  ngs_exclusive <- ngs_exclusive[
    !duplicated(lapply(
      seq_len(nrow(ngs_exclusive)),
      function(x) {
        tmp <- as.character(unlist(ngs_exclusive[x, 2:(ncol(ngs_exclusive) - 1)]))
        tmp <- tmp[tmp != ""]
        return(sort(unlist(tmp)))
      }
    )),
  ]

  # adds in _SV, _AMP, and _DEL
  data_mutations_me <- data_mut
  for (j in seq_len(nrow(ngs_exclusive))) {
    data_mutations_me <- add_me(
      data = data_mutations_me,
      names = as.character(unlist(ngs_exclusive[j, 2:(ncol(ngs_exclusive) - 1)])),
      logic = "or"
    )
  }
  data_mutations_me <- data_mutations_me[, !duplicated(colnames(data_mutations_me))]
  save(data_mutations_me, file = "metadata/data_mutations_me_original.RData")

  # replace all ; in the end with nothing, cleaner for plotting
  ngs_exclusive$name <- gsub(";*$", "", ngs_exclusive$name)

  message(
    "Number of single genes: ",
    as.character(ncol(data_mut))
  )

  message(
    "Number of single plus modules: ",
    as.character(ncol(data_mutations_me))
  )

  return(list(
    "data_mutations_me" = data_mutations_me,
    "ngs_exclusive" = ngs_exclusive
  ))
}



#' Filter
#'
#' Filtering exclusivity dataframe
#'
#' @param exclusivities dataframe of somatic mutations found with the
#' mutual exclusivity calculation
#' @param data_mut A binary event matrix, here mutation: yes, no
#' @param found_modules character vector of predictive
#' modules of somatic mutations
#' @param prior_module prior knowledge on existing multi-gene biomarker
#'
#' @return updated dataframe
#'
#' @noRd
filter_exclusivities <- function(data_mutations,
                                 ngs_exclusive,
                                 exclusivities,
                                 data_mut,
                                 found_modules,
                                 prior_module) {
  ngs_alts <- c()
  for (x in seq_len(nrow(exclusivities))) {
    # adds in _SV, _AMP, and _DEL
    ngs_alts[x] <- add_me(
      data = data_mutations,
      names = as.character(unlist(ngs_exclusive[x, 2:(ncol(ngs_exclusive) - 1)])),
      logic = "or"
    )[, ncol(data_mutations) + 1, drop = FALSE] %>% colnames()
  }

  exclusivities$identified <- (ngs_alts %in% found_modules)
  exclusivities$without_prior <- factor(unlist(lapply(
    exclusivities$identifier,
    function(x) paste0("", paste(x[!(x %in% prior_module)], collapse = "/"))
  )))
  exclusivities$which.with.prior <- factor(unlist(lapply(
    exclusivities$identifier,
    function(x) all(prior_module %in% x)
  )))
  return(exclusivities)
}



#' redistributed2
#'
#' second helper function that checks if mutually exclusive modules in subtypes
#' contain enough reclassified mutations for each added
#' (currently unused)
#'
#' @param alteration column name of mutations_df
#' @param mutations data_mutations dataframe
#' @param clinical data_clinical dataframe
#' @param subtype_column strat in calc_asso_cox
#' @param subtype which subpopulation
#' @param treatment what treatment
#'
#' @noRd
redistributed2 <- function(alteration, # column name of mutations_df
                           mutations, # data_mutations dataframe
                           clinical, # data_clinical dataframe
                           subtype_column, # strat in calc_asso_cox
                           subtype, # which subpopulation
                           treatment # what treatment
) {
  names <- na.omit(clinical$sample[(clinical[, as.character(subtype_column)] == subtype) &
    (clinical$treatment == as.character(treatment))])
  local_muts <- mutations[mutations$sample %in% names, ]
  local_mut <- local_muts[, alteration]

  containing_alterations <- unlist(strsplit(alteration, "_or_"))
  local_mut_derivatives <- local_muts[, containing_alterations, drop = FALSE]

  # returns the redistributed elements of ME modules compared to its baseline
  redistributed <- lapply(seq_len(ncol(local_mut_derivatives)), function(x) {
    tmp <- length(which(local_mut_derivatives[, x] == "1"))
    return(tmp)
  })
  names(redistributed) <- colnames(local_mut_derivatives)
  if (length(containing_alterations) == 1) {
    redistributed[redistributed == 0] <- Inf
  }

  # how many were distributed by MEs compared to the most
  # similar one gene baseline
  names(redistributed) <- unlist(lapply(
    names(redistributed),
    function(x) strsplit(x, "_")[[1]][1]
  ))
  redistributed <- sort(unlist(redistributed), decreasing = TRUE)
  return(redistributed)
}



##############################################################
#' Helper for calculate mutual exclusive
#'
#' @param save path for ranked-groups.txt
#'
#' @return imported output of the mutex algorithm
#'
#' @noRd
read_ranked_groups <- function(save) {
  max_elements <- 0
  message("Looking for ", save, "/ranked-groups.txt", " ...")
  count_members <- read.delim(paste0(save, "/ranked-groups.txt"),
    header = FALSE, sep = "\n"
  )
  for (i in seq_len(nrow(count_members))) {
    if (length(unlist(stringr::str_split(count_members[i, ],
      pattern = "\t"
    ))) > max_elements) {
      max_elements <- length(unlist(stringr::str_split(count_members[i, ],
        pattern = "\t"
      )))
    }
  }
  ngs_exclusive <- read.delim(paste0(save, "/ranked-groups.txt"),
    header = FALSE,
    sep = "\t",
    col.names = seq_len(max_elements)
  )

  return(ngs_exclusive)
}


##############################################################
#' Helper for calculating n-th lowest number in vector
#'
#' @param x numeric vector
#' @param n integer for n-th lowest in x
#'
#' @return numeric
#'
#' @noRd
which_nth_lowest <- function(x, n) {
  if (length(na.omit(unique(x))) < n) { # get the extrema closest
    n <- length(na.omit(unique(x)))
    if (n == 0) {
      n <- 1
    }
  }

  for (i in seq_len(n - 1L)) x[x == min(x, na.rm = TRUE)] <- Inf
  res <- which(x == min(x, na.rm = TRUE))

  if (length(x) < n) {
    res <- 0
  }
  return(res)
}


##############################################################
#' Helper for calculating n-th lowest number in vector
#'
#' @param x numeric vector
#' @param n integer for n-th highest in x
#'
#' @return numeric
#'
#' @noRd
which_nth_highest <- function(x, n) {
  if (length(na.omit(unique(x))) < n) { # get the extrema closest
    n <- length(na.omit(unique(x)))
    if (n == 0) {
      n <- 1
    }
  }

  for (i in seq_len(n - 1L)) x[x == max(x, na.rm = TRUE)] <- -Inf
  res <- which(x == max(x, na.rm = TRUE))

  if (length(x) < n) {
    res <- 0
  }
  return(res)
}
