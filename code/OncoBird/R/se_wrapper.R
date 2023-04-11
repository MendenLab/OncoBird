########################################################################
# WRAPPING FUNCTIONS for SUMMARIZEDEXPERIMENT



#' Calculate subtype enrichment
#'
#' Calculating hypergeometric tests for two subtype columns to check for
#' associations.
#'
#' @param se SE object
#' @param col_label name of subtype 1
#' @param row_label name of subtype 2
#' @param digits how many digits for the reported p-value
#'
#' @return SE object with metadata containing subtype enrichment 
#' tests 'subtype_enrichment'
#' @export
#'
#' @examples {
#' 
#' simed_data <- sim_data()
#'   
#' se <- prepare_data(
#'     data = simed_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(simed_data)[1:5])
#'     
#' se <- cl_subtype_enrichment(
#'     se = se,
#'     col_label = "s",
#'     row_label = "r",
#'     digits = 3)
#'   
#' }
#'
cl_subtype_enrichment <- function(se,
                                  col_label,
                                  row_label,
                                  digits) {

  # Check input
  assertthat::assert_that(is.character(col_label))
  assertthat::assert_that(is.character(row_label))
  assertthat::assert_that(length(rowData(se)[[row_label]]) > 0)
  assertthat::assert_that(length(rowData(se)[[col_label]]) > 0)
  assertthat::assert_that(
    length(rowData(se)[[row_label]]) == length(rowData(se)[[col_label]])
  )
  assertthat::assert_that(length(assay(se)) > 0)

  # Calculate Enrichment
  subtype_enrichment <- enrichment(
    tab = table(
      rowData(se)[[row_label]],
      rowData(se)[[col_label]]
    ),
    col_label = col_label,
    row_label = row_label,
    digits = digits
  )

  # Save to metadata(se)
  metadata(se)$subtype_enrichment <- subtype_enrichment
  metadata(se)$wf_meta <- list(
    col_subtype_enrichment = col_label,
    row_subtype_enrichment = row_label
  )

  return(se)
}



########################################################################
#' Plot subtype enrichment
#'
#' Plotting results after 'cl_subtype_enrichment'
#'
#' @param se SE object
#' @param limits plot color limits
#'
#' @return ggplot object
#' @export
#'
#' @examples {
#'   
#'   simed_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = simed_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(simed_data)[1:5])
#'     
#'   se <- cl_subtype_enrichment(
#'     se = se,
#'     col_label = "s",
#'     row_label = "r",
#'     digits = 3
#'   )
#'   
#'   pl_subtype_enrichment(se)
#'   
#' }
#'
pl_subtype_enrichment <- function(se, limits = c(-4,4)) {
  # Check input
  assertthat::assert_that(is.data.frame(metadata(se)$subtype_enrichment))
  col_label <- metadata(se)$wf_meta[["col_subtype_enrichment"]]
  row_label <- metadata(se)$wf_meta[["row_subtype_enrichment"]]

  assertthat::assert_that(
    col_label %in% colnames(metadata(se)$subtype_enrichment)
  )
  assertthat::assert_that(
    row_label %in% colnames(metadata(se)$subtype_enrichment)
  )

  # Plot Enrichment
  plot_subtype_enrichment(
    subtype_enrichment_df = metadata(se)$subtype_enrichment,
    col_label = col_label,
    row_label = row_label,
    limits = limits
  )
}



########################################################################
#' Calculate enrichment for mutations in subtypes
#'
#' Running hypergeometric tests for mutations in specified subtypes
#'
#' @param se SE object
#' @param sample_column subtype column
#' @param min_mutants minimum number of mutants for conducting a
#' statistical test
#' @param clin_alterations binary clinical columns to join with mutations
#' @param subtype character vector containing subtypes to screen for enrichments
#'
#' @return SE object with metadata containing subtype enrichment tests
#' 'enrichment_genomics'
#' @export
#'
#' @examples {
#'   
#'  sim_data <- sim_data()
#'   
#'  se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'  se <- cl_subtype_enrichment(
#'     se = se,
#'     col_label = "s",
#'     row_label = "r",
#'     digits = 3
#'  )
#'   
#'  se <- cl_enrichment_genomics(se,
#'     sample_column = "sample",
#'     min_mutants = 10,
#'     subtype = c("s", "r")
#'   )
#' }
#'
cl_enrichment_genomics <- function(se,
                                   sample_column,
                                   min_mutants,
                                   clin_alterations = NULL,
                                   subtype) {
  assertthat::assert_that(length(subtype) < 3 && length(subtype) > 0)
  if (length(subtype) == 1) {
    # name new object according to site of interest
    enrichment_name <- paste0(subtype, "_mutations_enrichment")

    site_mutation_enrichment <- enrichment_genomics(
      bem = assay(se),
      clin = rowData(se),
      sample_column = sample_column,
      min_mutants = min_mutants,
      clin_alterations = clin_alterations,
      subtype = subtype
    )

    metadata(se)[[enrichment_name]] <- site_mutation_enrichment

    workflow_metadata <- metadata(se)$wf_meta
    workflow_metadata[["enrichment_genomics"]] <- enrichment_name
    metadata(se)$wf_meta <- workflow_metadata

    return(se)
  } else {
    i <- 0
    for (type in subtype) {
      i <- i + 1
      enrichment_name <- paste0(type, "_mutations_enrichment")

      site_mutation_enrichment <- enrichment_genomics(
        bem = assay(se),
        clin = rowData(se),
        sample_column = sample_column,
        min_mutants = min_mutants,
        clin_alterations = clin_alterations,
        subtype = type
      )

      metadata(se)[[enrichment_name]] <- site_mutation_enrichment

      workflow_metadata <- metadata(se)$wf_meta
      multiple_name <- paste0(i, "_enrichment_genomics")
      workflow_metadata[[multiple_name]] <- enrichment_name
      metadata(se)$wf_meta <- workflow_metadata
    }
    return(se)
  }
}



########################################################################
#' Plot enrichment for mutations in subtypes
#'
#' Plotting results after 'cl_enrichment_genomics'
#'
#' @param se SE object
#' @param p_value p-value threshold for plotting results
#' @param fdr false discovery rate threshold for reporting results
#'
#' @return ggplot object
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'   se <- cl_subtype_enrichment(
#'     se = se,
#'     col_label = "s",
#'     row_label = "r",
#'     digits = 3
#'   )
#'   
#'   se <- cl_enrichment_genomics(se,
#'     sample_column = "sample",
#'     min_mutants = 10,
#'     subtype = c("s", "r")
#'   )
#'   
#'   # no siginficants expected
#'   pl_enrichment_genomics(se, p_value = 0.5, fdr = 1)
#'   
#' }
#'
pl_enrichment_genomics <- function(se,
                                   p_value = 0.05,
                                   fdr = 0.2) {
  colors <- c(
    "orange", "darkblue", "pink",
    "darkgreen", "mediumturquoise", "salmon"
  )

  if ("1_enrichment_genomics" %in% names(metadata(se)$wf_meta) &&
    "2_enrichment_genomics" %in% names(metadata(se)$wf_meta)) {
    enrichment_name_1 <- metadata(se)$wf_meta[["1_enrichment_genomics"]]
    site_mutation_enrichment_1 <- metadata(se)[[enrichment_name_1]]

    enrichment_name_2 <- metadata(se)$wf_meta[["2_enrichment_genomics"]]
    site_mutation_enrichment_2 <- metadata(se)[[enrichment_name_2]]

    label_1 <- as.character(unique(site_mutation_enrichment_1$alteration)) 
    label_2 <- as.character(unique(site_mutation_enrichment_2$alteration))


    plot_subtype_mutations_enrichment(
      enr_list = list(
        site_mutation_enrichment_1, site_mutation_enrichment_2
      ),
      p_value_threshold = p_value,
      fdr_value_threshold = fdr,
      clin = as.data.frame(rowData(se)),
      bem = as.data.frame(assay(se)),
      labels_list = list(label_1, label_2),
      colors_list = list(
        colors[seq_len(length(label_1))],
        colors[(1 + length(label_1)):(length(label_1) + length(label_2))]
      )
    )
  } else {
    enrichment_name <- metadata(se)$wf_meta[["enrichment_genomics"]]
    site_mutation_enrichment <- metadata(se)[[enrichment_name]]

    label <- as.character(unique(site_mutation_enrichment$alteration))

    plot_subtype_mutations_enrichment(
      enr_list = list(site_mutation_enrichment),
      p_value_threshold = p_value,
      fdr_value_threshold = fdr,
      clin = as.data.frame(rowData(se)),
      bem = assay(se),
      labels_list = list(label),
      colors_list = list(colors[seq_len(length(label))])
    )
  }
}



########################################################################
#' Calculate mutual exclusivity
#'
#' Run the mutex algorithm on the mutational data
#'
#' @param se SE object
#' @param mutex_output_exists logical specifying if mutex output already exists
#' @param mutex_path path for mutex software
#' @param min_variants minimum number of patients mutated for analysing for
#' mutual exclusivity
#' @param prior_modules character vector of prior module to add to the modules
#' @param save path for existing mutex output
#'
#' @return SE object with metadata containing subtype enrichment tests
#' 'me_modules' and 'ngs_exclusive'
#' @export
#'
#' @examples \dontrun{
#'   sim_data <- sim_data()
#'
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = c("OS", "DFS"),
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "PatientID",
#'     treatment_column = "Adj",
#'     mutation_columns = colnames(sim_data)[14:ncol(sim_data)]
#'   )
#'
#'
#'   se <- cl_mututal_exclusivity(se,
#'     min_variants = 10,
#'     mutex_output_exists = TRUE,
#'     save = "metadata/",
#'     mutex_path = "path/to/mutex/installation"
#'   )
#' }
#'
cl_mututal_exclusivity <- function(se,
                                   mutex_output_exists,
                                   mutex_path = "/",
                                   min_variants = 10,
                                   prior_modules = NULL,
                                   save = "/home/metadata") {
  data_mut <- as.data.frame(SummarizedExperiment::assay(se))
  bem_mutex <- data_mut %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample")
  ### get only variants that more abundant than 10
  bem_mutex <- bem_mutex[, as.logical(apply(
    bem_mutex, 2, function(x) sum(as.numeric(x))) >= min_variants)]
  ngs <- make_mutex_df(bem_mutex)
  ngs <- ngs %>% tibble::rownames_to_column("Symbol")
  utils::write.table(
    x = ngs,
    file = paste0(save, "/OncoBird_mutex_ngs_alterations_minx.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  if (mutex_output_exists == FALSE) {
    if (!file.exists(paste0(save, "/OncoBird_mutex_ngs_alterations_minx.csv"))) {
      stop(save,
        "OncoBird_mutex_ngs_alterations_minx.csv does not exist."
      )
    }
    if (!file.exists(paste0(save, "/Network.sif"))) {
      stop(save, "Network.sif does not exist.")
    }
    if (!file.exists(paste0(save, "/parameters.txt"))) {
      stop(save, "parameters.txt does not exist.")
    }

    message("Calculating mutual exclusivity")
    system(paste0("java -jar ", mutex_path, "mutex/target/mutex.jar ", save),
      intern = TRUE)
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
      alts <- unlist(lapply(((ngs_exclusive[i, 2:(ncol(ngs_exclusive))])), as.character))
      if (any(know %in% alts)) {
        newalts <- c(alts[alts != ""], know[!know %in% alts])
        newalts <- c(newalts, rep("", ncol(ngs_exclusive) - 1 - length(newalts)))
        ngs_exclusive[j, ] <- c(1, newalts)
        j <- j + 1
      }
    }
    #
    newalts <- c(alts[alts != ""], know[!know %in% alts])
    newalts <- c(newalts, rep("", ncol(ngs_exclusive) - 1 - length(newalts)))
    ngs_exclusive[j, ] <- c(1, know, rep("", ncol(ngs_exclusive) - 1 - length(know)))
  }
  ngs_exclusive$name <- unlist(lapply(
    seq_len(nrow(ngs_exclusive)),
    function(x) {
      tmp <- as.character(unlist(ngs_exclusive[x, 2:ncol(ngs_exclusive)]))
      tmp <- tmp[tmp != "" & tmp != "NA"]
      return(paste(tmp, collapse = ";"))
    }
  ))

  # remove duplicate modules
  ngs_exclusive <- ngs_exclusive[
    !duplicated(lapply(
      seq_len(nrow(ngs_exclusive)),
      function(x) {
        tmp <- as.character(unlist(ngs_exclusive[
          x,
          2:(ncol(ngs_exclusive) - 1)
        ]))
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
  save(data_mutations_me, file = paste0(save, "/OncoBird_data_mutations_me.RData"))

  # replace all ; in the end with nothing, cleaner for plotting
  ngs_exclusive$name <- gsub(";*$", "", ngs_exclusive$name)

  message("Number of single genes: ", as.character(ncol(data_mut)))
  message(
    "Number of single plus modules: ",
    as.character(ncol(data_mutations_me))
  )

  metadata(se)[["me_modules"]] <- data_mutations_me
  metadata(se)[["ngs_exclusive"]] <- ngs_exclusive

  return(se)
}



########################################################################
#' Plot mutual exclusivity
#'
#' Plot the results of the mutex algorithm
#'
#' @param se SE object
#'
#' @return ggplot object
#' @export
#'
#' @examples \dontrun{
#'   sim_data <- sim_data()
#'
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = c("OS", "DFS"),
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "PatientID",
#'     treatment_column = "Adj",
#'     mutation_columns = colnames(sim_data)[14:ncol(sim_data)]
#'   )
#'
#'   pl_mutual_exclusivity(cl_mututal_exclusivity(se,
#'     min_variants = 10,
#'     mutex_output_exists = TRUE,
#'     save = "metadata/",
#'     mutex_path = "path/to/mutex/installation"
#'   ))
#' }
#'
pl_mutual_exclusivity <- function(se) {
  ngs_exclusive_plot <- metadata(se)$ngs_exclusive
  ngs_exclusive_plot$name <- gsub(";", "/", ngs_exclusive_plot$name)
  ngs_exclusive_plot$identifier <-
    lapply(
      seq_len(nrow(metadata(se)$ngs_exclusive)),
      function(x) {
        tmp <- as.character(unlist(metadata(se)$ngs_exclusive[
          x,
          2:(ncol(ngs_exclusive_plot) - 1)
        ]))
        tmp <- tmp[tmp != ""]
        return(sort(unlist(tmp)))
      }
    )
  ngs_exclusive_plot$Score <- as.numeric(gsub(
    ",", ".",
    gsub("E", "e", as.character(ngs_exclusive_plot[, 1]))
  ))
  ngs_exclusive_plot <-
    ngs_exclusive_plot[!duplicated(ngs_exclusive_plot$identifier), ]
  ggplot(ngs_exclusive_plot) +
    geom_bar(aes(
      x = reorder(name, -log10(Score)),
      y = -log10(Score)
    ), stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    xlab("Mutually exclusive module") +
    ylab("Exclusivity score") +
    theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm"))
}

########################################################################
#' Calculate treatment-specific biomarkers
#'
#' Calculate biomarkers for each readout, subtype and treatment arm
#'
#' @param se SE object
#' @param include_covariates character vector with confounding factors
#' @param min_samples minimum number of mutants for a mutation to be
#' included in statistical tests
#' @param min_redistribution minimum number of redistributed tumours per
#' gene module for being included
#' @param treatment character for treatment arm column
#' @param subtypes character for subtype column
#' @param readouts character for readout (OS,PFS,DFS,RFS,ORR)
#'
#' @return SE object with metadata containing treatment-specific statistical
#' tests 'treatment_specific_biomarkers'
#' @export
#'
#' @examples {
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#' }
#'
cl_treatment_specific_biomarkers <- function(se,
                                             include_covariates = NULL,
                                             min_samples = 10,
                                             min_redistribution = 3,
                                             treatment = NULL,
                                             subtypes = NULL,
                                             readouts = c("OS")) {

  # Check input

  data_clinical <- as.data.frame(rowData(se))

  if (!any(subtypes == "NGS.probe")) {
    subtypes <- c("NGS.probe", subtypes)
  }

  conditions <- create_conditions(
    treatment = treatment,
    strat = subtypes,
    survival = readouts
  )

  assertthat::assert_that("treatment" %in% colnames(data_clinical))
  assertthat::assert_that(is.character(treatment))
  assertthat::assert_that(is.character(subtypes))
  assertthat::assert_that(is.character(readouts))

  for (i in subtypes) {
    assertthat::assert_that(i %in% colnames(data_clinical))
  }

  if (any(is.null(include_covariates))) {
    include_covariates <- c("1")
    message("Using default c(1)")
  } else {
    if (all(include_covariates %in% colnames(data_clinical))) {
      message("\nUsing covariates: ", include_covariates)
    } else {
      stop("covariates needs at least one match with colnames")
    }
  }

  # Update se metadata
  workflow_metadata <- metadata(se)$wf_meta
  workflow_metadata[["treatment"]] <- treatment
  workflow_metadata[["subtypes"]] <- subtypes
  workflow_metadata[["readouts"]] <- readouts
  metadata(se)$wf_meta <- workflow_metadata

  feat <- metadata(se)$me_modules
  if (all(is.null(feat))) {
    warning("Mutually exclusive modules were not analyzed yet..")
    feat <- assay(se) %>% as.data.frame()
  }
  treatment_specific_biomarkers <- list()
  for (i in seq_len(nrow(conditions))) {
    treatment_specific_biomarkers[[i]] <- calc_asso_cox(
      response = data_clinical,
      features = feat,
      survival = conditions$survival[i],
      strat = conditions$strat[i],
      treatment = conditions$treatment[i],
      covariates = include_covariates,
      N = min_samples,
      min_redistribution = min_redistribution
    )
  }
  treatment_specific_biomarkers <- do.call(rbind, treatment_specific_biomarkers)

  metadata(se)[["treatment_specific_biomarkers"]] <- treatment_specific_biomarkers

  return(se)
}



########################################################################
#' Calculate predictive biomarkers
#'
#' Calculate biomarkers for each readout and subtype using interaction tests
#'
#' @param se SE object
#' @param include_covariates character vector with confounding factors
#' @param min_samples minimum number of mutants for a mutation to be
#' included in statistical tests
#' @param min_redistribution minimum number of redistributed tumours per
#' gene module for being included
#' @param compare_treatments logical with default 'True' (deprecated)
#'
#' @return SE object with metadata containing treatment-specific statistical
#' tests 'predictive_biomarkers'
#' @export
#'
#' @examples {
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'     
#' }
#'
cl_predictive_biomarkers <- function(se,
                                     include_covariates = NULL,
                                     min_samples = 10,
                                     min_redistribution = 3,
                                     compare_treatments = TRUE) {
  data_clinical <- as.data.frame(rowData(se))
  workflow_metadata <- metadata(se)$wf_meta
  treatment <- workflow_metadata[["treatment"]]
  subtypes <- workflow_metadata[["subtypes"]]
  readouts <- workflow_metadata[["readouts"]]


  # Fix bug, only interaction between two treatments possible
  data_clinical <- data_clinical[as.character(data_clinical$treatment) %in%
    treatment, ]
  if (length(levels(factor(data_clinical$treatment))) != 2) {
    stop("Number of treatment arms is not two !")
  }

  conditions <- create_conditions(
    treatment = treatment,
    strat = subtypes,
    survival = readouts
  )

  if (any(is.null(include_covariates))) {
    include_covariates <- c("1")
    message("Using default c(1)")
  } else {
    if (all(include_covariates %in% colnames(data_clinical))) {
      message("\nUsing covariates: ", include_covariates)
    } else {
      stop("covariates needs at least one match with colnames")
    }
  }

  feat <- metadata(se)$me_modules
  if (all(is.null(feat))) {
    warning("Mutually exclusive modules were not analyzed yet..")
    feat <- assay(se) %>% as.data.frame()
  }
  # Make sure selected treatment is the only factor levels of df
  data_clinical$treatment <- factor(data_clinical$treatment)


  predictive_biomarkers <- list()
  for (i in seq_len(nrow(conditions))) {
    # linear models for treatment interactions
    predictive_biomarkers[[i]] <- calc_asso_cox(
      response = data_clinical,
      features = feat,
      survival = conditions$survival[i],
      strat = conditions$strat[i],
      treatment = conditions$treatment[i],
      covariates = include_covariates,
      N = min_samples,
      compare_treatments = compare_treatments,
      min_redistribution = min_redistribution
    )
    if (nrow(predictive_biomarkers[[i]]) > 0) {
      predictive_biomarkers[[i]]$treatment <- "yes"
    } else {
      message("no pred biomarkers found for this condition")
    }
  }

  predictive_biomarkers <- do.call(rbind, predictive_biomarkers) %>%
    dplyr::filter(!duplicated(paste(survival, strat, cfe)))

  #
  #predictive_biomarkers$fdr <- predictive_biomarkers$fdr_int
  predictive_biomarkers$effectsize_old <- predictive_biomarkers$effectsize
  predictive_biomarkers$effectsize <- predictive_biomarkers$effInt

  metadata(se)[["predictive_biomarkers"]] <- predictive_biomarkers
  return(se)
}



########################################################################
#' Plot treatment-specific biomarkers
#'
#' Plots a volcano plot of treatment-specific biomarkers  for each
#' readout and treatment arm
#'
#' @param se SE object
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5]
#'     )
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS")
#'      )
#'
#'   pl_treatment_specific_biomarkers(se)
#' }
#'
pl_treatment_specific_biomarkers <- function(se) {
  workflow_metadata <- metadata(se)$wf_meta
  treatment <- workflow_metadata[["treatment"]]
  readouts <- workflow_metadata[["readouts"]]


  treatment_specific_biomarkers <-
    metadata(se)$treatment_specific_biomarkers %>%
    mutate(genelabels = fix_annotation(get_annotation(gene))) %>%
    mutate(annotation = paste(genelabels, treatment, survival, sep = "-")) %>%
    mutate(color = paste(strat, treatment, sep = "-"))

  # Plot treatment specific somatic biomarker
  conditions <- expand.grid(list(
    treatment = treatment,
    survival = readouts
  ))

  plot <- list()
  treatment_specific_biomarkers$treatment <- as.character(
    treatment_specific_biomarkers$treatment
  )

  for (i in seq_len(nrow(conditions))) {
    plot[[i]] <- plot_treatment_specific_biomarkers(
      df = treatment_specific_biomarkers %>%
        dplyr::filter(survival == conditions$survival[i] &
          strat == "available" & treatment == conditions$treatment[i]),
      color_highlight_con = 0.3,
      text_annotation = "genelabels",
      y.axis = "P.Value",
      nudge_pos = 1,
      nudge_neg = -1,
      text_size = 3,
      fdr_label_x = 2.4,
      max_text_annotation = 50,
      nudge_y = 2,
      color = "color"
    ) +
      scale_color_manual(values = "brown") +
      theme(legend.position = "bottom") +
      guides(color = FALSE) +
      labs(size = "Amount of altered patients", color = NULL) +
      labs(x = "log10(hazard ratio) - effect size") +
      xlim(c(-2.5, 2.5))
  }
  return(plot)
}



########################################################################
#' Plot treatment-specific and subtype-specific biomarkers
#'
#' Plots a forest plot of treatment-specific and subtype-specific biomarkers
#' for each readout and treatment arm
#'
#' @param se SE object
#' @param fdr_max FDR cutoff
#' @param colors list of color specifications
#' @param labels list of label specifications
#'
#' @return SE object with metadata 'plot$treatment_specific_biomarkers_subtype'
#' as ggplot object
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#' }
#'
pl_treatment_specific_biomarkers_subtype <- function(se,
                                                     fdr_max = 0.3,
                                                     colors = NULL,
                                                     labels = NULL) {
  readouts <- metadata(se)$wf_meta[["readouts"]]

  treatment_specific_biomarkers <-
    metadata(se)$treatment_specific_biomarkers %>%
    mutate(identifier = paste(survival, treatment, sep = "-"))

  treatment_specific_biomarkers_forest <- make_forest(
    treatment_specific_biomarkers
  )
  metadata(se)$treatment_specific_biomarkers_forest <- treatment_specific_biomarkers_forest

  plots <- list()
  for (j in readouts) {
    plots[[j]] <- plot_forest_treatment_specific_biomarkers(
      forest = treatment_specific_biomarkers_forest %>%
        dplyr::filter(survival == j),
      identifier = "identifier", colors = colors, labels = labels,
      fdr_threshold = fdr_max
    ) + theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm"))
  }

  metadata(se)$plot$treatment_specific_biomarkers_subtype <- plots
 # print(plots)
  return(se)
}



########################################################################
#' Calculate module visualisations for treatment-specific biomarkers
#'
#' Calculate forest plots, heatmaps and oncoprints for
#' treatment-specific biomarkers
#'
#' @param se SE object
#' @param p_value p-value threshold
#' @param fdr FDR threshold
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   metadata(se)$tsb_modules_oncoprint$OS$`chemo/resistant`
#' }
#'
cl_tsb_modules_oncoprint <- function(se,
                                     p_value,
                                     fdr) {
  readouts <- metadata(se)$wf_meta[["readouts"]]
  treatment_specific_biomarkers_modules_oncoprint <- list()

  feat <- metadata(se)$me_modules
  if (all(is.null(feat))) {
    warning("Mutually exclusive modules were not analyzed yet..")
    feat <- assay(se) %>% as.data.frame()
  }

  for (i in readouts) {
    treatment_specific_biomarkers_modules_oncoprint[[i]] <- plot_modules_oncoprint_treatment_specific_biomarkers(
      assos = metadata(se)$treatment_specific_biomarkers,
      strat = "available",
      survival = i,
      bem = feat,
      p_value = p_value, # works for smaller 0.1
      fdr = fdr
    )
  }
  metadata(se)$tsb_modules_oncoprint <- treatment_specific_biomarkers_modules_oncoprint
  return(se)
}



########################################################################
#' Plot predictive biomarkers
#'
#' Plot biomarkers for each readout across all tumours using
#' interaction tests in forest plots
#'
#' @param se SE object
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples {
#'
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#' 
#'   pl_predictive_biomarkers(se)
#'   
#' }
#'
pl_predictive_biomarkers <- function(se) {
  readouts <- metadata(se)$wf_meta[["readouts"]]

  predictive_biomarkers <-
    metadata(se)$predictive_biomarkers %>%
    mutate(genelabels = fix_annotation(get_annotation(gene))) %>%
    mutate(annotation = paste(genelabels, treatment, survival, sep = "-")) %>%
    mutate(color = paste(strat, treatment, sep = "-"))

  # Plot treatment specific somatic biomarker
  conditions <- expand.grid(list(
    treatment = c("yes"),
    survival = readouts
  ))

  plot <- list()
  for (i in seq_len(nrow(conditions))) {
    predictive_biomarkers$fdr <- predictive_biomarkers$fdr_int
    plot[[i]] <- plot_treatment_specific_biomarkers(
      df = predictive_biomarkers %>%
        dplyr::filter(survival == as.character(conditions$survival[i]) &
          strat == "available" &
          treatment == as.character(conditions$treatment[i])),
      color_highlight_con = 0.6,
      text_annotation = "genelabels",
      y.axis = "pInt",
      nudge_pos = 1.5,
      nudge_neg = -1,
      text_size = 3,
      fdr_label_x = 2.4,
      max_text_annotation = 50,
      nudge_y = 0.2,
      color = "color"
    ) +
      scale_color_manual(values = "brown") +
      theme(legend.position = "bottom") +
      guides(color = FALSE) +
      labs(size = "Amount of altered patients", color = NULL) +
      labs(x = "log10(hazard ratio) - effect size") #+
    xlim(c(-2.5, 2.5))
  }
  return(plot)
}



########################################################################
#' Plot subtype-specific predictive biomarkers
#'
#' Plot biomarkers for each readout and subtype using interaction
#' tests in forest plots
#'
#' @param se SE object
#' @param fdr_max FDR cutoff
#' @param colors list of color specifications
#' @param labels list of label specifications
#'
#' @return SE object with metadata 'plot$predictive_biomarkers_subtypes'
#' as ggplot object
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#' }
#'
pl_predictive_biomarkers_subtypes <- function(se,
                                              fdr_max = 0.6,
                                              colors = NULL,
                                              labels = NULL) {
  readouts <- metadata(se)$wf_meta[["readouts"]]

  predictive_biomarkers <-
    metadata(se)$predictive_biomarkers %>%
    mutate(identifier = paste(survival, treatment, sep = "-"))

  is_surv <- readouts %in% predictive_biomarkers$survival

  predictive_biomarkers_forest <- make_forest_interaction(predictive_biomarkers)
  metadata(se)$predictive_biomarkers_forest <- predictive_biomarkers_forest

  plots <- list()
  for (j in readouts[is_surv]) {
    plots[[j]] <- plot_forest_treatment_specific_biomarkers(
      forest = predictive_biomarkers_forest %>% dplyr::filter(survival == j),
      identifier = "identifier", colors = colors, labels = labels,
      fdr_column = "fdr_int",
      fdr_threshold = fdr_max
    ) + theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm"))
  }

  metadata(se)$plot$predictive_biomarkers_subtypes <- plots
  #print(plots)
  return(se)
}



########################################################################
#' Calculate treatment comparisons for predictive biomarkers
#'
#' Calculate treatment comparison tests for the significant interaction tests
#'
#' @param se SE object
#' @param subtypes character vector of subtype columns
#' @param readouts character vector of readouts
#' @param covariates character vector of confounding factors
#' @param fdr_i FDR cutoff for interaction tests
#'
#' @return SE object with metadata containing statistical tests in
#' 'pred_comparison' and 'conditions'
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#' }
#'
cl_predictive_comparison <- function(se,
                                     subtypes,
                                     readouts,
                                     covariates,
                                     fdr_i = 0.6) {
  readouts <- metadata(se)$wf_meta[["readouts"]]
  subtypes <- metadata(se)$wf_meta[["subtypes"]]

  subtype_column <- c()
  subtype <- c()
  for (i in seq_len(length(subtypes))) {
    len <- na.omit(unique(rowData(se)[subtypes[i]]))
    names <- as.character(len[, 1])
    subtype_column <- c(subtype_column, rep(subtypes[i], nrow(len)))
    subtype <- c(subtype, names)
  }

  if (!length(subtype) == length(subtype_column)) {
    stop("subtype and subtype_column need to have same length")
  }

  conditions <- reshape::expand.grid.df(
    data.frame(
      subtype = subtype,
      subtype_column = subtype_column
    ),
    data.frame(survival = readouts)
  )


  feat <- metadata(se)$me_modules
  if (all(is.null(feat))) {
    warning("Mutually exclusive modules were not analyzed yet..")
    feat <- assay(se) %>% as.data.frame()
  }


  comparisons <- list()
  # Check each subtype
  for (i in seq_len(nrow(conditions))) {
    message("Subtype: ", conditions$subtype[i])
    rowData(se)$costum <- as.character(rowData(se)[, as.character(conditions$subtype_column[i])]) == conditions$subtype[i]
    tmp <- compare_treatments(
      asso_cox_cet = metadata(se)$predictive_biomarkers %>%
        dplyr::filter(survival == conditions$survival[i] &
          strat == conditions$subtype[i] & fdr_int < fdr_i),
      asso_cox_bev = metadata(se)$predictive_biomarkers %>%
        dplyr::filter(survival == conditions$survival[i] &
          strat == conditions$subtype[i] & fdr_int < fdr_i),
      cms = TRUE,
      readout = conditions$survival[i],
      sav_cox_fig5 = rowData(se),
      bem_cox_fig4_me = feat,
      strat_variable = "costum",
      name = conditions$subtype[i],
      add_row = FALSE,
      covariates = covariates,
      fdr_threshold = 1, # for the selection of prognostic markers to test, deprecated
      verbose = FALSE
    )
    comparisons[[i]] <- tmp
  }

  metadata(se)[["pred_comparison"]] <- comparisons
  metadata(se)[["conditions"]] <- conditions

  return(se)
}



########################################################################
#' Plot treatment comparisons for predictive biomarkers
#'
#' Plot treatment comparison tests in forest plots for the significant
#' interaction tests
#'
#' @param se SE object
#' @param min_altered minimum number of mutant samples to be
#' tested for a treatment comparison
#' @param max_fdr FDR cutoff for subtype-specific tests
#' @param colors list of color specifications
#' @param labels list of label specifications
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#'     
#'   pl_predictive_comparison(se)
#' }
#'
pl_predictive_comparison <- function(se,
                                     min_altered = 10, # default
                                     max_fdr = 0.3, # default
                                     colors = NULL,
                                     labels = NULL) {
  conditions <- metadata(se)$conditions
  readouts <- unique(conditions$survival)
  predictive_comparison <- metadata(se)$pred_comparison

  treatment_specific_biomarkers <- metadata(se)$treatment_specific_biomarkers
  forest_comparison <- list()
  for (i in readouts) {
    forest_comparison[[i]] <- tryCatch(
      plot_forest_comparison(
        comp = predictive_comparison[which(conditions$survival == i)],
        min_altered = min_altered,
        treatment_biomarkers = treatment_specific_biomarkers[treatment_specific_biomarkers$survival == i, ],
        fdr_biomarker_threshold = max_fdr,
        colors = colors,
        labels = labels
      ),
      error = function(e) NULL
    )
    # plotting function is not called anymore
    # if (!is.null(forest_comparison[[i]])) {
    #   gridExtra::grid.arrange(grobs = list(
    #     forest_comparison[[i]]$p_extraction[[1]] +
    #       theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm")) +
    #       theme(legend.position = "none") +
    #       ggtitle("mutant"),
    #     forest_comparison[[i]]$p_extraction[[2]] +
    #       theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm")) +
    #       theme(legend.position = "right") +
    #       ggtitle("wild type")
    #   ), layout_matrix = rbind(c(1, 2)), top = i)
    # } else {
    #   message(
    #     "Predictive biomarkers for ", i,
    #     " are insignificant, increase theshold.."
    #   )
    # }
  }
  return(forest_comparison)
}


########################################################################
#' Plot examples
#'
#' Visualise a specific example for the presented analysis steps
#'
#' @param se SE object
#' @param mutations character vector of mutations
#' @param readout character readout
#' @param subtype_column character subtype name
#' @param subtype character subtype level name
#' @param treatment character treatment level
#' @param covariates character vector of confounding factors
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#'   
#'   pl_example(se,
#'     mutations = "V4_SV",
#'     readout = "OS",
#'     subtype_column = "r",
#'     subtype = "r_cond_1",
#'     treatment = c("chemo", "no_chemo"))
#' }
#'
pl_example <- function(se,
                       mutations,
                       readout,
                       subtype_column,
                       subtype,
                       treatment = NULL,
                       covariates = c(1)) {
  if (readout == "OS") {
    lab <- "Overall survival"
  }
  if (readout == "PFS") {
    lab <- "Progression-free survival"
  }
  if (readout == "DFS") {
    lab <- "Disease-free survival"
  }
  if (readout == "RFS") {
    lab <- "Relapse-free survival"
  }

  bem <- metadata(se)$me_modules

  if (all(is.null(bem))) {
    warning("Mutually exclusive modules were not analyzed yet..")
    bem <- assay(se) %>% as.data.frame()
  }

  clin <- rowData(se)

  if (!any(is.null(treatment))) {
    clin <- clin[clin$treatment %in% treatment, ]
    clin$treatment <- factor(clin$treatment)

    int <- intersect(clin$sample, bem$sample)
    clin <- clin[clin$sample %in% int, ]
    bem <- bem[bem$sample %in% int, ]
  }

  assertthat::assert_that(all(treatment %in% clin$treatment))
  assertthat::assert_that(readout %in% c("OS", "PFS", "DFS", "ORR", "RFS"))
  assertthat::assert_that(subtype %in% as.character(clin[, subtype_column]))
  assertthat::assert_that(subtype_column %in% colnames(clin))
  step <- plot_example(
    bem = bem,
    clin = clin,
    subtype_column = subtype_column,
    subtype = subtype,
    alteration = mutations,
    survival = readout, covariates = covariates
  )
  pt1 <- step$specific[[1]] + ylab(lab)
  pt2 <- step$specific[[2]] + ylab(lab)

  pa1 <- step$comparison$altered + ylab(lab)
  pa2 <- step$comparison$`wild type` + ylab(lab)
  
  pall <- step$all + ylab(lab)

  return(list(pt1, pt2, pa1, pa2, pall))
}



#' Plot oncoprint
#'
#' Plot oncoprint for interesting mutations
#'
#' @param se SE object
#' @param min_mutants mininum number of mutated tumours per mutation
#' to be included in the oncoprint
#' @param remove_mutations character vector of mutations to ignore for
#' oncoprint
#' @return None
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = c("OS", "DFS"),
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "PatientID",
#'     treatment_column = "Adj",
#'     mutation_columns = colnames(sim_data)[14:ncol(sim_data)]
#'   )
#'
#'   cl_subtype_enrichment(
#'     se = se,
#'     col_label = "Smoking.history",
#'     row_label = "EGFR.subtype",
#'     digits = 3
#'   )
#'
#'   se <- cl_enrichment_genomics(se,
#'     sample_column = "sample",
#'     min_mutants = 10,
#'     subtype = c("Smoking.history", "EGFR.subtype")
#'   )
#'   
#'   pl_oncoprint(se)
#' }
#'
pl_oncoprint <- function(se,
                         min_mutants = 5,
                         remove_mutations = NULL) {
  data_mutations <- assay(se) %>% as.data.frame()
  if (!any(is.null(remove_mutations))) {
    data_mutations %>% dplyr::select(-remove_mutations)
  }
  bem_oncoprint <- data_mutations %>% dplyr::select(-sample)
  bem_oncoprint <- bem_oncoprint[, order(apply(
    bem_oncoprint, 2,
    function(x) sum(as.numeric(x))
  ), decreasing = TRUE)]

  bem_oncoprint_type <- lapply(
    colnames(bem_oncoprint),
    function(x) strsplit(x, "_")[[1]][2]
  ) %>% unlist()
  bem_oncoprint_gene <- lapply(
    colnames(bem_oncoprint),
    function(x) strsplit(x, "_")[[1]][1]
  ) %>% unlist()
  for (x in seq_len(length(bem_oncoprint_type))) {
    bem_oncoprint[, x][bem_oncoprint[, x] == "1"] <- bem_oncoprint_type[x]
  }

  anno <- bem_oncoprint %>% t()
  anno[anno == "SV"] <- "MUT"
  anno[anno == "LOSS"] <- "LOSS"
  anno[anno == "AMP"] <- "AMP"
  anno[anno == "0"] <- ""
  anno[!(anno %in% c("MUT", "LOSS", "AMP", ""))] <- "FEAT"
  anno <- anno %>% as.data.frame()
  anno$gene <- bem_oncoprint_gene

  anno_long <- tidyr::pivot_longer(
    data = anno,
    colnames(anno)[-length(colnames(anno))]
  )
  anno_long <- anno_long %>%
    group_by(name, gene) %>%
    summarise(alterations = paste(unique(value), collapse = ";"))
  anno <- tidyr::pivot_wider(
    data = anno_long, names_from = gene,
    values_from = alterations, id_cols = name
  ) %>%
    tibble::column_to_rownames("name")
  # reduce to alteration with more than 12 patients
  anno <- anno[, apply(
    anno, 2,
    function(x) length(which(x != ""))
  ) > min_mutants]
  oncoprint <- ComplexHeatmap::oncoPrint(anno %>% t(),
    alter_fun = alter_fun,
    get_type = function(x) strsplit(x, ";")[[1]],
    col = c(
      "MUT" = "grey30", "AMP" = "darkred",
      "LOSS" = "darkblue", "NA" = "white", "FEAT" = "darkorange"
    ),
    heatmap_legend_param = list(
      title = "Alternations",
      at = c("MUT", "AMP", "LOSS", "NA", "FEAT"),
      labels = c("Loss", "Amplification", "Mutation", "NA", "other")
    )
  )
  oncoprint
}



#' Prepare datasets for run
#'
#' Initiate SE object for downstream analysis
#'
#' @param data dataframe containing clinical and mutational variables
#' @param data_clinical alternatively to 'data', dataframe of clinical
#' variables
#' @param data_mutations alternatively to 'data', dataframe of mutational
#' variables
#' @param vars character vector of clinical endpoints to be analysed
#' (OS/PFS/ORR/RFS/DFS)
#' @param med_impute_muts logical if missing clinical data is imputed by
#' the median
#' @param remove_clin character vector of clinical variables to be excluded
#' from the downstream analysis
#' @param treatment_column character to specify the treatment column
#' @param sample_column character to specify the sample identifier column
#' @param mutation_columns character vector to specifiy the mutational column,
#' if only 'data' is supplied
#'
#' @return SE object
#' @export
#'
#' @examples {
#'   sim_data <- sim_data()
#'
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = c("OS", "DFS"),
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "PatientID",
#'     treatment_column = "Adj",
#'     mutation_columns = colnames(sim_data)[14:ncol(sim_data)]
#'   )
#' }
#'
prepare_data <- function(data = NULL,
                         data_clinical = NULL,
                         data_mutations = NULL,
                         vars = c("OS", "PFS", "RFS", "DFS", "ORR"),
                         med_impute_muts = TRUE,
                         remove_clin = NA,
                         treatment_column = "treatment",
                         sample_column = "sample",
                         mutation_columns = NULL) {
  data <- as.data.frame(data)
  if (is.null(data_mutations) & is.null(data_clinical)) {

    # Check if mutations supplied
    if (is.null(mutation_columns)) {
      stop("No mutation columns supplied !")
    }

    # Check for sample and treatment column in data
    if (!(sample_column %in% colnames(data))) {
      stop("No sample column in data !\n")
    }
    if (!(treatment_column %in% colnames(data))) {
      stop("No treatment column in data !\n")
    }
    data[sample_column] <- make.names(data[, sample_column])
    colnames(data)[which(colnames(data) == sample_column)] <- "sample"
    data[treatment_column] <- make.names(data[, treatment_column])
    colnames(data)[which(colnames(data) == treatment_column)] <- "treatment"

    # Assign mutational data
    data_clinical <- as.data.frame(data[!colnames(data) %in% mutation_columns])
    data_mutations <- as.data.frame(cbind(
      data[, "sample", drop = FALSE],
      data[colnames(data) %in% mutation_columns]
    ))
  } else {

    # Check for sample and treatment column
    if (!(sample_column %in% colnames(data_mutations))) {
      stop("No sample column in mutational data !\n")
    }
    if (!(sample_column %in% colnames(data_clinical))) {
      stop("No sample column in clinical data !\n")
    }
    if (!(treatment_column %in% colnames(data_clinical))) {
      stop("No treatment column in clinical data !\n")
    }
    data_clinical[, sample_column] <- make.names(data_clinical[, sample_column])
    colnames(data_clinical)[which(colnames(data_clinical) == sample_column)] <- "sample"
    data_mutations[, sample_column] <- make.names(data_mutations[, sample_column])
    colnames(data_mutations)[which(colnames(data_mutations) == sample_column)] <- "sample"
    data_clinical[, treatment_column] <- make.names(data_clinical[, treatment_column])
    colnames(data_clinical)[which(colnames(data_clinical) == treatment_column)] <- "treatment"
  }

  # Convert mutations to characters
  data_mutations <- data_mutations %>%
    mutate_all(as.character)
  
  # Renaming mutations
  data_mutations <- data_mutations %>%
    select("sample", everything())
  gene <- unlist(lapply(
    colnames(data_mutations),
    function(x) strsplit(x, "_")[[1]][1]
  ))
  other <- (lapply(
    colnames(data_mutations),
    function(x) strsplit(x, "_")[[1]]
  ))
  # map names to the mutations types
  other <- unlist(lapply(other, function(x) {
    return <- "_SV"
    if (any(grepl("GAIN", toupper(x))) | any(grepl("AMP", toupper(x)))) {
      return <- "_AMP"
    }
    if (any(grepl("DEL", toupper(x))) | any(grepl("LOSS", toupper(x)))) {
      return <- "_LOSS"
    }
    if (any(grepl("CNV", toupper(x))) & (return == "_SV")) {
      warning("Unspecificed CNV appeared in mutation names,
              assuming it is a deletion....")
      return <- "_LOSS"
    }
    return(return)
  }))

  if (gene[1] == "sample") {
    other[1] <- ""
  } else {
    stop("No sample column ...")
  }
  name <- paste0(gene, sep = "", other)
  colnames(data_mutations) <- name

  # summarize duplicates
  dups <- name[duplicated(name)]
  for (dup in dups) {
    dup_df <- data_mutations[, colnames(data_mutations) %in% dup]
    data_mutations <- data_mutations[, !colnames(data_mutations) %in% dup]
    data_mutations[dup] <- apply(
      dup_df, 1,
      function(x) as.character(as.numeric(any(as.numeric(x) == 1)))
    )
  }

  # remove non-binaries
  binary <- as.logical(apply(
    data_mutations, 2,
    function(x) {
      (("1" %in% unique(x)) |
         ("0" %in% unique(x)))
    }
  ) |
    (colnames(data_mutations) %in% c("sample")))
  
  if(length(which(!binary))>0){
    mut_factor_message <- paste0("Columns of somatic alterations ",
          paste(colnames(data_mutations)[!binary], collapse = " / "),
          "contain non-binary ('0','1') values and are removed...")
    message(mut_factor_message)
    message("If those are intended to be kept, please reformat...")
  }
  data_mutations <- data_mutations[
    , binary
  ]
    
  # make the unstratified column
  data_clinical$NGS.probe <- "available"

  # make treatment a factor
  data_clinical$treatment <- factor(data_clinical$treatment)


  # make intersection between data
  int <- intersect(data_clinical$sample, data_mutations$sample)
  data_clinical <- data_clinical[data_clinical$sample %in% int, ] %>%
    arrange(sample)
  data_mutations <- data_mutations[data_mutations$sample %in% int, ] %>%
    arrange(sample)
  no_overlap <- paste0(as.character(length(int)), " samples are overlapping !\n")
  message(no_overlap)

  # Check names
  colnames(data_mutations) <- make.names(colnames(data_mutations))
  colnames(data_clinical) <- make.names(gsub("_", ".", colnames(data_clinical)))


  # get survival or response data
  readout <- c()
  for (var in vars) {
    resp <- unlist(lapply(
      colnames(data_clinical),
      function(x) (grepl(var, toupper(x)))
    ))
    resp <- colnames(data_clinical)[resp]


    # make the clinical response numeric
    data_clinical <- data_clinical %>%
      mutate_if(colnames(data_clinical) %in% resp, as.numeric)
    resp.event <- apply(
      data_clinical[, resp, drop = FALSE], 2,
      function(x) any(as.numeric(x %>% na.omit()) %in% c(0, 1))
    )
    resp.event <- names(which(resp.event))
    resp.time <- resp[!(resp %in% resp.event)]



    message("Taking ", resp.event, " as ", var, " event !\n")
    readout <- c(readout, resp.event)
    colnames(data_clinical)[which(colnames(data_clinical) == resp.event)] <- paste0(var, ".event")
    readout <- c(readout, resp.time)
    message("Taking ", resp.time, " as ", var, " time !\n")
    colnames(data_clinical)[which(colnames(data_clinical) == resp.time)] <- paste0(var, ".months")

    # Check if survival is given in months
    if (any(colnames(data_clinical) == paste0(var, ".months"))) {
      maxtime <- data_clinical[paste0(var, ".months")] %>%
        na.omit() %>%
        max() / 12
      if (maxtime > 30) {
        message("Time for ", var, " is likely in units of days,
                so converting to months...\n")
        data_clinical[paste0(var, ".months")] <- data_clinical[paste0(
          var,
          ".months"
        )] / 30
      }
    }

    if (var == "ORR") {
      colnames(data_clinical)[which(colnames(data_clinical) == paste0(resp.event, ".event"))] <- paste0(var)
      colnames(data_clinical)[which(colnames(data_clinical) == resp.event)] <- paste0(var)
    }
  }

  # get readouts to ignore for NA filtering
  readout_message <- paste0("Readouts are: ", paste(readout, collapse = " / "))
  message(readout_message)
  readout_names <- readout


  # median impute data mutation NAs
  if (med_impute_muts) {
    imp_muts <- names(which(apply(
      data_mutations, 2,
      function(x) length(which(is.na(x))) > 0
    )))
    imputed_mutation_message <- paste0(
      "Imputed mutations are: ",
      paste(imp_muts, collapse = " / ")
    )
    message(imputed_mutation_message)
    for (idx in imp_muts) {
      withna <- data_mutations[, idx]
      withna[is.na(withna)] <- as.character(ceiling(median(
        na.omit(as.numeric(withna))
      )))
      data_mutations[idx] <- withna
    }
  } else {
    imputed_mutation_message <- paste0("No imputed mutations")
  }

  # Ignore columns with many NA's
  imp_clin <- names(which(apply(
    data_clinical, 2,
    function(x) length(which(is.na(x))) > 0
  )))
  imp_clin <- imp_clin[!(imp_clin %in% readout_names)]
  clin_factor_message <- paste0("Clinical factors with nonzero amount of
                          NAs are: ", paste(imp_clin, collapse = " / "))
  message(clin_factor_message)
  message("Consider removing those beforehand...")


  # remove clins
  if (!any(is.na(remove_clin))) {
    message(paste0(
      "Removing clinical factors: ",
      paste(remove_clin, collapse = " / ")
    ))
    data_clinical <- data_clinical[, !(colnames(data_clinical) %in% remove_clin)]
  }

  # Check NAs
  which.na <- apply(data_mutations, 1, function(x) length(which(is.na(x))) > 0)
  # going from characters to logical vector
  readout <- colnames(data_clinical) %in% readout
  #| apply(data_clinical[,!readout], 1, function(x) length(which(is.na(x)))>0)
  # do not filter clinical data
  which.na <- which.na
  data_mutations <- data_mutations[!which.na, ]
  data_clinical <- data_clinical[!which.na, ]
  no_filtered <- paste0(
    "Filtered ", as.character(length(which(which.na))),
    " samples due to missing values in mutational data...\n"
  )
  message(no_filtered)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(data_mutations = as.matrix(data_mutations)),
    rowData = data_clinical
  )

  # write in metadata
  metadata(se)$shiny_message <- list(
    readouts = readout_message,
    imputed_mut = imputed_mutation_message,
    non_zero_clins = clin_factor_message,
    overlaps = no_overlap,
    filtered = no_filtered
  )
  return(se)
}


#' Make summary table
#'
#' make the summary table from the previous tests
#'
#' @param cond specify all readouts and subtypes
#' @param preds association tests for interactions after getting effect sizes
#'  and confidence intervals
#' @param tsbiomarkers association tests for treatment-specific biomarkers
#' after getting effect sizes and confidence intervals
#' @param clin clinical dataframe
#' @param mut mutational dataframe
#' @param mut_me mutational dataframe containing mutually exclusive modules
#' @param fdr_int_threshold interaction FDR threshold
#' @param fdr_biomarker_threshold teatment-specific FDR threshold
#' @param round_digit rounding digit
#' @param eps numerical tolerance of numerical values
#' @param sample_column column in clinical and mutational data specifying the
#' individual
#' @param covariates confounding factors
#' @param treatment treatment column
#'
#' @return summary table dataframe
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#'     
#'   summary_table <- make_table_summary(
#'     cond = metadata(se)$conditions,
#'     preds = metadata(se)$predictive_biomarkers_forest,
#'     tsbiomarkers = metadata(se)$treatment_specific_biomarkers_forest,
#'     clin = rowData(se),
#'     mut = assay(se),
#'     mut_me = metadata(se)$me_modules,
#'     fdr_int_threshold = 0.4,
#'     fdr_biomarker_threshold = 0.4,
#'     round_digit = 2,
#'     sample_column = "sample",
#'     covariates = c(1),
#'     treatment = metadata(se)$wf_meta$treatment
#' )
#' }
#'
make_table_summary <- function(cond = conditions,
                               preds = predictive_biomarkers_forest,
                               tsbiomarkers = treatment_specific_biomarkers_forest,
                               clin = data_clinical,
                               mut = data_mutations,
                               mut_me = data_mutations_me,
                               fdr_int_threshold = 0.1,
                               fdr_biomarker_threshold = 0.1,
                               round_digit = 2,
                               eps = 0.00001,
                               sample_column = "sample",
                               covariates = NULL,
                               treatment = NULL) {
  if (is.null(preds)) {
    stop("pl_predictive_biomarkers_subtypes has not been run yet...")
  }
  if (is.null(tsbiomarkers)) {
    stop("pl_treatment_specific_biomarkers_subtype has not been run yet...")
  }

  # fixing multiple treatment arms
  if (!is.null(treatment)) {
    del <- which((clin$treatment %in% treatment))
    clin <- clin[del, ]
    mut_me <- mut_me[del, ]
    clin$treatment <- factor(clin$treatment)
  }
  if (length(levels(factor(clin$treatment))) != 2) {
    stop("Not two treatments supplied")
  }

  if (is.null(mut_me)) {
    message("Mutually exclusive modules were not analysed yet...")
    mut_me <- mut %>% as.data.frame()
  }

  comps <- list()

  # Check each subtype
  for (i in seq_len(nrow(cond))) {
    clin$costum <- as.character(clin[, as.character(cond$subtype_column[i])]) == cond$subtype[i]
    tmp <- compare_treatments(
      asso_cox_cet = preds %>%
        dplyr::filter(strat == cond$subtype[i] & fdr_int < fdr_int_threshold),
      asso_cox_bev = preds %>%
        dplyr::filter(strat == cond$subtype[i] & fdr_int < fdr_int_threshold),
      cms = TRUE,
      readout = cond$survival[i],
      sav_cox_fig5 = clin,
      bem_cox_fig4_me = mut_me,
      strat_variable = "costum",
      name = cond$subtype[i],
      sample_column = sample_column,
      covariates = covariates,
      add_row = FALSE,
      fdr_threshold = 1, # for the selection of prognostic markers to test
      verbose = FALSE
    )
    comps[[i]] <- tmp
  }

  for (i in seq_len(length(comps))) {
    if (!is.null(comps[[i]])) {
      comps[[i]]$readout <- cond$survival[i]
    }
  }
  summary <- do.call(rbind, comps)


  # insert filtering for treatment_specific biomarkers
  if (!is.null(tsbiomarkers)) {
    which_biomarkers <- apply(
      tsbiomarkers[tsbiomarkers$fdr < fdr_biomarker_threshold, c("cfe", "strat")],
      1, function(x) paste(x, collapse = "---")
    )
    summary$help <- apply(
      summary[, c("gene", "strat_variable")],
      1, function(x) paste(x, collapse = "---")
    )
    summary <- summary[summary$help %in% unlist(which_biomarkers), , drop = FALSE]
    summary <- summary %>% dplyr::select(-help)
  }
  ####################################################


  summary <- make_forest_comp(df = summary)
  summary <- add_medians(summary)

  summary <- summary[, c(
    "hazard", "upper", "lower", "gene", "alteration",
    "mutantorwildtype", "strat_variable", "readout", "mean1",
    "mean2", "n1", "n2", "p"
  )]
  summary <- summary %>% mutate_at(
    c("hazard", "upper", "lower"),
    function(e) {
      e[e > 10] <- NA
      e[e < (-10)] <- NA
      return(round(e, digits = round_digit))
    }
  )

  # make characters
  summary <- summary %>%
    mutate(hazard = paste0(hazard, " [", upper, ";", lower, "]")) %>%
    dplyr::select(-c(upper, lower))

  summary <- summary %>%
    tidyr::pivot_wider(
      values_from = c("hazard", "mean1", "mean2", "n1", "n2", "p"),
      names_from = c("readout", "mutantorwildtype")
    )
  for (i in colnames(summary)) {
    summary[, i] <- unlist(lapply(
      as.data.frame(summary)[, i],
      function(x) {
        if (is.null(x)) {
          NA
        } else {
          unique(x)
        }
      }
    ))
  }


  preds <- preds[, c("gene", "strat", "pInt", "fdr_int", "survival")] %>%
    tidyr::pivot_wider(
      values_from = c("pInt", "fdr_int"),
      names_from = c("survival")
    ) %>%
    dplyr::rename(strat_variable = strat, gene = gene)
  summary <- merge(preds, summary, by = c("gene", "strat_variable"))

  tsbiomarkers <- tsbiomarkers[, c(
    "gene", "strat", "hazard", "upper",
    "lower", "P.Value", "fdr", "survival", "treatment"
  )] %>%
    mutate_at(c("hazard", "upper", "lower"), function(e) {
      e[e > 10] <- NA
      e[e < (-10)] <- NA
      return(round(e, digits = round_digit))
    }) %>%
    mutate(hazard = paste0(hazard, " [", upper, ";", lower, "]")) %>%
    dplyr::select(-c(upper, lower)) %>%
    tidyr::pivot_wider(
      values_from = c("P.Value", "fdr", "hazard"),
      names_from = c("survival", "treatment")
    ) %>%
    dplyr::rename(strat_variable = strat, gene = gene)
  summary <- merge(summary, tsbiomarkers, by = c("gene", "strat_variable"))

  # change names
  summary$geneanno <- summary$gene
  summary$gene <- fix_annotation(get_annotation(summary$gene))

  summary <- summary %>%
    mutate_if(
      function(x){is.numeric(x) & all(na.omit(x<=1))},
      function(x) format.pval(x, eps = eps, digits = round_digit)
    ) %>%
    dplyr::select(-alteration)
  return(summary)
}


#' Perform permutation-based multiplicity correction of treatment effects
#'
#' Derive multiplicity-corrected treatment effect p-values
#'
#' @param se SE object
#' @param treatment levels for treatment column
#' @param subtype character for subtype column
#' @param readout character for readout (OS,PFS,DFS,RFS,ORR)
#' @param include_covariates  character vector with confounding factors
#' @param min_samples Specifies minimum amount of mutants per alteration
#' @param min_redistribution how many patients get different status when
#' using MEs
#' @param meta_load logical if analysis exists already
#' @param meta_path path of metadata object
#' @param n_permutations number of permutations for multiplicity adjustment
#' @param fdr_int_threshold interaction FDR threshold
#' @param fdr_threshold treatment-specific FDR threshold
#' 
#' @return a table of significant gene modules
#'
#' @export
#'
#' @examples {
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#'     
#'     
#'   permuts <- cl_permutations(
#'     se = se,
#'     treatment = metadata(se)$wf_meta$treatment,
#'     subtype = "r",
#'     readout = c("OS"),
#'     include_covariates = NULL,
#'     min_samples = 10,
#'     min_redistribution = 10,
#'     meta_load = FALSE,
#'     meta_path = "metadata",
#'     n_permutations = 10,
#'     fdr_int_threshold = 0.4,
#'     fdr_threshold = 0.4
#' )
#' 
#' permuts
#' 
#' }
#'
cl_permutations <- function(se,
                            treatment,
                            subtype,
                            readout,
                            include_covariates = NULL,
                            min_samples = 10,
                            min_redistribution = 10,
                            meta_load = TRUE,
                            meta_path = "metadata",
                            n_permutations = 100,
                            fdr_int_threshold = 0.1,
                            fdr_threshold = 0.1) {
  assertthat::assert_that(length(subtype) == 1)
  assertthat::assert_that(length(readout) == 1)
  seh <- se
  if (meta_load) {
    if (file.exists(paste0(meta_path, "/permutations1.RData")) &
      file.exists(paste0(meta_path, "/permutations2.RData")) &
      file.exists(paste0(meta_path, "/permutations3.RData"))) {
      message("Loading metadata for permutations...")
      load(paste0(meta_path, "/permutations1.RData"))
      load(paste0(meta_path, "/permutations2.RData"))
      load(paste0(meta_path, "/permutations3.RData"))
      n_permutations <- length(resamplings1) - 1
    } else {
      stop("Metadata for permutations does not exist !")
    }
  } else {
    resamplings1 <- list()
    resamplings2 <- list()
    resamplings3 <- list()
    pb <- txtProgressBar(min = 0, max = n_permutations + 1, style = 3)
    for (k in seq_len((n_permutations + 1))) {
      seh <- se
      # data_mutations_r <- assay(se)
      data_clinical_r <- rowData(se)
      # set.seed(k)
      setTxtProgressBar(pb, k)
      res <- sample(seq_len(nrow(data_clinical_r)))
      if (k != (n_permutations+1)) {
        data_clinical_r$treatment <- data_clinical_r$treatment[res]
      }
      rowData(seh) <- data_clinical_r

      # Step 1
      seh <- cl_treatment_specific_biomarkers(seh,
        include_covariates = include_covariates,
        min_samples = min_samples,
        min_redistribution = min_redistribution,
        treatment = treatment,
        subtypes = subtype,
        readouts = readout
      )
      resamplings1[[k]] <- metadata(seh)$treatment_specific_biomarkers %>%
        dplyr::select(-c(model, medians))

      # Step 2
      seh <- cl_predictive_biomarkers(seh,
        include_covariates = include_covariates,
        min_samples = min_samples,
        min_redistribution = min_redistribution,
        compare_treatments = TRUE
      )
      resamplings2[[k]] <- metadata(seh)$predictive_biomarkers %>%
        dplyr::select(-c(model, medians))

      # Step 3
      seh <- cl_predictive_comparison(
        se = seh,
        subtypes = subtype,
        readouts = readout,
        covariates = include_covariates,
        fdr_i = 1.1
      )
      resamplings3[[k]] <- do.call(rbind, metadata(seh)$pred_comparison) %>%
        dplyr::select(-c(model, matrix))
    }
    # Save runs
    save(resamplings1, file = paste0(meta_path,"/permutations1.RData"))
    save(resamplings2, file = paste0(meta_path,"/permutations2.RData"))
    save(resamplings3, file = paste0(meta_path,"/permutations3.RData"))
  }

  k <- n_permutations + 1
  bio1 <- resamplings1[[k]][, c("fdr", "gene", "strat")]
  bio2 <- resamplings2[[k]][, c("fdr_int", "gene", "strat")]
  bio3 <- resamplings3[[k]][, c(
    "gene", "strat_variable",
    "treatmentp", "mutantorwildtype", "eff"
  )] %>%
    dplyr::rename(strat = strat_variable)
  bio <- merge(bio1, bio2, all = TRUE, by = c("gene", "strat"))
  bio <- merge(bio, bio3, all = TRUE, by = c("gene", "strat"))

  bioo <- distinct(bio[(bio$fdr < fdr_threshold) &
    (bio$fdr_int < fdr_int_threshold), ])
  bioo$id <- paste(bioo$strat, bioo$gene)

  message("Calculating null distributions...")
  p_null <- c()
  for (k in seq_len(n_permutations)) {
    bio1 <- resamplings1[[k]][, c("fdr", "gene", "strat")]
    bio2 <- resamplings2[[k]][, c("fdr_int", "gene", "strat")]
    bio3 <- resamplings3[[k]][, c(
      "gene", "strat_variable", "treatmentp",
      "mutantorwildtype", "eff"
    )] %>%
      dplyr::rename(strat = strat_variable)
    bio <- merge(bio1, bio2, all = TRUE, by = c("gene", "strat"))
    bio <- merge(bio, bio3, all = TRUE, by = c("gene", "strat"))
    bio$id <- paste(bio$strat, bio$gene)

    bio <- bio[bio$strat %in% na.omit(levels(factor(rowData(seh)[, subtype]))), ] # adding factor()
    bio <- distinct(bio[(bio$fdr < fdr_threshold) &
      (bio$fdr_int < fdr_int_threshold), ])

    p_null[k] <- min(c(bio$treatmentp[which.min(bio$treatmentp)], 1),
      na.rm = TRUE
    )
  }

  bioob <- bioo %>%
    group_by(gene, strat) %>%
    dplyr::slice(which.min(treatmentp))
  bioob$treatmentp_adj <- unlist(lapply(
    bioob$treatmentp,
    function(x) length(which((p_null < x))) / n_permutations
  ))
  bioob <- bioob %>% arrange(desc(eff))

  return(bioob)
}


#' Perform bootstrap for estimation of treatment effect confidence intervals
#'
#' Derive treatment effect adjusted confidence intervals
#'
#' @param se SE object
#' @param treatment levels for treatment column
#' @param subtype character for subtype column
#' @param readout character for readout (OS,PFS,DFS,RFS,ORR)
#' @param include_covariates  character vector with confounding factors
#' @param min_samples Specifies minimum amount of mutants per alteration
#' @param min_redistribution how many patients get different status when
#' using MEs
#' @param meta_load logical if analysis exists already
#' @param meta_path path of metadata object
#' @param n_bootstraps number of permutations for multiplicity adjustment
#' @param fdr_int_threshold interaction FDR threshold
#' @param fdr_threshold treatment-specific FDR threshold
#' 
#' @return a table of significant gene modules
#'
#' @export
#'
#' @examples \dontrun{
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#'     
#'   cl_bootstrap(
#'     se = se,
#'     treatment = metadata(se)$wf_meta$treatment,
#'     subtype = "r",
#'     readout = c("OS"),
#'     include_covariates = NULL,
#'     min_samples = 10,
#'     min_redistribution = 10,
#'     meta_load = FALSE,
#'     meta_path = "metadata",
#'     n_bootstraps = 10,
#'     fdr_int_threshold = 0.1,
#'     fdr_threshold = 0.1
#' )
#' }
#'
cl_bootstrap <- function(se,
                         treatment,
                         subtype,
                         readout,
                         include_covariates = NULL,
                         min_samples = 10,
                         min_redistribution = 10,
                         meta_load = TRUE,
                         meta_path = "metadata",
                         n_bootstraps = 100,
                         fdr_int_threshold = 0.1,
                         fdr_threshold = 0.1) {
  # assertthat::assert_that(length(subtype) == 1)
  assertthat::assert_that(length(readout) == 1)
  seh <- se
  if (meta_load) {
    if (file.exists(paste0(meta_path, "/bootstrap1.RData")) &
      file.exists(paste0(meta_path, "/bootstrap2.RData")) &
      file.exists(paste0(meta_path, "/bootstrap3.RData"))) {
      message("Loading metadata for bootstrap...")
      load(paste0(meta_path, "/bootstrap1.RData"))
      load(paste0(meta_path, "/bootstrap2.RData"))
      load(paste0(meta_path, "/bootstrap3.RData"))
      n_permutations <- length(resamplings1) - 1
    } else {
      stop("Metadata for bootstrap does not exist !")
    }
  } else {
    resamplings1 <- list()
    resamplings2 <- list()
    resamplings3 <- list()
    pb <- txtProgressBar(min = 0, max = n_bootstraps + 1, style = 3)
    for (k in seq_len((n_bootstraps + 1))) {
      # set.seed(k+111)
      setTxtProgressBar(pb, k)
      seh <- se
      data_mutations_r <- assay(se) %>%
        as.data.frame() %>%
        mutate(sample = rowData(se)$sample)
      data_clinical_r <- rowData(se)
      data_mutations_me_r <- metadata(se)$me_modules %>% as.data.frame()
      if (all(is.null(data_mutations_me_r))) {
        warning("Mutually exclusive modules were not analyzed yet..")
        data_mutations_me_r <- assay(se) %>% as.data.frame()
      }

      res <- sample(seq_len(nrow(data_clinical_r)), replace = TRUE)
      if (k != (n_bootstraps + 1)) {
        data_clinical_r <- data_clinical_r[res, ]
        data_clinical_r$sample <- make.names(data_clinical_r$sample,
          unique = TRUE
        )
        data_mutations_r <- data_mutations_r[res, ]
        # data_mutations_r$sample <- make.names(data_mutations_r$sample, unique = T)
        data_mutations_me_r <- data_mutations_me_r[res, ] 
        data_mutations_me_r$sample <- make.names(data_mutations_me_r$sample,
          unique = TRUE
        )
      }
      rowData(seh) <- data_clinical_r
      metadata(seh)$me_modules <- data_mutations_me_r
      assay(seh, withDimnames = FALSE) <- data_mutations_r

      # Step 1
      seh <- cl_treatment_specific_biomarkers(seh,
        include_covariates = include_covariates,
        min_samples = min_samples,
        min_redistribution = min_redistribution,
        treatment = treatment,
        subtypes = subtype,
        readouts = readout
      )
      resamplings1[[k]] <- metadata(seh)$treatment_specific_biomarkers %>%
        dplyr::select(-c(model, medians))

      # Step 2
      seh <- cl_predictive_biomarkers(seh,
        include_covariates = include_covariates,
        min_samples = min_samples,
        min_redistribution = min_redistribution,
        compare_treatments = TRUE
      )
      resamplings2[[k]] <- metadata(seh)$predictive_biomarkers %>%
        dplyr::select(-c(model, medians))

      # Step 3
      seh <- cl_predictive_comparison(
        se = seh,
        subtypes = subtype,
        readouts = readout,
        covariates = include_covariates,
        fdr_i = 1.1
      )
      resamplings3[[k]] <- do.call(rbind, metadata(seh)$pred_comparison) %>%
        dplyr::select(-c(model, matrix))
    }
    # Save runs
    save(resamplings1, file = paste0(meta_path,"/bootstrap1.RData"))
    save(resamplings2, file = paste0(meta_path,"/bootstrap2.RData"))
    save(resamplings3, file = paste0(meta_path,"/bootstrap3.RData"))
  }

  k <- (n_bootstraps + 1)
  bio1 <- resamplings1[[k]][, c("fdr", "gene", "strat")]
  bio2 <- resamplings2[[k]][, c("fdr_int", "gene", "strat")] # ,"effInt"
  bio3 <- resamplings3[[k]][, c(
    "gene", "strat_variable",
    "treatmentp", "mutantorwildtype", "eff"
  )] %>%
    dplyr::rename(strat = strat_variable)
  bio <- merge(bio1, bio2, all = TRUE, by = c("gene", "strat"))
  bio <- merge(bio, bio3, all = TRUE, by = c("gene", "strat"))

  biooo <- distinct(bio)
  biooo$id <- paste(biooo$strat, biooo$gene)

  bioo <- bio

  # Initiate the final selection of biomarkers
  bioob <- bioo[(bioo$fdr < fdr_threshold) &
    (bioo$fdr_int < fdr_int_threshold), ] %>%
    group_by(gene, strat) %>%
    dplyr::slice(which.min(treatmentp)) %>%
    arrange(desc(eff))

  # Process for resampling
  bioo$id <- paste(bioo$strat, bioo$gene)
  bioo <- bioo %>%
    dplyr::select(-fdr) %>%
    distinct() %>%
    na.omit()

  message("Calculating bootstrapped distributions...")
  cis <- list()
  cisb <- list()
  test <- c()
  for (l in seq_len(nrow(bioob))) {
    eff <- (bioob$eff[l])
    trte <- list()
    trteq <- list()
    tid <- list()
    for (k in seq_len(n_bootstraps)) {
      bio1 <- resamplings1[[k]][, c("fdr", "gene", "strat")]
      bio2 <- resamplings2[[k]][, c("fdr_int", "gene", "strat")]
      bio3 <- resamplings3[[k]][, c(
        "gene", "strat_variable",
        "treatmentp", "mutantorwildtype", "eff"
      )] %>%
        dplyr::rename(strat = strat_variable)
      bio <- merge(bio1, bio2, all = TRUE, by = c("gene", "strat"))
      bio <- merge(bio, bio3, all = TRUE, by = c("gene", "strat"))
      bio$id <- paste(bio$strat, bio$gene)
      # filter for associations where a stat. test was done, only small biases
      bio <- bio[bio$id %in% bioo$id, ] %>%
        dplyr::select(-fdr) %>%
        distinct()
      # print(nrow(bio)); test[k] <- min(bio$fdr_int, na.rm = T)
      # bio <- distinct(bio[(bio$fdr < fdr_threshold) &
      # (bio$fdr_int < fdr_int_threshold),])

      # choose same effect sign markers
      tmp_df <- bio[((bio$eff * sign(eff)) > 0), , drop = FALSE] %>% na.omit()
      tmp_df <- tmp_df %>%
        group_by(gene, strat) %>%
        dplyr::slice(which.max(abs(eff))) # dplyr::slice(which.min(treatmentp))
      tmp_df <- bio[bio$strat == bioob$strat[l], ] %>% na.omit()


      ll <- as.numeric(ave(row.names(bioob), bioob$strat, FUN = seq_along))[l]

      if (nrow(tmp_df) < ll) {
        ll <- nrow(tmp_df)
      }

      what_biomarker_o <- (tmp_df %>%
        arrange(fdr_int, treatmentp))[ll, c("id", "mutantorwildtype")]
      what_biomarker <- which(tmp_df$id == what_biomarker_o$id &
        tmp_df$mutantorwildtype == what_biomarker_o$mutantorwildtype)

      what_id_o <- tmp_df$id[what_biomarker]
      what_mt_o <- tmp_df$mutantorwildtype[what_biomarker]
      tmp_df <- tmp_df[what_biomarker, , drop = FALSE]
      tmp <- unique(tmp_df$eff) # should be already unique

      if (length(tmp) == 0) {
        tmp <- 0
      }
      trte[[k]] <- tmp

      tmp_df_o <- biooo[(biooo$id %in% what_id_o) &
        (biooo$mutantorwildtype %in% what_mt_o), ]
      # get HR from new subgroup but with original data
      ohr <- min(unique(tmp_df_o$eff[which.min(tmp_df_o$eff)]))
      if (length(ohr) == 0) {
        ohr <- 0
      }
      if (is.na(ohr)) {
        ohr <- 0
      }
      trteq[[k]] <- ohr

      tid[[k]] <- unique(tmp_df_o$id[which.min(tmp_df_o$eff)])
    }
    B <- (unlist(lapply(
      seq_len((n_bootstraps)),
      function(o) eff - (trte[[o]] - trteq[[o]])
    )))
    cis[[l]] <- c(original = eff, mean = mean(B), quantile(B, c(0.025, 0.975)))
    cisb[[l]] <- list(
      B = B,
      trte = unlist(trte),
      trteq = unlist(trteq),
      tid = unlist(tid)
    )
  }

  bioob$B <- cisb
  bioob$CIs <- cis
  return(bioob)
}


#' Get small test dataset for OncoBird
#'
#' Derive treatment effect adjusted confidence intervals
#'
#' @param path Path to the dummy dataset
#'
#' @export
#'
#' @return If path = NULL then it returns all datasets in extdata,
#' if path = filename then it returns paths to file
#' @examples {
#' 
#'   OncoBird_example()
#'   sim_data <- sim_data()
#'   str(sim_data)
#'   
#' }
OncoBird_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "OncoBird"))
  } else {
    system.file("extdata", path, package = "OncoBird", mustWork = TRUE)
  }
}


#' Simulate Dataset
#' 
#' Simulate a dataset containig mutational and clinical data 
#'
#' @param n.obs number of observations
#' @param n.vars number of mutation variables
#' @param n.subs number of substitutions
#'
#' @return a dataframe which can be provided to PrepareData()
#' @export
#'
#' @examples {
#' 
#' data <- sim_data()
#' 
#' data
#' }
sim_data <- function(n.obs = 400,
                     n.vars = 5, 
                     n.subs = 3
){
  x <- matrix(stats::rbinom(n.obs * n.vars, 1, prob = 0.15), n.obs, n.vars)
  s <- factor(stats::rbinom(n.obs, n.subs -1, 0.5) + 1)
  r <- factor(stats::rbinom(n.obs, n.subs -1, 0.5) + 1)
  # simulate randomized treatment
  trt <- stats::rbinom(n.obs, 1, prob = 0.5)
  
  # simulate delta
  delta <- (0.05 + 1*x[,1]*(s == 1) + 1*x[,2]*(s == 2)) 
  #(0.5 + x[,2] - 0.5 * x[,3] - 1 * x[,11] + 1 * x[,1] * x[,12] )
  
  # simulate main effects g(X)
  xbeta <- x[,5] - 1 * x[,4] + x[,3] + 0.5 * x[,1]
  xbeta <- xbeta + delta * (2 * trt - 1)
  
  
  surv.time <- exp(- xbeta + stats::rnorm(n.obs, sd = 0.5))
  cens.time <- exp(stats::rnorm(n.obs, sd = 3))
  OS  <- pmin(surv.time, cens.time)
  OS.status <- 1 * (surv.time <= cens.time)
  
  data <- as.data.frame(x) %>% 
    mutate(trt = trt) %>% 
    mutate(OS = OS) %>%
    mutate(OS.status = OS.status) %>%
    mutate(s = s) %>%
    mutate(r = r)
  data$ID <- seq.int(nrow(data))
  
  data <- data %>% 
    mutate(r = ifelse(r == 1 , "r_cond_1", 
               ifelse(r == 2 , "r_cond_2", "r_cond_3")))
  
  data <- data %>% 
    mutate(s = ifelse(s == 1 , "s_cond_1", 
               ifelse(s == 2 , "s_cond_2", "s_cond_3")))
  
  data <- data %>% 
    mutate(trt = ifelse(trt == 0 , "chemo", "no_chemo"))
  
  return(data)
}


#' Perform cross validation for performance
#'
#' performance of oncobird based on cross-validation for detecting interactions
#'
#' @param se SE object
#' @param treatment levels for treatment column
#' @param subtype character for subtype column
#' @param readout character for readout (OS,PFS,DFS,RFS,ORR)
#' @param include_covariates  character vector with confounding factors
#' @param min_samples Specifies minimum amount of mutants per alteration
#' @param min_redistribution how many patients get different status when
#' using MEs
#' @param meta_load logical if analysis exists already
#' @param meta_path path of metadata object
#' @param n_bootstraps number of permutations for multiplicity adjustment
#' @param fdr_int_threshold interaction FDR threshold
#' @param fdr_threshold treatment-specific FDR threshold
#' 
#' @return a table of significant gene modules
#'
#' @export
#'
#' @examples \dontrun{
#' 
#'   sim_data <- sim_data()
#'   
#'   se <- prepare_data(
#'     data = sim_data,
#'     vars = "OS",
#'     med_impute_muts = TRUE,
#'     remove_clin = NA,
#'     sample_column = "ID",
#'     treatment_column = "trt",
#'     mutation_columns = colnames(sim_data)[1:5])
#'     
#'     
#'   se <- cl_treatment_specific_biomarkers(se, 
#'      min_samples = 10,
#'      min_redistribution = 3,
#'      treatment = c("chemo", "no_chemo"),
#'      subtypes = c("s", "r"),
#'      readouts = c("OS"))
#'   
#'   se <- pl_treatment_specific_biomarkers_subtype(se, fdr_max = 0.1)
#'   
#'   se <- cl_tsb_modules_oncoprint(se, p_value = 0.05, fdr = 0.1)
#'   
#'   se <- cl_predictive_biomarkers(
#'     se = se,
#'     min_samples = 10,
#'     min_redistribution = 5)
#'   
#'   se <- pl_predictive_biomarkers_subtypes(se,
#'     fdr_max = 0.1,
#'     colors = NULL,
#'     labels = NULL)
#'   
#'   se <- cl_predictive_comparison(
#'     se = se,
#'     subtypes = c("s", "r"),
#'     readouts = c("OS"),
#'     covariates = FALSE,
#'     fdr_i = 0.5)
#'     
#'   cl_cv(
#'     se = se,
#'     treatment = metadata(se)$wf_meta$treatment,
#'     subtype = "r",
#'     readout = c("OS"),
#'     include_covariates = NULL,
#'     min_samples = 10,
#'     min_redistribution = 10,
#'     meta_load = FALSE,
#'     meta_path = "metadata",
#'     n_bootstraps = 10,
#'     fdr_int_threshold = 0.1,
#'     fdr_threshold = 0.1
#' )
#' }
#'
cl_cv <- function(se,
                         treatment,
                         subtype,
                         readout,
                         include_covariates = NULL,
                         min_samples = 10,
                         min_redistribution = 10,
                         meta_load = TRUE,
                         meta_path = "metadata",
                         nfold = 5,
                         initialisations = 5,
                         fdr_int_threshold = 0.1,
                         fdr_threshold = 0.1) {
  # assertthat::assert_that(length(subtype) == 1)
  assertthat::assert_that(length(readout) == 1)
  seh <- se
  if (meta_load) {
    if (file.exists(paste0(meta_path, "/cv1.RData")) &
        file.exists(paste0(meta_path, "/cv2.RData")) &
        file.exists(paste0(meta_path, "/cv3.RData"))) {
      message("Loading metadata for cross validation...")
      load(paste0(meta_path, "/cv1.RData"))
      load(paste0(meta_path, "/cv2.RData"))
      load(paste0(meta_path, "/cv3.RData"))
      load(paste0(meta_path, "/cv_results.RData"))
      #n_permutations <- length(resamplings1) - 1
    } else {
      stop("Metadata for cross validation does not exist !")
    }
  } else {
    resamplings1 <- list()
    resamplings2 <- list()
    resamplings3 <- list()
    results <- list()
    pb <- txtProgressBar(min = 0, max = initialisations + 1, style = 3)
    for (k in seq_len((initialisations))) {
      # set.seed(k+111)
      setTxtProgressBar(pb, k)
      seh <- se
      data_mutations_r <- assay(se) %>%
        as.data.frame() %>%
        mutate(sample = rowData(se)$sample)
      data_clinical_r <- rowData(se)
      data_mutations_me_r <- metadata(se)$me_modules %>% as.data.frame()
      if (all(is.null(data_mutations_me_r))) {
        warning("Mutually exclusive modules were not analyzed yet..")
        data_mutations_me_r <- assay(se) %>% as.data.frame()
      }
      
      res <- sample(seq_len(nrow(data_clinical_r)), replace = FALSE)
      
      
      resamplings1[[k]] <- list()
      resamplings2[[k]] <- list()
      resamplings3[[k]] <- list()
      results[[k]] <- list()
      for (l in seq_len(nfold)){
        data_clinical_r <- data_clinical_r[res, ]
        data_clinical_r$sample <- make.names(data_clinical_r$sample,
                                             unique = TRUE
        )
        data_mutations_r <- data_mutations_r[res, ]
        # data_mutations_r$sample <- make.names(data_mutations_r$sample, unique = T)
        data_mutations_me_r <- data_mutations_me_r[res, ] 
        data_mutations_me_r$sample <- make.names(data_mutations_me_r$sample,
                                                 unique = TRUE
        )
        
        # assign train/test
        fold_ids <- cut(1:nrow(data_mutations_r), breaks = nfold, labels = FALSE)
        data_mutations_r_train <- data_mutations_r[fold_ids != l,, drop = FALSE]
        data_mutations_r_test <- data_mutations_r[fold_ids == l,, drop = FALSE]
        data_mutations_me_r_train <- data_mutations_me_r[fold_ids != l,, drop = FALSE]
        data_mutations_me_r_test <- data_mutations_me_r[fold_ids == l,, drop = FALSE]
        data_clinical_r_train <- data_clinical_r[fold_ids != l,, drop = FALSE]
        data_clinical_r_test <- data_clinical_r[fold_ids == l,, drop = FALSE]
        
      #}
        #rowData(seh) <- data_clinical_r_train
        #metadata(seh)$me_modules <- data_mutations_me_r_train
        #assay(seh, withDimnames = FALSE) <- data_mutations_r_train
        
        seh <- SummarizedExperiment::SummarizedExperiment(
          assays = list(data_mutations = as.matrix(data_mutations_r_train)),
          rowData = data_clinical_r_train
        )
        metadata(seh)$me_modules <- data_mutations_me_r_train
        
        # Step 1
        seh <- cl_treatment_specific_biomarkers(seh,
                                                include_covariates = include_covariates,
                                                min_samples = min_samples,
                                                min_redistribution = min_redistribution,
                                                treatment = treatment,
                                                subtypes = subtype,
                                                readouts = readout
        )
        resamplings1[[k]][[l]] <- metadata(seh)$treatment_specific_biomarkers %>%
          dplyr::select(-c(model, medians))
        
        # Step 2
        seh <- cl_predictive_biomarkers(seh,
                                        include_covariates = include_covariates,
                                        min_samples = min_samples,
                                        min_redistribution = min_redistribution,
                                        compare_treatments = TRUE
        )
        resamplings2[[k]][[l]] <- metadata(seh)$predictive_biomarkers %>%
          dplyr::select(-c(model, medians))
        
        # Step 3
        seh <- cl_predictive_comparison(
          se = seh,
          subtypes = subtype,
          readouts = readout,
          covariates = include_covariates,
          fdr_i = 1.1
        )
        resamplings3[[k]][[l]] <- do.call(rbind, metadata(seh)$pred_comparison) %>%
          dplyr::select(-c(model, matrix))
      # end runs ############################################################  
        
      
      ## get scores from runs
      bio1 <- resamplings1[[k]][[l]][, c("fdr", "gene", "strat")]
      bio2 <- resamplings2[[k]][[l]][, c("fdr_int", "gene", "strat")] # ,"effInt"
      bio3 <- resamplings3[[k]][[l]][, c(
        "gene", "strat_variable",
        "treatmentp", "mutantorwildtype", "eff"
      )] %>%
        dplyr::rename(strat = strat_variable)
      bio <- merge(bio1, bio2, all = TRUE, by = c("gene", "strat"))
      bio <- merge(bio, bio3, all = TRUE, by = c("gene", "strat"))
      
      biooo <- distinct(bio)
      biooo$id <- paste(biooo$strat, biooo$gene)
      
      bioo <- bio  
      # Initiate the final selection of biomarkers
      bioob <- bioo[(bioo$fdr < fdr_threshold) &
                      (bioo$fdr_int < fdr_int_threshold), ] %>%
        group_by(gene, strat) %>%
        dplyr::slice(which.min(treatmentp)) %>%
        arrange(desc(eff))
      
      if(nrow(bioob) == 0){
        df <- rep(TRUE, nrow(data_clinical_r_test))
      }else{
        df <- as.data.frame(matrix(ncol = nrow(bioob), nrow = nrow(data_clinical_r_test)))
        for(j in 1:nrow(bioob)){
          df[,j] <- apply((as.data.frame(data_clinical_r_test[,subtype]) %>% dplyr::mutate_all(as.character)), 1, function(x) bioob[j, "strat"] %in% x) & 
            data_mutations_me_r_test[,unlist(bioob[j,"gene"])] == ifelse(unlist(bioob[j,"eff"]) < 0, as.character(bioob[j,"mutantorwildtype"]), as.character(abs(1-unlist(bioob[j,"mutantorwildtype"]))))
        }; df <- apply(df, 1, any)
      }
      
      if(readout == "OS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(OS.months, OS.event) ~ df*treatment",
          "", #covariates
          sep = ""
        )), data = cbind(data_mutations_me_r_test, data_clinical_r_test))
        trt <- survival::coxph(as.formula(paste(
          "survival::Surv(OS.months, OS.event) ~ treatment",
          "", #covariates
          sep = ""
        )), data = cbind(data_mutations_me_r_test, data_clinical_r_test)[df,])
      }
      if(readout == "PFS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(PFS.months, PFS.event) ~ df*treatment",
          "", #covariates
          sep = ""
        )), data = cbind(data_mutations_me_r_test, data_clinical_r_test)) 
        trt <- survival::coxph(as.formula(paste(
          "survival::Surv(PFS.months, PFS.event) ~ treatment",
          "", #covariates
          sep = ""
        )), data = cbind(data_mutations_me_r_test, data_clinical_r_test)[df,]) 
      }
      if(readout == "DFS") {
        cox <- survival::coxph(as.formula(paste(
          "survival::Surv(OS.months, OS.event) ~ df*treatment",
          "", #covariates
          sep = ""
        )), data = cbind(data_mutations_me_r_test, data_clinical_r_test)) 
        trt <- survival::coxph(as.formula(paste(
          "survival::Surv(DFS.months, DFS.event) ~ treatment",
          "", #covariates
          sep = ""
        )), data = cbind(data_mutations_me_r_test, data_clinical_r_test)[df,]) 
      }
      if (readout == "ORR") {
        cox <- glm(as.formula(paste(
          "ORR ~ df*treatment", "",
          sep = ""
        )),
        data = cbind(data_mutations_me_r_test, data_clinical_r_test), family = binomial(link = "logit")
        )
        trt <- glm(as.formula(paste(
          "ORR ~ treatment", "",
          sep = ""
        )),
        data = cbind(data_mutations_me_r_test, data_clinical_r_test)[df,], family = binomial(link = "logit")
        )
      }
      results[[k]][[l]] <- list()
      results[[k]][[l]]$markers <- bioob
      results[[k]][[l]]$validation <- cox
      results[[k]][[l]]$se <- seh
      results[[k]][[l]]$mut_me_test <- data_mutations_me_r_test
      results[[k]][[l]]$mut_test <- data_mutations_r_test
      results[[k]][[l]]$clin_test <- data_clinical_r_test
        
      }
      # Save runs
      save(resamplings1, file = paste0(meta_path,"/cv1.RData"))
      save(resamplings2, file = paste0(meta_path,"/cv2.RData"))
      save(resamplings3, file = paste0(meta_path,"/cv3.RData"))
      save(results, file = paste0(meta_path,"/cv_results.RData"))
    }
  }
  
  return(results)
}
