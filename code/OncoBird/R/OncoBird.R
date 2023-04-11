#' The OncoBird package
#'
#' @description
#' \code{OncoBird} identifies candidates for predictive biomarkers in oncology
#' trials by leveraging mutually exclusive somatic mutations in tumour subtypes
#' for subgroup analysis. It systematically screens for predictive effects of
#' somatic alterations in pre-defined tumour subtypes. It allows users to gain
#' insights about the biomarker landscape in their trial and explore putative
#' cancer subtypes and how they may be leveraged for precision oncology.
#' In our manuscript, we showcased \code{OncoBird} for a clinical trial for
#' FOLFIRI plus either cetuximab and bevacizumab in metastatic colorectal
#' cancer. Here, we showcase \code{OncoBird} by using a publicly available
#' clinical trial for gefitinib in non-small cell lung cancer.
#'
#' @importFrom dplyr mutate %>% group_by summarise mutate_at mutate_all
#' mutate_if bind_rows arrange filter select everything distinct desc
#' @importFrom stats p.adjust reorder sd na.omit phyper coef qt median
#' as.formula binomial glm model.matrix predict contrasts
#' @importFrom utils write.table read.delim help
#' @importFrom SummarizedExperiment assay assay<- rowData rowData<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom ggplot2 ggplot aes xlab ylab coord_flip labs geom_errorbar
#' geom_tile position_dodge annotate theme_minimal theme coord_cartesian unit
#' scale_color_manual guides scale_fill_manual scale_fill_continuous
#' geom_pointrange element_blank element_text scale_y_continuous element_blank
#' geom_text geom_hline geom_point sym geom_bar xlim ylim ggtitle theme_classic
#' element_rect scale_fill_distiller facet_wrap autoplot
#' @importFrom grDevices rainbow
#' @importFrom graphics points text
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats ave quantile
#' @importFrom data.table fread
#' @importFrom survival Surv
#' @importFrom ggfortify ggbiplot
#'
#' @author Alexander J Ohnmacht \email{ohnmachtalexander@gmail.com} and
#' Michael P Menden \email{michael.menden@helmholtz-muenchen.de}
#'
#' @docType package
#' @name OncoBird
NULL

# Quiet R CMD check
utils::globalVariables(c(
  ".", ".data", "Enrichment", "FDR", "Freq", "Order",
  "P.Value", "Score", "Var1", "abundance", "alteration", "alterations",
  "biomarker", "cfe", "clin", "conditions", "confidence", "data_clinical",
  "data_mutations", "data_mutations_me", "direction", "facet_wrap", "fdr",
  "fdr_int", "fdr_new", "fdr_new_int", "gene", "gene_annotation", "genelabels",
  "group", "hazard", "label", "list_glm", "lower", "marker", "mutantorwildtype",
  "name", "ngs_exclusive", "object", "predictive_biomarkers_forest",
  "status", "strat", "strat_variable", "survival", "treatment",
  "treatment_specific_biomarkers", "treatment_specific_biomarkers_forest",
  "upper", "value", "values", "fdr_old", "effectsize", "model", "medians",
  "treatmentp", "eff", "distinct", "assay"
))
