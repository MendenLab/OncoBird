#' Plot Subtype Mutations Enrichment
#'
#' With a binary event matrix investige enrichment of somatic mutations
#' in subtypes
#'
#' @param enr_list list of enrichment tests for different subtypes
#' from 'enrichment_genomics'
#' @param p_value_threshold p-value threshold for plotting
#' @param fdr_value_threshold fdr threshold for plotting
#' @param clin clinical dataframe
#' @param bem mutational dataframe
#' @param sort_by what level in subtypes to sort in
#' @param plot_only_one_if_two_subgroups binary if one plot for 2 subgroups
#' @param labels_list list of labels
#' @param colors_list list of colors
#'
#' @return ggplot object
#'
#' @noRd
plot_subtype_mutations_enrichment <- function(enr_list,
                                              p_value_threshold = 1,
                                              fdr_value_threshold = 0.95,
                                              clin,
                                              bem,
                                              sort_by = "2",
                                              plot_only_one_if_two_subgroups = TRUE,
                                              labels_list = NULL,
                                              colors_list = NULL) {
  bem_cox_fig4 <- bem
  if (length(enr_list) > 2) {
    stop("Currently, more than two subtypes are not supported for plotting.")
  }

  subgroups <- lapply(enr_list, function(x) as.character(unique(x$alteration)))

  subgroups <- lapply(
    subgroups,
    function(x) {
      if (plot_only_one_if_two_subgroups & length(x) == 2) {
        sort(x)[1]
      } else {
        x
      }
    }
  )

  tab <- data.table::rbindlist(enr_list)
  tab$alteration <- as.character(tab$alteration)
  tab <- tab %>%
    group_by(status) %>%
    filter(any(abs(value) <= p_value_threshold) &
      any(abs(fdr) <= fdr_value_threshold))

  if (nrow(tab) == 0) {
    stop("No significant enrichments found,
         please try a less stringent cutoff.")
  }

  sav_cox <- 2 - clin[, colnames(clin) %in% unique(tab$status), drop = FALSE] %>%
    mutate_all(factor) %>%
    mutate_all(as.numeric)
  if (ncol(sav_cox) > 0) {
    alts <- lapply(
      seq_len(unique(table(tab$alteration))),
      function(x) {
        cbind(
          bem_cox_fig4, sav_cox
        )[, as.character(unique(tab$status)[x])]
      }
    )
  } else {
    alts <- lapply(
      seq_len(unique(table(tab$alteration))),
      function(x) bem_cox_fig4[, as.character(unique(tab$status)[x])]
    )
  }
  # only add text to one, otherwise multiple text appear
  tab$text[tab$alteration == (sort(names(table(tab$alteration)))[1])] <-
    unlist(lapply(alts, function(x) {
      paste0(round(length(which(x == 1)) /
        length(which(x %in% c(0, 1))) * 100), "%")
    }))

  tab$status <- unlist(lapply(
    tab$status,
    function(x) gsub("_SV", " mutation", x)
  ))
  tab$status <- unlist(lapply(
    tab$status,
    function(x) gsub("_AMP", " amplification", x)
  ))
  tab$status <- unlist(lapply(
    tab$status,
    function(x) gsub("_LOSS", " loss", x)
  ))

  tab <- tab %>% # tell order
    group_by(status) %>%
    mutate(Order = min(sign(value[alteration == sort_by]) *
      log(abs(value[alteration == sort_by]))))

  enr_site_fig4 <- ggplot(
    data = tab[tab$alteration %in% subgroups[[1]], ],
    aes(
      x = reorder(`status`, Order), fill = `alteration`,
      y = sign(value) * -log10(abs(value))
    )
  ) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    labs(
      y = "signed log10(p-value)",
      x = "",
      fill = "Enrichment in\nrespective\nsubpopulation",
      tag = ""
    ) +
    theme_minimal() +
    coord_flip() +
    geom_hline(
      yintercept = log10(0.05), linetype = "dotted",
      color = "grey40", size = 0.8
    ) +
    geom_hline(
      yintercept = -log10(0.05), linetype = "dotted",
      color = "grey40", size = 0.8
    ) +
    geom_text(aes(label = text, x = reorder(status, abundance))) +
    theme(axis.text.y = element_text(size = 12))

  if (!any(is.null(labels_list)) & !any(is.null(colors_list))) {
    enr_site_fig4 <- enr_site_fig4 +
      scale_fill_manual(labels = labels_list[[1]], values = colors_list[[1]])
  }

  if (length(enr_list) == 2) { # if two enrichment test supplied
    enr_cms_fig4 <- ggplot(
      data = tab[tab$alteration %in% subgroups[[2]], ],
      aes(
        x = reorder(`status`, Order), fill = `alteration`,
        y = sign(value) * -log10(abs(value))
      )
    ) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      labs(
        y = "signed log10(p-value)",
        x = "",
        fill = "Enrichment in\nrespective\nsubpopulation",
        tag = ""
      ) +
      theme_minimal() +
      coord_flip() +
      geom_hline(
        yintercept = log10(0.05), linetype = "dotted",
        color = "grey40", size = 0.8
      ) +
      geom_hline(
        yintercept = -log10(0.05), linetype = "dotted",
        color = "grey40", size = 0.8
      ) +
      geom_text(aes(label = text, x = reorder(status, abundance)), y = 14) +
      theme(
        axis.text.y = element_blank()
      ) +
      ylim(c(-8, 15))

    if (!is.null(labels_list) & !is.null(colors_list)) {
      enr_cms_fig4 <- enr_cms_fig4 +
        scale_fill_manual(labels = labels_list[[2]], values = colors_list[[2]])
    }

    plt <- gridExtra::grid.arrange(
      grobs = list(
        enr_site_fig4 + theme(legend.position = "top"),
        enr_cms_fig4 + theme(legend.position = "top")
      ),
      layout_matrix = rbind(c(1, 2), c(1, 2), c(1, 2))
    )
  } else {
    plt <- gridExtra::grid.arrange(
      grobs = list(
        enr_site_fig4 + theme(legend.position = "top")
      ),
      layout_matrix = rbind(c(1), c(1), c(1))
    )
  }

  return(plt)
} #################################################



#' Plot Subtype Enrichments
#'
#'
#'
#' @param subtype_enrichment_df subtype enrichtment df from se 
#' @param col_label name of column 
#' @param row_label name of row
#' @param limits
#'
#' @return None
#' @details Plot the result of enrichment() e.g. right vs left sided enrichment
#' in subtypes
#'
#' @noRd
plot_subtype_enrichment <- function(subtype_enrichment_df,
                                    col_label,
                                    row_label,
                                    limits) {
  ggplot(
    subtype_enrichment_df,
    aes(
      x = .data[[col_label]], y = .data[[row_label]],
      fill = sign(value) * -log10(abs(value))
    )
  ) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu", limits = limits) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    labs(fill = "signed\n-log10(p-value)") +
    geom_text(aes(label = label, color = Enrichment)) +
    theme_minimal()
}



#' Plot Treatment Specific Biomarkers
#'
#'
#' @param df object of treatment specific biomarkers
#' @param color_highlight_con fdr cutoff
#' @param text_annotation what is the text column of df
#' @param max_text_annotation how many top significant text points
#' are highlighted
#' @param y.axis  which column gives y.axis
#' @param x.axis the column of effect size to plot
#' @param color color column
#' @param size size of points
#' @param fdr_label_x fdr label
#' @param text_size text size
#' @param nudge_pos nudge positive 
#' @param nudge_neg nudge negative
#' @param nudge_y nudge y
#'
#' @return list of n(row) plots
#'
#' @noRd
plot_treatment_specific_biomarkers <- function(df,
                                               color_highlight_con,
                                               text_annotation,
                                               max_text_annotation = 20,
                                               y.axis = "P.Value",
                                               x.axis = "effectsize",
                                               color = "strat",
                                               size = "confidence",
                                               fdr_label_x = NA,
                                               text_size = 5,
                                               nudge_pos = 2.7,
                                               nudge_neg = -2.7,
                                               nudge_y = 0) {
  color_highlight_con <- df$fdr < color_highlight_con

  max_p <- 1.3 * -log10(min(abs(df[, y.axis])))
  max_eff <- if (is.na(fdr_label_x)) {
    max(abs(df[, x.axis]))
  } else {
    max_eff <- fdr_label_x
  }
  if (length(which(color_highlight_con)) == 0) {
    howmany <- 1
  } else {
    howmany <- length(which(color_highlight_con))
  }
  gene_annotation_logical <-
    ((df[, y.axis]) <= sort(df[, y.axis],
      decreasing = FALSE
    )[min(max_text_annotation, howmany)])
  df$gene_annotation <- ifelse(gene_annotation_logical &
    color_highlight_con, df[, text_annotation], "")
  if (length(which(sort(na.omit(df[, "fdr"])) < 0.05)) > 0) {
    h0 <- sort(na.omit(df[, y.axis]))[max(which(sort(na.omit(
      df[, "fdr"]
    )) < 0.05))]
  } else {
    h0 <- 999999
  }
  if (length(which(sort(na.omit(df[, "fdr"])) < 0.1)) > 0) {
    h1 <- sort(na.omit(df[, y.axis]))[max(which(sort(na.omit(
      df[, "fdr"]
    )) < 0.1))]
  } else {
    h1 <- 999999
  }
  if (length(which(sort(na.omit(df[, "fdr"])) < 0.2)) > 0) {
    h2 <- sort(na.omit(df[, y.axis]))[max(which(sort(na.omit(
      df[, "fdr"]
    )) < 0.2))]
  } else {
    h2 <- 999999
  }
  if (length(which(sort(na.omit(df[, "fdr"])) < 0.3)) > 0) {
    h3 <- sort(na.omit(df[, y.axis]))[max(which(sort(na.omit(
      df[, "fdr"]
    )) < 0.3))]
  } else {
    h3 <- 999999
  }

  nudge_pos <- -(subset(df, color_highlight_con & (df[, x.axis] > 0))[, x.axis]) +
    nudge_pos
  if (length(nudge_pos) == 0) {
    nudge_pos <- 0
  }
  nudge_neg <- -(subset(df, color_highlight_con & (df[, x.axis] < 0))[, x.axis]) +
    nudge_neg
  if (length(nudge_neg) == 0) {
    nudge_neg <- 0
  }
  p2 <- ggplot() +
    geom_point(
      data = subset(df, !color_highlight_con),
      aes(
        x = !!(sym(x.axis)), y = -log10(!!sym(y.axis)),
        size = !!sym(size)
      ), color = "grey40", alpha = 0.5
    ) +
    geom_point(data = subset(df, color_highlight_con), aes(
      x = !!(sym(x.axis)),
      y = -log10(!!sym(y.axis)), color = factor(!!sym(color)),
      size = !!sym(size)
    ), alpha = 0.5) +
    labs(
      title = paste("", df$survival, " ", df$treatment, sep = ""),
      colour = "Stratification", size = "Confidence score"
    ) +
    xlab("Effect Size") +
    ylab(if (y.axis == "P.Value") {
      "-log10(p)"
    } else {
      "-log10(q)"
    }) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(0, max_p)) +
    xlim(c(-max_eff, max_eff)) +
    theme_minimal() +
    labs(color = color, size = size) +
    ggrepel::geom_text_repel(
      data = subset(df, color_highlight_con &
        (df[, x.axis] > 0)), aes(
        label = gene_annotation,
        x = !!(sym(x.axis)),
        y = -log10(!!sym(y.axis))
      ),
      force = 2,
      size = text_size,
      nudge_x = nudge_pos,
      direction = "y",
      hjust = 0,
      nudge_y = nudge_y,
      segment.size = 0.2
    ) +
    ggrepel::geom_text_repel(
      data = subset(df, color_highlight_con &
        (df[, x.axis] < 0)), aes(
        label = gene_annotation,
        x = !!(sym(x.axis)),
        y = -log10(!!sym(y.axis))
      ),
      force = 2,
      size = text_size,
      nudge_x = nudge_neg,
      direction = "y",
      hjust = 1,
      nudge_y = nudge_y,
      segment.size = 0.2
    )

  if (TRUE) {
    p2 <- p2 +
      geom_hline(
        yintercept = -log10(h0), linetype = "dashed",
        color = "grey66", size = 0.3
      ) +
      geom_hline(
        yintercept = -log10(h1), linetype = "dashed",
        color = "grey66", size = 0.3
      ) +
      geom_hline(
        yintercept = -log10(h2), linetype = "dashed",
        color = "grey66", size = 0.3
      ) +
      geom_hline(
        yintercept = -log10(h3), linetype = "dashed",
        color = "grey66", size = 0.3
      ) +
      annotate("text", x = -max_eff, y = if (-log10(h0) > -log10(0.05)) {
        -log10(h0)
      } else {
        99999
      }, label = "FDR=5%", hjust = 0, vjust = 0) +
      annotate("text", x = -max_eff, y = if (-log10(h1) > -log10(0.05)) {
        -log10(h1)
      } else {
        99999
      }, label = "FDR=10%", hjust = 0, vjust = 0) +
      annotate("text", x = -max_eff, y = if (-log10(h2) > -log10(0.05)) {
        -log10(h2)
      } else {
        99999
      }, label = "FDR=20%", hjust = 0, vjust = 0) +
      annotate("text", x = -max_eff, y = if (-log10(h3) > -log10(0.05)) {
        -log10(h3)
      } else {
        99999
      }, label = "FDR=30%", hjust = 0, vjust = 0)
  }

  return(p2)
}
######################################################



#' Get Annotation
#'
#' helper function
#'
#' @param annotation list of names for mutations for changing the naming
#'
#' @return list of updated names
#'
#' @noRd
get_annotation <- function(annotation) {
  tmp <- unlist(lapply(annotation, function(x) {
    gsub("_SV", "mut", x)
  }))
  tmp <- unlist(lapply(tmp, function(x) {
    gsub("_AMP", "amp", x)
  }))
  tmp <- unlist(lapply(tmp, function(x) {
    gsub("_LOSS", "loss", x)
  }))
  tmp <- unlist(lapply(tmp, function(x) {
    gsub("_or_", "/", x)
  }))
  tmp_res <- tmp
  if (is.factor(annotation)) {
    levels <- levels(annotation)
    tmp <- unlist(lapply(levels, function(x) {
      gsub("_SV", "mut", x)
    }))
    tmp <- unlist(lapply(tmp, function(x) {
      gsub("_AMP", "amp", x)
    }))
    tmp <- unlist(lapply(tmp, function(x) {
      gsub("_LOSS", "loss", x)
    }))
    tmp <- unlist(lapply(tmp, function(x) {
      gsub("_or_", "/", x)
    }))
    tmp_res <- factor(tmp_res, levels = tmp)
  }
  return(tmp_res)
}


#' Fix Annotation
#'
#' helper functions
#'
#' @param annotation list of names for mutations for changing the naming
#'
#' @return list of updated names
#'
#' @noRd
fix_annotation <- function(annotation) {
  tmp_anno <- lapply(annotation, function(x) strsplit(x, "/")[[1]])
  tmp <- lapply(tmp_anno, function(x) {
    gsub("amp", "", x)
  })
  tmp <- lapply(tmp, function(x) {
    gsub("loss", "", x)
  })
  tmp <- lapply(tmp, function(x) {
    gsub("mut", "", x)
  })
  tmp <- lapply(tmp, function(x) paste0(unlist(x[duplicated(x)]), "alt"))
  tmp_anno <- lapply(
    seq_len(length(tmp_anno)),
    function(x) {
      if (any(tmp[[x]] != "alt")) {
        c(tmp[[x]], tmp_anno[[x]])
      } else {
        tmp_anno[[x]]
      }
    }
  )
  keep <- lapply(
    tmp_anno,
    function(x) {
      !duplicated(gsub("alt", "", gsub(
        "loss", "",
        gsub("amp", "", gsub("mut", "", x))
      )))
    }
  )
  tmp_anno <- unlist(lapply(
    seq_len(length(tmp_anno)),
    function(x) paste((tmp_anno[[x]])[keep[[x]]], collapse = "/")
  ))
  return(tmp_anno)
}



#' Plot Forest Treatment Specific Biomarkers
#'
#' make forest plots for treatment- and subtype-specific biomarkers
#'
#' @param forest dataframe for association tests
#' @param identifier "identifier"
#' @param fdr_threshold fdr threshold for plotting significant markers
#' @param fdr_column column name in the forest dataframe for fdr threshold
#' plotting
#' @param colors list of colors
#' @param labels list of labels
#' 
#'
#' @return ggplot object
#'
#' @noRd
plot_forest_treatment_specific_biomarkers <- function(forest,
                                                      identifier = "identifier",
                                                      fdr_threshold = 0.3,
                                                      fdr_column = "fdr",
                                                      colors = NULL,
                                                      labels = NULL) {
  fig4new <- forest
  fig4new$FDR <- fig4new[, fdr_column]
  fig4new$identifier <-
    factor(fig4new[, identifier], levels = unique(fig4new[, identifier]))
  fig4new$gene <-
    fix_annotation(get_annotation(annotation = fig4new$gene))

  Levels <-
    (fig4new %>% filter(strat == "available") %>%
      arrange(identifier, hazard))$gene %>%
    unique()
  fig4new$gene <- factor(fig4new$gene, levels = Levels)
  label <- unique(fig4new$identifier)
  fig4new <- fig4new %>%
    group_by(gene) %>%
    filter(any(FDR <= fdr_threshold))
  if (nrow(fig4new) == 0) {
    warning("Increase FDR threshold for ", unique(label))
    p_extraction <- NULL
  } else {
    fig4new$strat_variable <-
      factor(fig4new$strat, levels = fig4new$strat %>% unique())
    fig4new$strat_variable <- as.character(fig4new$strat_variable)
    fig4new$strat_variable[fig4new$FDR > fdr_threshold] <- "not significant"

    if (any(abs(fig4new$upper) == Inf | abs(fig4new$lower) == Inf)) {
      w <- !(abs(fig4new$upper) == Inf | abs(fig4new$lower) == Inf)
      fig4new <- fig4new[w, ]
      warning(
        "Removed ", length(which(w)),
        " entry with non-numerical intervals."
      )
    }

    p_extraction <- ggplot() +
      geom_pointrange(
        data = subset(fig4new),
        aes(
          x = gene, y = hazard, ymin = lower, ymax = upper,
          group = strat,
          color = factor(strat_variable)
        ),
        position = position_dodge(width = 0.5),
        alpha = 0.6
      ) +
      facet_wrap(~identifier, ncol = 1) +
      geom_hline(yintercept = 1, lty = 2) +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10)
      ) +
      labs(color = "Subgroup", y = "hazard ratio") +
      theme(
        plot.margin = unit(c(0.2, 0, 0.35, 0), "inch"),
        legend.box.spacing = unit(0.5, "inch")
      ) +
      scale_y_continuous(trans = "log10") +
      xlab("")

    if (!any(is.null(colors)) & !any(is.null(labels))) {
      p_extraction <-
        p_extraction <-
        p_extraction + scale_color_manual(values = colors, labels = labels)
    }
  }
  return(p_extraction)
}



######################################################
#' Plot Forest Comparison
#'
#' make treatment forest plot from predictive biomarkers
#'
#' @param comp list of comparisons from compare_treatments
#' @param min_altered minimum of cell lines in a cohort
#' @param fdr_threshold (deprecated)
#' @param fdr_biomarker_threshold fdr cutoff for treatment-specific biomarkers
#' @param treatment_biomarkers treatment-specific biomarkers for filtering of
#' interactions tests
#' @param colors list of colors
#' @param labels list of labels
#'
#' @return ggplot object for treatment comparisons
#'
#' @noRd
plot_forest_comparison <- function(comp,
                                   min_altered = 12,
                                   fdr_threshold = 1,
                                   fdr_biomarker_threshold = 1,
                                   treatment_biomarkers = NULL,
                                   colors = NULL,
                                   labels = NULL) {
  res_f <- comp

  if (all(unlist(lapply(res_f, is.null)))) {
    stop("All predictive biomarkers are insignificant..., increase threshold..")
  }
  res_f <- make_forest_comp(do.call(rbind, res_f))
  # insert filtering for treatment_specific biomarkers
  if (!is.null(treatment_biomarkers)) {
    which_biomarkers <-
      apply(treatment_biomarkers[treatment_biomarkers$fdr <
        fdr_biomarker_threshold, c("cfe", "strat")], 1, function(x) {
        paste(x, collapse = "---")
      })
    res_f$help <-
      apply(res_f[, c("gene", "strat_variable")], 1, function(x) {
        paste(x, collapse = "---")
      })
    res_f <-
      res_f[res_f$help %in% unlist(which_biomarkers), , drop = FALSE]
    res_f <- res_f %>% dplyr::select(-help)
  }
  ####################################################

  if (nrow(res_f) == 0) {
    stop("All predictive biomarkers are insignificant..., increase threshold..")
  }

  res_f$gene <- fix_annotation(get_annotation(annotation = res_f$gene))
  res_f$mutantorwildtype <-
    res_f$mutantorwildtype %>%
    as.character() %>%
    dplyr::recode(
      "0" = "wild type",
      "1" = "altered"
    )
  res_f$genemutantorwildtype <- paste0(res_f$gene, res_f$mutantorwildtype)
  Levels <- (res_f %>% arrange(hazard))$gene %>% unique()
  res_f$gene <- factor(res_f$gene, levels = Levels)
  res_f_all <- res_f
  res_f <- res_f[res_f$confidence >= min_altered, ]

  if (nrow(res_f) == 0) {
    stop("All predictive biomarkers are insignificant..., increase threshold..")
  }

  res_f <-
    res_f %>%
    group_by(gene, strat_variable) %>%
    filter(any(fdr <= fdr_threshold))
  res_f$strat_variable <-
    factor(res_f$strat_variable, levels = res_f$strat_variable %>% unique())
  res_f <-
    res_f[!duplicated(paste0(res_f$genemutantorwildtype, res_f$strat_variable)), ]
  res_f <-
    res_f[res_f$gene %in% names(which(!apply(as.data.frame.matrix(table(
      res_f$gene, res_f$mutantorwildtype
    )), 1, function(x) {
      any(x == 0)
    }))), ]
  res_f$p_os <- res_f$p
  res_f$eff_os <- res_f$eff

  if (nrow(res_f) == 0) {
    stop("All predictive biomarkers are insignificant..., increase threshold..")
  }

  which.marker <- res_f_all
  which.marker <- which.marker[which.marker$confidence >= min_altered, ]
  which.marker <- which.marker %>%
    group_by(gene, strat_variable) %>%
    filter(any(fdr <= fdr_threshold))

  if (nrow(which.marker) == 0) {
    stop("All predictive biomarkers are insignificant..., increase threshold..")
  }

  which.marker$strat_variable <- factor(which.marker$strat_variable,
    levels = which.marker$strat_variable %>% unique()
  )
  which.marker <- which.marker[which.marker$gene %in% names(which(!apply(
    as.data.frame.matrix(table(
      which.marker$gene, which.marker$mutantorwildtype
    )), 1,
    function(x) any(x == 0)
  ))), ]
  which.marker <- which.marker[which.marker$fdr <= fdr_threshold, ]

  if (nrow(which.marker) == 0) {
    stop("All predictive biomarkers are insignificant..., increase threshold..")
  }

  which.marker <- which.marker %>%
    mutate(biomarker = ifelse((mutantorwildtype == "altered") &
      (direction == "sensitive"),
    paste0(strsplit(strat, "_")[[1]][1], "-", "sensitivity"),
    paste0(strsplit(strat, "_")[[1]][1], "-", "resistance")
    ))
  which.biomarkers <- ggplot(which.marker, aes(x = gene, y = biomarker),
    fill = "grey"
  ) +
    geom_tile() +
    coord_flip() +
    theme_classic() +
    scale_fill_continuous(na.value = "grey22")

  res_f$group <- res_f$strat_variable
  res_f$group <- as.character(res_f$group)
  res_f$group[res_f$fdr > fdr_threshold] <- "not significant"
  res_f$group <- factor(res_f$group)
  ###


  p_extraction <- list()
  for (i in c("altered", "wild type")) {
    p_extraction[[i]] <- ggplot(
      data = subset(res_f, mutantorwildtype == i),
      aes(x = gene, y = hazard, ymin = lower, ymax = upper)
    ) +
      geom_pointrange(aes(
        color = group,
        group = strat_variable
      ),
      position = position_dodge(width = 0.5),
      alpha = 0.6
      ) +
      coord_flip(clip = "off") +
      geom_hline(yintercept = 1, lty = 2) +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = if (i == "altered") {
          element_text(angle = 45)
        } else {
          element_blank()
        }, axis.ticks.y = element_blank(),
        plot.title = element_text(size = 10)
      ) +
      labs(color = "Subgroup", y = "hazard ratio") +
      theme(
        plot.margin = unit(c(0.2, 0, 0.35, 0), "inch"),
        legend.box.spacing = unit(0.5, "inch")
      ) +
      scale_y_continuous(trans = "log10") +
      xlab("") +
      theme(legend.position = "none")

    if (!any(is.null(colors)) & !any(is.null(labels))) {
      p_extraction[[i]] <- p_extraction[[i]] <- p_extraction[[i]] +
        scale_color_manual(values = colors, labels = labels)
    }
  }

  return(list(
    p_extraction = p_extraction,
    which.biomarkers = which.biomarkers
  ))
}



#' Plot Example
#'
#' plot any example for the presented analysis
#'
#' @param bem mutational dataframe
#' @param clin clinical dataframe
#' @param sample_column column specifying individuals
#' @param subtype_column column name specifying subtype
#' @param subtype factor level of subtype
#' @param treatment_column treatment column name
#' @param survival readout (survival or binary)
#' @param alteration investigated mutations
#' @param covariates confounding factors
#'
#' @return list of ggplot objects
#'
#' @noRd
plot_example <- function( ######################################################
                         ######################################################
                         bem,
                         clin,
                         sample_column = "sample",
                         subtype_column,
                         subtype,
                         treatment_column = "treatment",
                         survival,
                         alteration,
                         covariates = c(1)) {
  intersecting_samples <- intersect(bem[, sample_column], clin[, sample_column])
  bem <- bem[bem[, sample_column] %in% intersecting_samples, ]
  clin <- clin[clin[, sample_column] %in% intersecting_samples, ]
  data <- cbind(bem, clin)
  tlevels <- rev(levels(factor(unlist(data[, treatment_column]))))
  if (length(tlevels) > 2) {
    stop("Currently, more than 2 treatment arms are not supported,
         so reduce treatment arms to two !")
  }
  alt <- expand.grid(
    status = c("altered", "wild type"),
    alt = alteration, treatment = tlevels
  )
  message(alt)
  alt$color <- c("darkblue", "blue", "darkred", "red")
  subtypeadj <- subtype # label in plot lab

  list_1 <- list()
  for (treatment in tlevels) {
    model <- plot_cox(
      response = data,
      features = data,
      feature_name = alteration,
      survival = survival,
      strat = subtype_column,
      treatment = treatment,
      covariates = covariates,
      N = 2,
      rna_asso = FALSE
    )[[subtype]]

    if (is.null(model$cox) == FALSE) {
      ps <- summary(model$cox)$coefficients[1, 5]
      message("Specific:")
      message(ps)
    } else {
      message("Conditions cannot be compared,
              most likely too little data points available")
    }
    alt1 <- alt
    alt1$color <- c("darkred", "red","darkblue", "blue")
    list_1[[treatment]] <- ggplot2::autoplot(model$fit) +
      labs(
        fill = paste(
          fix_annotation(get_annotation(alteration)), "in", subtypeadj
        ),
        color = paste(fix_annotation(get_annotation(alteration)), "in", subtypeadj)
      ) +
      scale_color_manual(
        labels = rev(alt1$status[alt$treatment == treatment]),
        values = rev(alt1$color[alt$treatment == treatment])
      ) +
      scale_fill_manual(
        labels = rev(alt1$status[alt$treatment == treatment]),
        values = rev(alt1$color[alt$treatment == treatment])
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "vertical"
      ) +
      xlab("Months since randomisation") +
      ylab(survival)
  }

  list_2 <- list()
  for (status in c("altered", "wild type")) {
    data_here <- data
    data_here$treatment_strat <- data_here$treatment
    data_here$treatment <- "yes"
    if (all(c("0", "1") %in% data_here[, alteration])) {
      if (status == "altered") {
        which <- "1"
      }
      if (status == "wild type") {
        which <- "0"
      }
    } else {
      if (status == "altered") {
        which <- factor(data_here[, alteration])[1]
      }
      if (status == "wild type") {
        which <- factor(data_here[, alteration])[2]
      }
    }
    data_here <- data_here[data_here[, alteration] == which, ]
    model <- plot_cox(
      response = data_here,
      features = data_here,
      feature_name = "treatment_strat",
      survival = survival,
      strat = subtype_column,
      treatment = "yes",
      covariates = covariates,
      N = 2,
      rna_asso = TRUE
    )[[subtype]]

    if (is.null(model$cox) == FALSE) {
      ps <- summary(model$cox)$coefficients[1, 5]

      message("Comparison:")
      message(ps)
    }


    list_2[[status]] <- ggplot2::autoplot(model$fit) +
      labs(
        fill = paste(paste(
          "", fix_annotation(get_annotation(alteration)), status
        ), "in", subtypeadj),
        color = paste(paste(
          "", fix_annotation(get_annotation(alteration)), status
        ), "in", subtypeadj)
      ) +
      scale_color_manual( # labels= alt$treatment[alt$status == status]
        values = (alt$color[alt$status == status])
      ) +
      scale_fill_manual( # labels= alt$treatment[alt$status == status]
        values = (alt$color[alt$status == status])
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "vertical"
      ) +
      xlab("Months since randomisation") +
      ylab(survival)
  }

  list_3 <- list()
  data_here <- data
  data_here$treatment_strat <- data_here$treatment
  data_here$treatment <- "yes"
  data_here <- data_here[!is.na(data_here[, alteration]), ]
  data_here$info <-
    factor(
      paste(
        data_here$treatment_strat,
        "and",
        fix_annotation(get_annotation(alteration)),
        data_here[, alteration] %>%
          dplyr::recode("0" = "wild type", "1" = "altered")
      )
    )

  model <- plot_cox(
    response = data_here,
    features = data_here,
    feature_name = "info",
    survival = survival,
    strat = subtype_column,
    treatment = "yes",
    covariates = covariates,
    N = 2,
    rna_asso = TRUE
  )[[subtype]]

  list_3 <- ggplot2::autoplot(model$fit) +
    labs(
      fill = paste("treatment and patients", "in", subtypeadj),
      color = paste("treatment and patients", "in", subtypeadj)
    ) +
    scale_color_manual(
      values = alt$color
    ) +
    scale_fill_manual(
      values = alt$color
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical"
    ) +
    xlab("Months since randomisation") +
    ylab(survival)

  return(list(specific = list_1, comparison = list_2, all = list_3))
}
######################################################


#' Plot Modules treatment Specific Biomarkers
#'
#' plot treatment-specific modules
#'
#' @param assos dataframe of treatment-specific biomarkers
#' @param strat subtype
#' @param survival readout
#' @param p_value p-value threshold
#' @param fdr fdr threshold
#' @param bem mutational dataframe
#'
#' @return list of ggplot objects
#'
#' @noRd
plot_modules_oncoprint_treatment_specific_biomarkers <- function(assos = treatment_specific_biomarkers,
                                                                 strat = "available",
                                                                 survival = "OS",
                                                                 p_value = 0.05,
                                                                 fdr = 0.3,
                                                                 bem = data_mutations_me
                                                                 ######
) {
  assos <- assos[assos$strat == strat & assos$survival == survival, ]
  indexing <- expand.grid(
    levels(unlist(assos["treatment"])),
    c("sensitive", "resistant")
  )
  fig4_plots <- list()
  for (J in seq_len(nrow(indexing))) {
    sor <- indexing[J, "Var2"]
    asso <- assos[assos$treatment == indexing[J, "Var1"], ]

    if (nrow(asso) == 0) {
      asso_cox_cet_interactions_anno <- asso # assign empty df that triggers next if
    } else {
      asso_cox_cet_interactions_anno <- cbind(asso, get_anno_combinations(asso))

      asso_cox_cet_interactions_anno <- asso_cox_cet_interactions_anno %>%
        mutate_at("gene", function(x) {
          factor(
            x,
            levels = (asso_cox_cet_interactions_anno[order(
              -log10(asso_cox_cet_interactions_anno$P.Value) *
                sign(asso_cox_cet_interactions_anno$effectsize) *
                if (sor == "sensitive") {
                  -1
                } else {
                  -1
                },
              decreasing = TRUE
            ), "gene"])
          )
        })
      asso_cox_cet_interactions_anno <- asso_cox_cet_interactions_anno[
        asso_cox_cet_interactions_anno$P.Value < p_value &
          !is.na(asso_cox_cet_interactions_anno$P.Value) &
          asso_cox_cet_interactions_anno$fdr < fdr,
      ]
    }


    if (nrow(asso_cox_cet_interactions_anno) == 0) {
      fig4_plots[[J]] <- list(NULL, NULL, NULL)
    } else {
      asso_cox_cet_interactions_anno <-
        get_forest_new(asso_cox_cet_interactions_anno)
      asso_cox_cet_interactions_anno <-
        asso_cox_cet_interactions_anno[order(
          asso_cox_cet_interactions_anno$P.Value,
          decreasing = FALSE
        ), ]

      # only sensitive or resistant
      asso_cox_cet_interactions_anno <- asso_cox_cet_interactions_anno[
        if (sor == "sensitive") {
          asso_cox_cet_interactions_anno$effectsize < 0
        } else {
          asso_cox_cet_interactions_anno$effectsize > 0
        },
      ]

      if (nrow(asso_cox_cet_interactions_anno) == 0) {
        fig4_plots[[J]] <- list(NULL, NULL, NULL)
      } else {
        if (survival == "ORR") {
          asso_cox_cet_interactions_anno$hazard <-
            exp(asso_cox_cet_interactions_anno$hazard)
          asso_cox_cet_interactions_anno$lower <-
            exp(asso_cox_cet_interactions_anno$lower)
          asso_cox_cet_interactions_anno$upper <-
            exp(asso_cox_cet_interactions_anno$upper)
        }

        bar1 <- ggplot(
          data = asso_cox_cet_interactions_anno,
          aes(
            x = gene,
            y = hazard,
            ymin = lower,
            ymax = upper,
            color = confidence
          )
        ) +
          geom_pointrange() +
          coord_cartesian(clip = "off") +
          geom_hline(yintercept = 1, lty = 2) +
          theme_minimal() +
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 10)
          ) +
          labs(color = "log10(number of\naltered patients)", y = "hazard ratio") +
          theme(
            plot.margin = unit(c(0.2, 0, 0.35, 0), "inch"),
            legend.box.spacing = unit(0.9, "inch")
          ) +
          scale_y_continuous(trans = "log10") +
          geom_text(
            aes(
              label = paste0("p=", as.character(rstatix::p_round(P.Value))),
              x = gene,
              y = max(lower) + 0.8
            ),
            hjust = 0, angle = 45
          )


        hm <-
          colnames(asso_cox_cet_interactions_anno)[!colnames(asso_cox_cet_interactions_anno) %in%
            c(
              "hazard",
              "upper",
              "lower",
              "effectsize",
              "P.Value",
              "cfe",
              "strat",
              "npos",
              "nneg",
              "model",
              "medians",
              "fdr",
              "fdr_old",
              "gene",
              "confidence",
              "survival",
              "treatment",
              "genelabels",
              "annotation",
              "color",
              "identifier"
            )]
        data <- asso_cox_cet_interactions_anno[, c(hm, "gene")] %>%
          tidyr::pivot_longer(cols = hm)
        Order <- rev(data$name[data$value == 1 & !(data$name %in% c("NA"))] %>%
          unique())
        data <- data %>%
          dplyr::group_by(name) %>%
          dplyr::filter(
            sum(value, na.rm = TRUE) > 0 & name != "NA" & !is.na(gene)
          )
        data$gene <- factor(data$gene)
        levels <- levels(data$gene)

        is.single <- lapply(
          as.character(data$gene),
          function(x) {
            length(gsub("\\)", "", gsub(
              "\\(", "",
              strsplit(x, c("_and_", "_or_"))
            ))) == 1 &
              !(x %in% c("RAS_SV_or_BRAF_SV", "RAS_SV"))
          }
        ) %>% unlist()
        data$value[is.single & (data$value == 1)] <- 2
        data <- data[!is.na(data$name), ]

        data$gene <- factor(data$gene, levels = levels)
        data$name <- factor(data$name, levels = Order)
        bar2 <- ggplot(
          data = data,
          aes(
            x = get_annotation(name),
            y = gene,
            fill = factor(value)
          )
        ) +
          geom_tile() +
          scale_fill_manual(values = c("grey80", "grey50", "black")) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1.1),
            plot.title = element_text(size = 10),
            legend.position = "none"
          ) +
          xlab("") +
          ylab("") +
          theme(
            plot.margin = unit(c(0.2, 0, 0.1, 0.6), "inch"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 9),
            panel.background = element_rect(color = "grey80", fill = "grey80"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          ) +
          coord_flip()
        message(rev(levels(data$name)))
        arg1 <- lapply(rev(levels(data$name)), function(x) bem[, x])
        bem_here <- bem[do.call(order, arg1), rev(levels(data$name)), drop = FALSE]
        rownamess <- row.names(bem_here)

        ### oncoprint snip
        bem_oncoprint_type <-
          lapply(colnames(bem_here), function(x) {
            strsplit(x, "_")[[1]][2]
          }) %>% unlist()

        bem_oncoprint_type[bem_oncoprint_type == "any"] <- "AMP"
        bem_oncoprint_gene <-
          lapply(colnames(bem_here), function(x) {
            strsplit(x, "_")[[1]][1]
          }) %>% unlist()

        for (x in seq_len(length(bem_oncoprint_type))) {
          bem_here[, x][bem_here[, x] == 1] <- bem_oncoprint_type[x]
        }
        anno <- bem_here %>% t()
        anno[anno == "SV"] <- "MUT"
        anno[anno == "LOSS"] <- "LOSS"
        anno[anno == "AMP"] <- "AMP"
        anno[anno == "0"] <- ""
        anno[!(anno %in% c("MUT", "LOSS", "AMP", ""))] <- "FEAT"
        anno <- anno %>% as.data.frame()
        anno$gene <- colnames(bem_here)

        anno_long <-
          tidyr::pivot_longer(
            data = anno,
            colnames(anno)[-length(colnames(anno))]
          )
        anno_long <- anno_long %>%
          group_by(name, gene) %>%
          summarise(alterations = paste(unique(value), collapse = ";"))
        anno <-
          tidyr::pivot_wider(
            data = anno_long,
            names_from = gene,
            values_from = alterations,
            id_cols = name
          ) %>% tibble::column_to_rownames("name")

        bar3 <- ComplexHeatmap::oncoPrint(
          anno %>% t(),
          alter_fun = alter_fun,
          get_type = function(x) {
            strsplit(x, ";")[[1]]
          },
          row_order = match(rev(Order), colnames(anno)),
          remove_empty_columns = TRUE,
          show_row_names = FALSE,
          col = c(
            "MUT" = "grey30",
            "AMP" = "darkred",
            "LOSS" = "darkblue",
            "NA" = "white",
            "FEAT" = "yellow"
          ),
          heatmap_legend_param = list(
            title = "Alternations",
            at = c("AMP", "LOSS", "MUT", "NA", "FEAT"),
            labels = c("Amplification", "Loss", "Mutation", "NA", "other")
          )
        )
        ##################

        index <- paste(as.character(unlist(indexing[J, ])), collapse = "/")
        fig4_plots[[index]] <- list(bar3, bar2, bar1)
      }
    }
  }

  return(fig4_plots)
}


##' Plot mutual exclusivity
#'
#' Plot mutex scores of mutually exclusive modules
#'
#' @param ngs_exclusive dataframe from mutex algorithm after preparation
#'
#' @return None
#'
#' @noRd
plot_mutual_exclusivity <- function(ngs_exclusive) {
  ngs_exclusive_plot <- ngs_exclusive
  ngs_exclusive_plot$name <- gsub(";", "/", ngs_exclusive_plot$name)
  ngs_exclusive_plot$identifier <-
    lapply(seq_len(nrow(ngs_exclusive)), function(x) {
      tmp <-
        as.character(unlist(ngs_exclusive[x, 2:(ncol(ngs_exclusive_plot) - 1)]))
      tmp <- tmp[tmp != ""]
      return(sort(unlist(tmp)))
    })
  ngs_exclusive_plot$Score <-
    as.numeric(gsub(",", ".", gsub(
      "E", "e", as.character(ngs_exclusive_plot[, 1])
    )))
  ngs_exclusive_plot <-
    ngs_exclusive_plot[!duplicated(ngs_exclusive_plot$identifier), ]
  ggplot(ngs_exclusive_plot) +
    geom_bar(aes(
      x = reorder(name, -log10(Score)),
      y = -log10(Score)
    ), stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    # coord_cartesian(y = c(0,1.1))+
    xlab("Mutually exclusive module") +
    ylab("Exclusivity score") +
    theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm"))
}
