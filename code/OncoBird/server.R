server <- shinyServer(function(input, output, server) {

  library(OncoBird)

  dataInputMut <- reactive({
    req(input$data_mutations)
    inFile <- input$data_mutations
    data_mutations <- as.data.frame(readr::read_csv(inFile$datapath))
  })

  dataInputClin <- reactive({
    req(input$data_clinical)
    inFile2 <- input$data_clinical
    data_clinical <- as.data.frame(readr::read_csv(inFile2$datapath))
  })

  surv <- reactive(input$select_oncology)
  treatment_col <- reactive(input$select_treatment)
  treatment_choice <- reactive(input$Select_treat)
  patid <- reactive(input$select_patientid)
  subtype <- reactive(input$Select_clin)

  ############## PREPARE DATA###################################################
  prepare_se <- function(datamut, dataclin, surv, treatment_col, patid) {
    prepare_data(
      data_mutations = datamut, # sample and treatment names can be supplied here!
      data_clinical = dataclin,
      med_impute_muts = T,
      vars = surv,
      treat = treatment_col,
      sample_column = patid
    )
  }

  prepareData <- eventReactive(input$go, {
    req(surv(), treatment_col(), patid(), dataInputMut(), dataInputClin())
    se <- prepare_se(datamut = dataInputMut(), dataclin = dataInputClin(), surv = surv(), treatment_col = treatment_col(), patid = patid())
  })

  ##################################################################################
  subtypeEnrichment <- reactive({
    req(prepareData())
    if (length(subtype()) > 1) {
      se <- cl_subtype_enrichment(
        se = prepareData(),
        col_label = subtype()[1],
        row_label = subtype()[2],
        digits = 3
      )
    } else {
      se <- NULL
    }
  })

  enrichmentGenomics <- reactive({
    if (!is.null(subtypeEnrichment())) {
      se <- cl_enrichment_genomics(subtypeEnrichment(),
        sample_column = "sample",
        min_mutants = input$min_mutants,
        subtype = subtype()
      )
    } else {
      se <- NULL
    }
  })

  # earliest point to compute this? after prepare data?
  mutualExclusive <- reactive({
    req(prepareData())
    req(input$prior_modules)
    withProgress(
      {
        if ("no_prior" %in% input$prior_modules) {
          prior <- NULL
        } else {
          prior <- input$prior_modules
        }
        se <- cl_mututal_exclusivity(se = prepareData(),
          min_variants = input$min_vars_exc,
          mutex_output_exists = !input$recalc_mutex,
          prior_modules = c(prior), # add mutex path
          mutex_path = "/",
          save = "./metadata"
        )
      },
      message = "Calculating Mutual Exclusive Modules"
    )
  })



  treatmentSpecificBiomarkers <- reactive({
    withProgress(
      {
        req(mutualExclusive())
        if (any(input$select_covariates == "no_covariates")) {
          se <- cl_treatment_specific_biomarkers(mutualExclusive(),
            include_covariates = NULL,
            min_samples = input$min_samples_tsb,
            min_redistribution = input$min_redis_tsb,
            treatment = treatment_choice(),
            subtypes = subtype(),
            readouts = surv()
          )
        } else {
          se <- cl_treatment_specific_biomarkers(
            se = mutualExclusive(),
            include_covariates = input$select_covariates,
            min_samples = input$min_samples_tsb,
            min_redistribution = input$min_redis_tsb,
            treatment = treatment_choice(), # select 2 values from treatment
            subtypes = subtype(),
            readouts = surv()
          )
        }
      },
      message = "Calculating Treatment Specific Biomarkers"
    )
  })

  treatmentSpecificBiomarkerSubtypes <- reactive({
    req(treatmentSpecificBiomarkers())
    se <- pl_treatment_specific_biomarkers_subtype(treatmentSpecificBiomarkers(),
      fdr_max = input$min_q_tsbs,
      colors = NULL,
      labels = NULL
    )
  })


  predictiveBiomarkers <- reactive({ # what does predictive Biomarkers depend on?
    withProgress(
      {
        req(treatmentSpecificBiomarkerSubtypes())
        req(input$select_covariates)
        if (any(input$select_covariates == "no_covariates")) {
          se <- cl_predictive_biomarkers(
            se = treatmentSpecificBiomarkerSubtypes(),
            include_covariates = NULL, # order matters here when more than 1 covariate! Intended? has to be reversed:  c("surgery","metastasis") --> c("metastasis","surgery")
            min_samples = input$min_samples_pb,
            min_redistribution = input$min_redis_pb,
            compare_treatments = T
          )
        } else if (length(input$select_covariates) > 1) {
          se <- cl_predictive_biomarkers(
            se = treatmentSpecificBiomarkerSubtypes(),
            include_covariates = c(input$select_covariates[1], c(input$select_covariates[2])), # order matters here when more than 1 covariate! Intended? has to be reversed:  c("surgery","metastasis") --> c("metastasis","surgery")
            min_samples = input$min_samples_pb,
            min_redistribution = input$min_redis_pb,
            compare_treatments = T
          )
        } else {
          se <- cl_predictive_biomarkers(
            se = treatmentSpecificBiomarkerSubtypes(),
            include_covariates = c(input$select_covariates[1]), # order matters here when more than 1 covariate! Intended? has to be reversed:  c("surgery","metastasis") --> c("metastasis","surgery")
            min_samples = input$min_samples_pb,
            min_redistribution = input$min_redis_pb,
            compare_treatments = T
          )
        }
      },
      message = "Calculating Predictive Biomarkers"
    )
  })


  predictiveBiomarkerSubtypes <- reactive({
    req(predictiveBiomarkers())
    se <- pl_predictive_biomarkers_subtypes(predictiveBiomarkers(),
      fdr_max = input$min_q_pbs,
      colors = NULL,
      labels = NULL
    )
  })

  predictiveComp <- reactive({ # what does predictive Biomarkers depend on?
    withProgress(
      {
        req(predictiveBiomarkerSubtypes(), input$min_q_predcomp)
        se <- cl_predictive_comparison(
          se = predictiveBiomarkerSubtypes(),
          subtypes = subtype(),
          readouts = surv(),
          covariates = F,
          fdr_i = input$min_q_predcomp
        )
      },
      message = "Calculating Predictive Comparison"
    )
  })

  plotExample <- reactive({
    req(predictiveComp(), input$select_example_mutation, input$select_example_readout, input$select_example_subtype_value, input$select_example_subtype)

    se <- pl_example(predictiveComp(),
      mutations = input$select_example_mutation,
      readout = input$select_example_readout,
      subtype_column = input$select_example_subtype,
      subtype = input$select_example_subtype_value,
      treatment = c(input$Select_treat)
    )
  })
  
  ############## Helpers ######################
  valuesTreatment <- reactive({
    req(prepareData())
    treatmentvalues <- unique(rowData(prepareData())$treatment) # should be the selected treatment column in the end
  })

  namesDataClin <- reactive({
    req(prepareData())
    if (is.null(prepareData())) {
      return()
    }
    col_choice <- colnames(rowData(prepareData()))
    final_choice <- setNames(col_choice, col_choice)
  })

  namesDataMut <- reactive({
    req(prepareData())
    if (is.null(prepareData())) {
      return()
    }
    col_choice <- colnames(assay(prepareData()))
    col_choice <- c("no_prior", unlist(lapply(col_choice, function(x) strsplit(x, "_")[[1]][1]))) # renaming mutations
    final_choice <- setNames(col_choice, col_choice)
  })

  namesSurvival <- reactive({
    req(dataInputClin())
    col_choice <- colnames(dataInputClin())
    vars <- c("OS", "PFS", "RFS", "DFS", "ORR")
    col_choice <- col_choice[col_choice %in% vars]
    final_choice <- setNames(col_choice, col_choice)
  })

  interestingOncology <- reactive({
    req(dataInputClin())
    vars <- c("OS", "PFS", "RFS", "DFS", "ORR")
    col_choice <- colnames(dataInputClin())
    final_choice <- setNames(vars, vars)
  })

  interestingTreatment <- reactive({
    req(dataInputClin())
    col_choice <- colnames(dataInputClin())
    final_choice <- setNames(col_choice, col_choice)
  })

  choices_tsb <- reactive({
    req(treatmentSpecificBiomarkers())
    plots <- pl_treatment_specific_biomarkers(treatmentSpecificBiomarkers())
    plot_names <- list()
    for (i in 1:length(plots)) {
      plot_names[[unique(plots[[i]]$labels$title)]] <- as.integer(i)
    }
    plot_names
  })

  choices_tsb_subtype <- reactive({
    req(treatmentSpecificBiomarkerSubtypes())
    plots <- metadata(treatmentSpecificBiomarkerSubtypes())$plot
    tsb_plots <- plots$treatment_specific_biomarkers_subtype
    plot_names <- list()
    for (i in 1:length(tsb_plots)) {
      plot_names[[names(tsb_plots)[i]]] <- names(tsb_plots)[i]
    }
    plot_names
  })

  choices_pb <- reactive({
    req(predictiveBiomarkers())
    plots <- pl_predictive_biomarkers(predictiveBiomarkers())
    plot_names <- list()
    for (i in 1:length(plots)) {
      plot_names[[unique(plots[[i]]$labels$title)]] <- as.integer(i)
    }
    plot_names
  })


  choices_pb_subtype <- reactive({
    req(predictiveBiomarkerSubtypes())
    plots <- metadata(predictiveBiomarkerSubtypes())$plot
    pb_plots <- plots$predictive_biomarkers_subtypes
    plot_names <- list()
    for (i in 1:length(pb_plots)) {
      plot_names[[names(pb_plots)[i]]] <- names(pb_plots)[i]
    }
    plot_names
  })

  choices_pc <- reactive({
    req(predictiveComp())
    plots <- pl_predictive_comparison(predictiveComp(),
      min_altered = 2,
      max_fdr = input$min_q_predcomp, 
      colors = NULL, labels = NULL
    )
    pc_choice <- names(plots)
    pc_choice
  })

  choices_example_mutation <- reactive({
    req(predictiveComp())
    mut_choice <- colnames(assay(predictiveComp())) # maybe order by significane
    mut_choice <- mut_choice[mut_choice != "sample"]
    names(mut_choice) <- mut_choice
    mut_choice
  })

  colnames_pc <- reactive({
    req(predictiveComp())
    colnames_pc <- colnames(rowData(predictiveComp())) # maybe order by significane
  })



  ############## Outputs ####################

  ############## User Selections ############

  ######## Prepare Data #####################
  output$select_oncology_column <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Select Clinical Endpoint of Interest")
    )
    selectizeInput("select_oncology", "Select Clinical Endpoint column", choices = interestingOncology(), multiple = T)
  })

  output$select_treatment_column <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Select Treatment of Interest")
    )
    selectizeInput("select_treatment", "Select treatment column", choices = interestingTreatment(), multiple = F)
  })

  output$select_patientid_column <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Select Patient ID column")
    )
    selectizeInput("select_patientid", "Select patientID column", choices = interestingTreatment(), multiple = F)
  })

  output$warnings <- renderText({
    req(prepareData())
    paste(metadata(prepareData())$shiny_message[["readouts"]],
      metadata(prepareData())$shiny_message[["imputed_mut"]],
      metadata(prepareData())$shiny_message[["non_zero_clins"]],
      metadata(prepareData())$shiny_message[["overlaps"]],
      metadata(prepareData())$shiny_message[["filtered"]],
      sep = "\n"
    )
  })


  ######## Subtype Enrichment #####################
  output$select_clinical <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Select up to 2 columns as tumour subtypes")
    )

    selectizeInput("Select_clin", "Choose multiple columns for as Tumour subtypes", choices = namesDataClin(), multiple = T, options = list(maxItems = 2))
  })

  ######## ME ##############################
  output$prior_modules <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Select mutations to group")
    )

    selectizeInput("Select_clin", "Choose 2 columns for grouping somatic mutations", choices = namesDataMut(), multiple = T, options = list(maxItems = 5))
  })

  ######## Biomarkers ##############################
  output$select_survival <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Choose survival strategy")
    )

    selectizeInput("Select_survival", "Choose survival strategy", choices = namesSurvival(), multiple = T, options = list(maxItems = 2))
  })

  output$select_covariates <- renderUI({
    if (is.null(input$data_clinical)) {
      return()
    }
    list(
      hr(),
      helpText("Select 0-2 Covariates columns for Biomarker Search")
    )
    all_options <- c(c("no_covariates"), namesDataClin())
    selectizeInput("select_covariates", "Choose up to 2 Covariates", choices = all_options, multiple = T, options = list(maxItems = 2))
  })

  output$select_treatment <- renderUI({
    req(prepareData())
    list(
      hr(),
      helpText("Select 2 Treatments")
    )

    selectizeInput("Select_treat", "Choose 2 Treatments", choices = valuesTreatment(), multiple = T, options = list(maxItems = 2))
  })

  ######## Example ##############################
  output$select_example_mutation <- renderUI({
    list(
      hr(),
      helpText("Choose a Mutation")
    )
    selectizeInput("select_example_mutation", "Choose a Mutation", choices = choices_example_mutation(), selected = 1, multiple = F)
  })

  output$select_example_readout <- renderUI({
    list(
      hr(),
      helpText("Select Survival Strategy")
    )
    selectizeInput("select_example_readout", "Select Survival Strategy", choices = surv(), selected = 1, multiple = F)
  })

  values_example <- reactive({
    req(input$select_example_subtype, predictiveComp())
    colname <- input$select_example_subtype
    data <- as.data.frame(rowData(predictiveComp()))
    values_example <- unique(data[[colname]])
  })

  output$select_example_subtype <- renderUI({
    list(
      hr(),
      helpText("Select Subtype Columnm")
    )
    add_ngs <- c(subtype(), "NGS.probe")
    selectizeInput("select_example_subtype", "Select Subtype Columnm", choices = add_ngs, selected = 1, multiple = F)
  })

  output$select_example_subtype_value <- renderUI({
    list(
      hr(),
      helpText("Select Subtype Columnm")
    )
    selectizeInput("select_example_subtype_value", "Select Subtype Columnm", choices = values_example(), selected = 1, multiple = F)
  })

  output$prior_modules <- renderUI({
    list(
      hr(),
      helpText("Select genes to group as prior")
    )
    selectizeInput("prior_modules", "Select Columns to group", choices = namesDataMut(), multiple = T)
  })

  ######## Plot Choices ###########################
  output$select_tsb_plot <- renderUI({
    tsb_plots <- choices_tsb()
    selectizeInput("Select_tsb_plot", "Select a Plot", choices = tsb_plots, selected = 1, multiple = F)
  })

  output$select_tsb_subtype_plot <- renderUI({
    tsb_subtype_plots <- choices_tsb_subtype()
    selectizeInput("select_tsb_subtype_plot", "Select a Plot", choices = tsb_subtype_plots, selected = 1, multiple = F)
  })

  output$select_pb_plot <- renderUI({
    selectizeInput("Select_pb_plot", "Select a Plot", choices = choices_pb(), selected = 1, multiple = F)
  })

  output$select_pb_subtype_plot <- renderUI({
    selectizeInput("select_pb_subtype_plot", "Select a Plot", choices = choices_pb_subtype(), selected = 1, multiple = F)
  })

  output$select_pc_plot <- renderUI({
    selectizeInput("select_pc_plot", "Select a Plot", choices = choices_pc(), selected = 1, multiple = F)
  })

 output$select_ml_barplot <- renderUI({
   selectizeInput("select_ml_barplot", "Select a Plot", choices = ml_barplotnames(), selected = 1, multiple = F)
  })
 
 output$select_permut_readout <- renderUI({
   selectizeInput("select_permut_readout", "Select Readout", choices = surv(), selected = 1, multiple = F)
 })

 output$select_boot_readout <- renderUI({
   selectizeInput("select_boot_readout", "Select Readout", choices = surv(), selected = 1, multiple = F)
 })
 
 output$select_permut_subtype <- renderUI({
   selectizeInput("select_permut_subtype", "Select Subtype", choices = subtype(), selected = 1, multiple = F)
 })
 
 output$select_boot_subtype <- renderUI({
   selectizeInput("select_boot_subtype", "Select Subtype", choices = subtype(), selected = 1, multiple = F)
 })



  ######## Data Table Output ###########################

  ########### Input Data Display ########################
  output$data_clin_table <- DT::renderDT({
    if (is.null(input$select_treatment)) {
      dataInputClin()
    } else {
      as.data.frame(rowData(prepareData()))
    }
  })


  ########### Summary Display & Download ################
  output$Summary_Table <- DT::renderDT({
    req(predictiveComp())
    se <- predictiveComp()
    
    if(any(input$select_covariates == "no_covariates")){
      cov <- c(1)
    }else{
      cov <- input$select_covariates
    }
    
    summary_table <- make_table_summary(
      cond = metadata(se)$conditions,
      preds = metadata(se)$predictive_biomarkers_forest,
      tsbiomarkers = metadata(se)$treatment_specific_biomarkers_forest,
      clin = rowData(se),
      mut = assay(se),
      mut_me = metadata(se)$me_modules,
      fdr_int_threshold = input$min_q_predcomp,
      fdr_biomarker_threshold = input$min_q_tsbs,
      round_digit = 2,
      sample_column = "sample",
      covariates = cov,
      treatment = metadata(se)$wf_meta$treatment
    )

    DT::datatable(summary_table,
      extensions = c("Buttons"),
      options = list(
        paging = TRUE,
        searching = TRUE,
        fixedColumns = TRUE,
        autoWidth = TRUE,
        ordering = TRUE,
        dom = "tB",
        buttons = list("copy", "print", list(
          extend = "collection",
          buttons = list(
            list(extend = "excel", filename = "OncoBird_Summary")
          ),
          text = "Download"
        ))
      ), class = "display"
    )
  })
 
 
 ########### Permutations Display & Download ################
 output$Permutation_Table <- DT::renderDT({

   req(predictiveComp(), input$select_permut_readout, input$select_permut_subtype)
   # subtpye always NGS.PROBE? 
   # Readout: Pick readout from first selections 
   
   subgroups_CIs  <- cl_permutations(se = predictiveComp(),
                                     treatment = treatment_choice(),
                                     subtype = input$select_permut_subtype, #"NGS.probe", input$select_permut_subtype
                                     readout = input$select_permut_readout,
                                     include_covariates = input$select_covariates,
                                     min_samples = input$min_samples_per,
                                     min_redistribution = input$min_redis_per,
                                     meta_load = FALSE,
                                     meta_path = "extdata",
                                     n_permutations = input$n_permut,
                                     fdr_int_threshold = input$fdr_int_permut,
                                     fdr_threshold = input$fdr_permut
   )

   DT::datatable(subgroups_CIs,
                 extensions = c("Buttons"),
                 options = list(
                   paging = TRUE,
                   searching = TRUE,
                   fixedColumns = TRUE,
                   autoWidth = TRUE,
                   ordering = TRUE,
                   dom = "tB",
                   buttons = list("copy", "print", list(
                     extend = "collection",
                     buttons = list(
                       list(extend = "excel", filename = "OncoBird_Permutations")
                     ),
                     text = "Download"
                   ))
                 ), class = "display"
   )
 })
 
 
 ########### Bootstrap Display & Download ################
 output$Bootstrap_Table <- DT::renderDT({
   # DEPENDS
   req(predictiveComp(), input$select_boot_readout, input$select_boot_subtype)
   
   subgroups <- cl_bootstrap(se = predictiveComp(),
                                 treatment = treatment_choice(),
                                 subtype = input$select_boot_subtype, 
                                 readout = input$select_boot_readout,
                                 include_covariates = input$select_covariates,
                                 min_samples = input$min_samples_boot,
                                 min_redistribution = input$min_redis_boot,
                                 meta_load = FALSE,
                                 meta_path = "metadata_paper",
                                 n_bootstraps = input$n_boot,
                                 fdr_int_threshold = input$fdr_boot,
                                 fdr_threshold = input$fdr_int_boot
   )
   
   
   DT::datatable(subgroups,
                 extensions = c("Buttons"),
                 options = list(
                   paging = TRUE,
                   searching = TRUE,
                   fixedColumns = TRUE,
                   autoWidth = TRUE,
                   ordering = TRUE,
                   dom = "tB",
                   buttons = list("copy", "print", list(
                     extend = "collection",
                     buttons = list(
                       list(extend = "excel", filename = "OncoBird_Bootstraps")
                     ),
                     text = "Download"
                   ))
                 ), class = "display"
   )
 })

  ######## Plot Outputs ###########################
  # 1
  output$Subtype_Enrichment <- renderPlot({
    validate(need(!is.null(subtypeEnrichment()), "Not possible with 1 subtype!"))
    pl_subtype_enrichment(subtypeEnrichment())
  })

  # 2
  output$Genomic_Enrichments <- renderPlot({
    validate(need(!is.null(enrichmentGenomics()), "Not possible with 1 subtype!"))
    pl_enrichment_genomics(enrichmentGenomics(),
      p_value = input$min_p_gen,
      fwe = input$min_q_gen
    )
  })
  # 3
  output$Mutual_Exclusive <- renderPlot({
    validate(need(!is.null(mutualExclusive()), "Choose Modules first"))
    pl_mutual_exclusivity(mutualExclusive())
  })

  # 4
  output$Oncoprint <- renderPlot({
    pl_oncoprint(mutualExclusive())
  })

  # 5
  output$Treatment_specific_biomarkers <- renderPlot({
    req(input$Select_tsb_plot, treatmentSpecificBiomarkers())
    pl_treatment_specific_biomarkers(treatmentSpecificBiomarkers())[as.integer(input$Select_tsb_plot)]
  })

  # 6
  output$Treatment_specific_biomarker_subtypes <- renderPlot({
    req(input$select_tsb_subtype_plot, treatmentSpecificBiomarkerSubtypes())
    plots <- metadata(treatmentSpecificBiomarkerSubtypes())$plot
    tsb_plots <- plots$treatment_specific_biomarkers_subtype
    tsb_plots <- tsb_plots[as.character(input$select_tsb_subtype_plot)]
    tsb_plots
  })

  # 7
  output$Predictive_Biomarkers <- renderPlot({
    req(input$Select_pb_plot)
    pl_predictive_biomarkers(predictiveBiomarkers())[as.integer(input$Select_pb_plot)]
  })

  # 8
  output$Predictive_Biomarker_Subtypes <- renderPlot({
    req(input$select_pb_subtype_plot, predictiveBiomarkerSubtypes())
    plots <- metadata(predictiveBiomarkerSubtypes())$plot
    pbs_plots <- plots$predictive_biomarkers_subtypes
    pbs_plots <- pbs_plots[as.character(input$select_pb_subtype_plot)]
    pbs_plots
  })

  # 9
  output$Predictive_Comp <- renderPlot({
    req(predictiveComp(), input$min_q_predcomp, input$select_pc_plot)
    plots <- pl_predictive_comparison(predictiveComp(),
      min_altered = 2,
      max_fdr = input$min_q_predcomp, colors = NULL, labels = NULL
    )

    gridExtra::grid.arrange(grobs = list(
      plots[[input$select_pc_plot]]$p_extraction[[1]] + theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm")) + theme(legend.position = "none") + ggtitle("mutant"),
      plots[[input$select_pc_plot]]$p_extraction[[2]] + theme(plot.margin = unit(c(1, 1, 3, 1.2), "cm")) + theme(legend.position = "right") + ggtitle("wild type")
    ), layout_matrix = rbind(c(1, 2)), top = input$select_pc_plot)
  })

  # 10
  output$Plot_Example <- renderPlot({
    req(plotExample())
    plotExample()
  })
})
