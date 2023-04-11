ui <- shinyUI(
  fluidPage(
    titlePanel("OncoBird"),
    tabsetPanel(
      tabPanel(
        "Upload File",
        titlePanel("Uploading Files"),
        sidebarLayout(
          sidebarPanel(
            fileInput("data_mutations", 'Choose "data_mutations"',
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"
              )
            ),
            tags$br(),
            fileInput("data_clinical", 'Choose "data_clinical"',
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"
              )
            ),
            uiOutput("select_oncology_column"),
            uiOutput("select_treatment_column"),
            uiOutput("select_patientid_column"),
            actionButton(inputId = "go", label = "Prepare data"),
            uiOutput("select_clinical"),
            uiOutput("select_treatment"),
            uiOutput("select_covariates"),
            submitButton("Submit"),
            tags$br(),
            checkboxInput("recalc_mutex", "Recalculate Mutex", FALSE)
          ),
          mainPanel(
            DT::DTOutput(outputId = "data_clin_table"),
            verbatimTextOutput("warnings")
          )
        ),
      ),
      ## SUBTYPE ENRICHMENT
      tabPanel(
        "Subtype Enrichment",
        headerPanel("Subtype Enrichment"),
        mainPanel(
          plotOutput(outputId = "Subtype_Enrichment")
        )
      ),
      ## GENOMIC ENRICHMENT
      tabPanel(
        "Genomic_Enrichments",
        headerPanel("Genomic Enrichments"),
        sidebarPanel(
          sliderInput(inputId = "min_mutants", label = "Select Minimum number of Mutants", min = 0, max = 100, value = 10, step = 1), # max dynamic?

          tags$br(),
          sliderInput(inputId = "min_p_gen", label = "Select P-Value threshhold", min = 0, max = 1, value = 0.1, step = 0.01),
          sliderInput(inputId = "min_q_gen", label = "Select FDR threshhold", min = 0, max = 1, value = 1, step = 0.01),
          submitButton("Submit")
        ),
        mainPanel(
          plotOutput(outputId = "Genomic_Enrichments")
        )
      ),
      ## MUTUAL EXCLUSIVE
      tabPanel(
        "Mutual Exclusivity",
        headerPanel("Mutual Exclusivity"),
        sidebarPanel(
          sliderInput(inputId = "min_vars_exc", label = "Select Minimum number of Variants", min = 0, max = 100, value = 10, step = 1), # maybe make this max valu flexible
          uiOutput("prior_modules"),
          submitButton("Submit")
        ),
        mainPanel(
          plotOutput(outputId = "Mutual_Exclusive")
        )
      ),
      ## ONCOPRINT
      tabPanel(
        "Oncoprint",
        headerPanel("Oncoprint"),
        mainPanel(
          plotOutput(outputId = "Oncoprint")
        )
      ),

      ## TREATMENT SPECIFIC BIOMARKERS
      tabPanel(
        "Treatment Specific Biomarkers",
        headerPanel("Treatment Specific Biomarkers"),
        sidebarPanel(
          sliderInput(inputId = "min_samples_tsb", label = "Minimum samples", min = 0, max = 100, value = 4, step = 1), # doesnt work if 0
          sliderInput(inputId = "min_redis_tsb", label = "Minimum redistributions", min = 0, max = 100, value = 3, step = 1),
          submitButton("Submit")
        ),
        mainPanel(
          uiOutput("select_tsb_plot"),
          plotOutput(outputId = "Treatment_specific_biomarkers")
        )
      ),

      ## TREATMENT SPECIFIC BIOMARKER SUBTPES
      tabPanel(
        "Treatment Specific Biomarker Subtyes",
        headerPanel("Treatment Specific Biomarker Subtypes"),
        sidebarPanel(
          sliderInput(inputId = "min_q_tsbs", label = "Select FDR threshhold", min = 0, max = 1, value = 1, step = 0.01),
          submitButton("Submit")
        ),
        mainPanel(
          uiOutput("select_tsb_subtype_plot"),
          plotOutput(outputId = "Treatment_specific_biomarker_subtypes")
        )
      ),
      tabPanel(
        "Predictive Biomarkers",
        headerPanel("Predictive Biomarkers"),
        sidebarPanel(
          sliderInput(inputId = "min_samples_pb", label = "Minimum samples", min = 0, max = 100, value = 4, step = 1), # doesnt work if 0
          sliderInput(inputId = "min_redis_pb", label = "Minimum redistributions", min = 0, max = 100, value = 3, step = 1),
          submitButton("Submit")
        ),
        mainPanel(
          uiOutput("select_pb_plot"),
          plotOutput(outputId = "Predictive_Biomarkers"),
        )
      ),
      tabPanel(
        "Predictive Biomarker Subtypes",
        headerPanel("Predictive Biomarker Subtypes"),
        sidebarPanel(
          sliderInput(inputId = "min_q_pbs", label = "Select FDR threshhold", min = 0, max = 1, value = 0.1, step = 0.01),
          submitButton("Submit")
        ),
        mainPanel(
          uiOutput("select_pb_subtype_plot"),
          plotOutput(outputId = "Predictive_Biomarker_Subtypes")
        )
      ),
      tabPanel(
        "Predictive Comparison",
        headerPanel("Predictive Comparison"),
        sidebarPanel(
          sliderInput(inputId = "min_q_predcomp", label = "Select FDR threshhold (treatment-specific)", min = 0, max = 1, value = 0.3, step = 0.01),
          submitButton("Submit")
        ),
        mainPanel(
          uiOutput("select_pc_plot"),
          plotOutput(outputId = "Predictive_Comp")
        )
      ),
      tabPanel(
        "Plot Example",
        headerPanel("Plot Example"),
        sidebarPanel(
          uiOutput("select_example_mutation"),
          uiOutput("select_example_readout"),
          uiOutput("select_example_subtype"),
          uiOutput("select_example_subtype_value"),
          submitButton("Submit")
        ),
        mainPanel(
          plotOutput(outputId = "Plot_Example"),
        )
      ),
      tabPanel(
        "Summary",
        headerPanel("Summary"),
        mainPanel(
          DT::DTOutput(outputId = "Summary_Table")
        )
      ),
      tabPanel(
        "Permutations",
        headerPanel("Permutations"),
        sidebarPanel(
          uiOutput("select_permut_readout"),
          uiOutput("select_permut_subtype"),
          sliderInput(inputId = "min_samples_per", label = "Select minimum number of samples", min = 1, max = 20, value = 10, step = 1),
          sliderInput(inputId = "min_redis_per", label = "Select minimum number of redistributions", min = 1, max = 20, value = 10, step = 1),
          sliderInput(inputId = "n_permut", label = "Select No permutations", min = 1, max = 2000, value = 10, step = 1),
          sliderInput(inputId = "fdr_permut", label = "Select fdr threshhold", min = 0, max = 1, value = 1, step = 0.01),
          sliderInput(inputId = "fdr_int_permut", label = "Select fdr int threshhold", min = 0, max = 1, value = 1, step = 0.01),
          submitButton("Submit")
        ),
        mainPanel(
          DT::DTOutput(outputId = "Permutation_Table")
        )
      ),
      tabPanel(
        "Bootstraps",
        headerPanel("Bootstraps"),
        sidebarPanel(
          uiOutput("select_boot_readout"),
          uiOutput("select_boot_subtype"),
          sliderInput(inputId = "min_samples_boot", label = "Select minimum number of samples", min = 1, max = 20, value = 10, step = 1),
          sliderInput(inputId = "min_redis_boot", label = "Select minimum number of redistributions", min = 1, max = 20, value = 10, step = 1),
          sliderInput(inputId = "n_boot", label = "Select No permutations", min = 1, max = 2000, value = 10, step = 1),
          sliderInput(inputId = "fdr_boot", label = "Select fdr threshhold", min = 0, max = 1, value = 1, step = 0.01),
          sliderInput(inputId = "fdr_int_boot", label = "Select fdr int threshhold", min = 0, max = 1, value = 1, step = 0.01),
          submitButton("Submit")
        ),
        mainPanel(
          DT::DTOutput(outputId = "Bootstrap_Table")
        )
      )
    )
  )
)
