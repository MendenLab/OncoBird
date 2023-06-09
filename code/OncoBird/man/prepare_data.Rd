% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/se_wrapper.R
\name{prepare_data}
\alias{prepare_data}
\title{Prepare datasets for run}
\usage{
prepare_data(
  data = NULL,
  data_clinical = NULL,
  data_mutations = NULL,
  vars = c("OS", "PFS", "RFS", "DFS", "ORR"),
  med_impute_muts = TRUE,
  remove_clin = NA,
  treatment_column = "treatment",
  sample_column = "sample",
  mutation_columns = NULL
)
}
\arguments{
\item{data}{dataframe containing clinical and mutational variables}

\item{data_clinical}{alternatively to 'data', dataframe of clinical
variables}

\item{data_mutations}{alternatively to 'data', dataframe of mutational
variables}

\item{vars}{character vector of clinical endpoints to be analysed
(OS/PFS/ORR/RFS/DFS)}

\item{med_impute_muts}{logical if missing clinical data is imputed by
the median}

\item{remove_clin}{character vector of clinical variables to be excluded
from the downstream analysis}

\item{treatment_column}{character to specify the treatment column}

\item{sample_column}{character to specify the sample identifier column}

\item{mutation_columns}{character vector to specifiy the mutational column,
if only 'data' is supplied}
}
\value{
SE object
}
\description{
Initiate SE object for downstream analysis
}
\examples{
{
  sim_data <- sim_data()

  se <- prepare_data(
    data = sim_data,
    vars = c("OS", "DFS"),
    med_impute_muts = TRUE,
    remove_clin = NA,
    sample_column = "PatientID",
    treatment_column = "Adj",
    mutation_columns = colnames(sim_data)[14:ncol(sim_data)]
  )
}

}
