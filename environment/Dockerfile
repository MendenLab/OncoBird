# hash:sha256:b2999a8c058033e8e6f09c97b0fa319fb79702b52f68e7d36819c7f967b729f4
FROM registry.codeocean.com/codeocean/r-studio:1.2.5019-r4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        apt-utils=1.6.14 \
        cmake=3.10.2-1ubuntu2.18.04.2 \
        r-base-dev=4.1.2-1.1804.0 \
        libnlopt0 \
        libnlopt-dev \
        pandoc \
        pandoc-citeproc \
        texlive-latex-base \
        lmodern \
        texlive-fonts-recommended \
        texlive-latex-recommended \
        texlive-latex-extra \
    && rm -rf /var/lib/apt/lists/*
    
RUN R -e "install.packages('nloptr')"

RUN R -e "install.packages('BiocManager')"

RUN R -e "BiocManager::install(c('tidyverse', 'data.table', 'survminer', 'reshape', 'ggfortify', 'pROC', 'formattable', 'gridExtra', 'h2o', 'glmnet', 'glmnetUtils', 'rjson', 'randomForest', 'caret', 'survival', 'Hmisc', 'randomForestSRC', 'here', 'survcomp', 'coin', 'partykit','DT','survminer','tinytex', 'rmarkdown'))"
RUN R -e "BiocManager::install(c('prettydoc','SummarizedExperiment'))"
RUN export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
RUN export CAIRO_INCLUDE_PATH=$CONDA_PREFIX/include/cairo
RUN export CAIRO_CFLAGS=-I$CONDA_PREFIX/include/cairo
RUN export CAIRO_LIBS=$CONDA_PREFIX/lib/cairo
RUN R -e "BiocManager::install(c('ComplexHeatmap'))"
#RUN R -e "tinytex::install_tinytex()"

RUN R -e "devtools::install_github(repo = 'aljoshoh/OncoBird', subdir = 'code/OncoBird')"
