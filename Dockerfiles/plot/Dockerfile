FROM bioconductor/bioconductor_docker:RELEASE_3_12  

# install argument parsing tools
RUN R -e "install.packages(c('argparse','optparse','this.path','here','import'))"
# install data processing R packages
RUN R -e "install.packages(c('R.utils','data.table','reshape2','anndata'))" 
# install visualization tools
RUN R -e "install.packages(c('ggplot2','GGally','ggpubr','egg','UpSetR','VennDiagram','tidyverse'))"