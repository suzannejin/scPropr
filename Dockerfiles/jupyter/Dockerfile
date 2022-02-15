FROM bioconductor/bioconductor_docker:RELEASE_3_12  

# =============== #
# INSTALL JUPYTER #
# =============== #

## Python 3 installations
RUN apt-get update \
	&& apt-get install -y software-properties-common \
	&& add-apt-repository universe \
	&& apt-get update \
	&& apt-get -y --no-install-recommends install python3 python3-dev \
	&& curl https://bootstrap.pypa.io/get-pip.py --output get-pip.py \
	&& python3 get-pip.py \
	&& pip3 install wheel \
	## Install sklearn, pandas, jupyter notebook on python
	&& pip3 install numpy \
	sklearn \
	pandas \
	pyyaml \
	cwltool \
	jupyter \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& rm -rf get-pip.py 

# Install all OS dependencies for fully functional notebook server
RUN apt-get update && apt-get install -yq --no-install-recommends \
    build-essential \
    inkscape \
    libsm6 \
    libxext-dev \
    libxrender1 \
    lmodern \
    netcat \
    # ---- nbconvert dependencies ----
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-plain-generic \
    # ----
    tzdata \
    unzip \
    nano-tiny \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install jupyter R kernel
RUN apt-get update \
    && apt-get install -y --no-install-recommends libzmq5
RUN R -e "install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'uuid', 'digest', 'IRkernel'))" && \
    R -e "IRkernel::installspec(user = FALSE)"

# Install jupyter bash kernel
RUN pip3 install bash_kernel && python3 -m bash_kernel.install

# install jupyter extensions
# and use useful extensions: Codefolding, ExecuteTime, Highlighter, limit output, Scratchpad, Skip-Traceback, Spellchecker, Table of Content (toc2), Variable inspector
RUN pip install jupyter_contrib_nbextensions
RUN pip install jupyter_nbextensions_configurator
RUN jupyter contrib nbextension install 
RUN jupyter nbextensions_configurator enable 
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension install https://rawgit.com/jfbercher/small_nbextensions/master/highlighter.zip 
RUN jupyter nbextension enable highlighter/highlighter

# install memory resource usage checker
RUN pip install jupyter-resource-usage


# ==================== #
# INSTALL REQUIREMENTS #
# ==================== #

# install argument parsing tools
RUN R -e "install.packages(c('argparse','optparse'))"
# install data processing R packages
RUN R -e "install.packages(c('R.utils','data.table','reshape2','anndata'))" 
# install visualization tools
RUN R -e "install.packages(c('ggplot2','GGally','ggpubr','egg','UpSetR','VennDiagram','tidyverse'))"
# install partial correlation packages
RUN R -e "install.packages(c('ppcor','corpcor'))"
# install zero imputation packages
RUN R -e "install.packages(c('zCompositions'))"
# install proportionality package
RUN R -e "devtools::install_github('tpq/propr@83dc0c54caac5f4e5b334fe31b953cc2469cf155')"
# install R single cell package
RUN R -e "install.packages(c('Seurat'))"
# install python single cell packages
RUN pip3 install scanpy \
    && pip3 install anndata
# install normalization packages
# scran package depends on many other packages, so it's a bit installation
RUN R -e "install.packages(c('sctransform'))"
RUN R -e "BiocManager::install(c('edgeR','scran'))"