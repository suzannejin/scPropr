FROM bioconductor/bioconductor_docker:RELEASE_3_12  

# install basic R tools
ADD basic.R /tmp/basic.R 
RUN Rscript /tmp/basic.R

# install propr & independencies
ADD propr.R /tmp/propr.R
RUN Rscript /tmp/propr.R

# install normalization methods
ADD normalize.R /tmp/normalize.R
RUN Rscript /tmp/normalize.R

# install sanity
RUN git clone https://github.com/jmbreda/Sanity.git \
    && apt-get update \
    && apt-get install -y libgomp1 \
    && cd Sanity/src \
    && make \
    && cp ../bin/Sanity /bin/.
