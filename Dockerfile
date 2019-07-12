FROM centos:7

RUN yum -y update && \
    yum install -y vim && \
    yum install -y gcc-c++ && \
    yum install -y zlib-devel && \
    yum install -y epel-release && \
    yum install -y python34 && \
    yum install -y python34-setuptools && \
    yum install -y wget && \
    yum install -y R && \
    yum install -y bzip2 && \
    yum install -y curl-devel && \
    yum install -y cairo-devel && \
    yum install -y openssl-devel && \
    yum install -y libstdc++-static

RUN easy_install-3.4 pip

WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN tar -xf samtools-1.9.tar.bz2
WORKDIR /opt/samtools-1.9
RUN ./configure --without-curses
RUN make
RUN cp samtools /usr/bin/

WORKDIR /opt/Anaquin
COPY . /opt/Anaquin
COPY data resources

WORKDIR /opt/Anaquin
RUN make

RUN Rscript -e 'install.packages(c("knitr", "ggplot2", "ROCR", "methods", "plyr", "data.table"), repos="https://cran.rstudio.com")'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org"); BiocManager::install("qvalue", version = "3.9", ask=F)'
RUN wget https://github.com/sequinstandards/RAnaquin/archive/master.zip
RUN unzip master.zip
RUN R CMD build RAnaquin-master
RUN R CMD INSTALL Anaquin_*.tar.gz

ENV PATH="/opt/Anaquin:${PATH}"
