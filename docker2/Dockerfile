# Set the base image to the BioPerl prebuilt prerequisites image
# for building from source
FROM bioperl/bioperl-deps

# # File Author / Maintainer
# LABEL maintainer="Hilmar Lapp <hlapp@drycafe.net>"

MAINTAINER dekkerlab

USER root

# # Set the working directory to where we will install bioperl
# # in the builder stage
# WORKDIR /bioperl

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN apt-get update --yes && \
    apt-get install --yes git

# install imagemagick
RUN apt-get install --yes imagemagick

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda3 -b && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda3/bin:${PATH}

# Install conda dependencies
ADD cworld_environment.yml /
ADD VERSION /
RUN pwd
RUN conda config --set always_yes yes --set changeps1 no && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --get && \
    conda update -q conda && \
    conda info -a && \
    conda env update -q -n root --file cworld_environment.yml && \
    conda clean --tarballs --index-cache --lock

# RUN conda install pysam


#export MKL OMP etc ...


# # get the version of the GDlib:
# RUN perl -MGD -e 'print $GD::VERSION ."\n";'
# # local isntall
# perl Build.PL
# ./Build
# ./Build install --install_base /your/custom/dir
# (ensure /your/custom/dir is added to your PERL5LIB path)
# e.g.
# ./Build install --install_base ~/perl5
# # then in .bashrc
# export PERL5LIB=${PERL5LIB}:/home/<yourusername>/perl5/lib/perl5

WORKDIR /cworld-dekker
ADD Build.PL .
ADD lib ./lib
ADD scripts ./scripts
ADD MANIFEST .

# global install ...
RUN perl ./Build.PL
RUN ./Build
RUN ./Build install
RUN ./Build install --install_base /perl5
ENV PERL5LIB=${PERL5LIB}:/perl5/lib/perl5

RUN ./Build distclean
# After installing the module, you should be free to run anything in scripts/ e.g.
# $ perl scripts/heatmap.pl

# WORKDIR /root
# it should be more civilized, but for now, let's hope it just works
