FROM ubuntu:16.04

MAINTAINER dekkerlab

#  a lot of it is stollen from :
# https://hub.docker.com/r/genomicpariscentre/bioperl/dockerfile/

# install a bunch of system-wide stuff
# and figure out later what's essential ...
RUN apt-get update --fix-missing && \
 apt-get install -q -y apt-utils wget \
 curl bzip2 pkg-config \
 libbz2-dev git build-essential \
 zlib1g-dev locales vim fontconfig \
 perl expat libexpat-dev\
 ttf-dejavu cpanminus libgd-perl

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8


# Install perl modules 
RUN apt-get install -y cpanminus

RUN cpanm CPAN::Meta \
     readline \
     Term::ReadKey \
     YAML \
     Digest::SHA \
     Module::Build \
     ExtUtils::MakeMaker \
     Test::More \
     Data::Stag \
     Config::Simple \
     Statistics::Lite \
     Statistics::Descriptive 

RUN apt-get install --yes \
 libarchive-zip-perl


# Install GD
RUN apt-get remove --yes libgd-gd2-perl

RUN apt-get install --yes \
 libgd2-noxpm-dev

RUN cpanm GD \
 GD::Graph \
 GD::Graph::smoothlines 

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

RUN conda install pysam


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



# global install ...
RUN perl ../Build.PL
RUN ../Build
RUN ../Build install


RUN ../Build distclean
# After installing the module, you should be free to run anything in scripts/ e.g.

# $ perl scripts/heatmap.pl