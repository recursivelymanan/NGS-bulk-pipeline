FROM ubuntu:22.04



# MAINTAINER Manan Chopra <manan.zchopra@gmail.com>

# RUN apt-get update && \
#     apt-get install -y wget bzip2

# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh && \
#     bash Miniconda3-latest-MacOSX-arm64.sh -b -p /opt/conda && \
#     rm Miniconda3-latest-MacOSX-arm64.sh && \
#     ~/miniconda3/bin/conda init bash