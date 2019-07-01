FROM continuumio/miniconda3
MAINTAINER Wanson Choi <wschoi.bhlab@gmail.com>

### Anaconda 
# Pandas
RUN conda install -y pandas=0.23.0
# R
RUN conda install -y -c conda-forge xorg-libx11
RUN conda install -y -c r r
RUN Rscript -e "install.packages(c('gplots', 'RColorBrewer', 'shape'), dependencies=TRUE, repos='http://healthstat.snu.ac.kr/CRAN/')"
# Java
RUN conda install -y -c anaconda openjdk
# RUN conda install -y -c cyclus java-jdk

### Main Sources
RUN mkdir /HATK
# ADD HATK.py /HATK
# ADD IMGT2Seq /HATK/IMGT2Seq
# ADD bMarkerGenerator /HATK/bMarkerGenerator
# ADD NomenCleaner /HATK/NomenCleaner
# ADD HLA2HPED /HATK/HLA2HPED
# ADD HLA_Analysis /HATK/HLA_Analysis
# ADD HLA_Manhattan /HATK/HLA_Manhattan
# ADD HLA_Heatmap /HATK/HLA_Heatmap
# ADD src /HATK/src
ADD dependency /HATK/dependency

# RUN git clone https://github.com/WansonChoi/HATK.git

WORKDIR /HATK


RUN echo "HATK_Docker has been prepared."


# CMD python MakeReference_v3.py
