FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB pacbio pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/pacbio-preprocess-1.1/bin:$PATH

