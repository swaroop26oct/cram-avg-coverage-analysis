FROM ubuntu:22.04

LABEL maintainer="Nagaswaroop"
LABEL description="Jupyter environment for CRAM coverage analysis"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        samtools \
        python3 \
        python3-pip \
        wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir jupyter matplotlib

WORKDIR /data
COPY coverage_analysis.ipynb .

EXPOSE 8888
CMD ["python3", "-m", "jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
