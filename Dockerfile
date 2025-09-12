FROM python:3.13-slim

ARG HTSLIB_VER=1.22.1

# unbuffered output
ENV PYTHONUNBUFFERED=1
# don't write .pyc files
ENV PYTHONDONTWRITEBYTECODE=1

RUN apt-get update \
    && apt-get install -y \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libpcre2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
	libdeflate-dev \
    bzip2 \
    gcc \
    g++ \
    make \
    curl \
    datamash \
    && pip install uv --upgrade \
    && rm -rf /var/lib/apt/lists/*

# htslib (tabix)
WORKDIR /opt/htslib
RUN curl -LO https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 && \
    tar -xvjf htslib-${HTSLIB_VER}.tar.bz2 && cd htslib-${HTSLIB_VER} && \
    ./configure --with-libdeflate && make && make install && cd .. && rm -rf htslib-${HTSLIB_VER}*

WORKDIR /tmp
COPY requirements.txt .
RUN uv pip install --system -r requirements.txt \
    && rm -f requirements.txt
