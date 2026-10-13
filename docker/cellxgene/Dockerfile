# cellxgene pipeline image: worker (aggregate.py), planner (plan.py), publisher (publish.py)
# and driver (driver.py). Build with the repo root as context:
#   docker build -f docker/cellxgene/Dockerfile .
# WORKER_IMAGE_DIGEST is deliberately not set here: a digest cannot be known before the image
# is pushed. The driver must set it on every job from the pushed digest; aggregate.py records
# it in _SUCCESS.json.
FROM python:3.13-slim

ARG HTSLIB_VER=1.22.1

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# build deps are removed in the same layer; the runtime libs htslib links against stay
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libbz2-dev liblzma-dev zlib1g-dev libssl-dev libcurl4-openssl-dev libdeflate-dev \
        bzip2 gcc make curl ca-certificates \
    && mkdir /opt/htslib && cd /opt/htslib \
    && curl -fsSLO https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VER}.tar.bz2 && cd htslib-${HTSLIB_VER} \
    && ./configure --with-libdeflate && make -j"$(nproc)" && make install && ldconfig \
    && cd / && rm -rf /opt/htslib \
    && apt-get purge -y --auto-remove gcc make curl bzip2 libbz2-dev liblzma-dev zlib1g-dev libssl-dev libcurl4-openssl-dev libdeflate-dev \
    && apt-get install -y --no-install-recommends libbz2-1.0 liblzma5 zlib1g libssl3 libcurl4 libdeflate0 \
    && rm -rf /var/lib/apt/lists/*

COPY docker/cellxgene/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

# rollup.py imports sumstat_utils from the parent of its own directory
WORKDIR /app
COPY scripts/sumstat_utils.py scripts/
COPY scripts/cellxgene scripts/cellxgene
COPY tests/cellxgene tests/cellxgene

RUN useradd --uid 10001 --create-home worker
USER worker
