# cellxgene pipeline image: worker (aggregate.py), planner (plan.py), publisher (publish.py)
# and driver (driver.py). Build with the repo root as context and a genetics-results-db checkout
# as the named context `db` (the driver's load step runs that repo's BigQuery loader):
#   docker build -f docker/cellxgene/Dockerfile --build-context db=../genetics-results-db .
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

# rollup.py and publish.py put scripts/ on sys.path for sumstat_utils and peak_utils (which
# itself imports sumstat_utils); every other import is stdlib, a pip package or scripts/cellxgene
WORKDIR /app
COPY scripts/sumstat_utils.py scripts/peak_utils.py scripts/
COPY scripts/cellxgene scripts/cellxgene
COPY tests/cellxgene tests/cellxgene

# only what db's load_cellxgene.sh runs; scripts/cellxgene/load_cellxgene.sh finds it here
COPY --from=db --chmod=755 scripts/load_cellxgene.sh /opt/genetics-results-db/scripts/
COPY --from=db scripts/load_data.py /opt/genetics-results-db/scripts/
COPY --from=db scripts/lib/common.sh /opt/genetics-results-db/scripts/lib/
ENV GENETICS_RESULTS_DB=/opt/genetics-results-db

# a missing module or dependency fails the build, not the first real run. --help runs argparse
# before any network or credential use; the loader is run by path as load_cellxgene.sh does
RUN for s in driver publish plan aggregate rollup; do python scripts/cellxgene/$s.py --help >/dev/null || exit 1; done \
    && python /opt/genetics-results-db/scripts/load_data.py --help >/dev/null

RUN useradd --uid 10001 --create-home worker
USER worker
