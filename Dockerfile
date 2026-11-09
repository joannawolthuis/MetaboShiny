# Tested dependency baseline; no floating package downloads during this build.
# The complete legacy source recipe is preserved in docker/Dockerfile.legacy.
ARG DEPENDENCY_IMAGE=jcwolthuis/metaboshiny@sha256:f5f49da426c1b31caa7674648dc91172cb32fd3eb734ce2dcc18ca4e815e852c
FROM ${DEPENDENCY_IMAGE}

# magick 2.6.0 in the baseline needs the Ubuntu Focal ImageMagick C++ runtime.
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libmagick++-6.q16-8 \
    && rm -rf /var/lib/apt/lists/*

ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    METABOSHINY_PORT=8080 METABOSHINY_HOST=0.0.0.0 METABOSHINY_PARALLEL=0
WORKDIR /opt/metaboshiny
COPY DESCRIPTION docker/check_dependencies.R ./
RUN Rscript --vanilla check_dependencies.R
COPY . ./
RUN R CMD INSTALL --no-multiarch . \
    && Rscript --vanilla scripts/test_startup.R \
    && Rscript --vanilla docker/smoke_ui.R \
    && Rscript --vanilla docker/smoke_data.R \
    && Rscript --vanilla -e 'write.csv(installed.packages()[,c("Package","Version")], "/opt/metaboshiny/package-versions.csv", row.names=FALSE); writeLines(capture.output(sessionInfo()), "/opt/metaboshiny/build-sessionInfo.txt")'
EXPOSE 8080
HEALTHCHECK --interval=15s --timeout=5s --start-period=120s --retries=4 \
    CMD curl --max-time 4 --fail --silent --output /dev/null "http://127.0.0.1:${METABOSHINY_PORT}/" || exit 1
ENTRYPOINT ["/bin/sh", "/opt/metaboshiny/docker/entrypoint.sh"]
