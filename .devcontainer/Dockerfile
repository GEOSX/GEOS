# Define `GEOS_TPL_TAG` (as defined in `.github/workflows/ci_tests.yml`) to use the proper version of the TPLs.
# If you are building the current `Dockerfile` to use codespaces in the cloud,
# you most likely want to define this variable in the `.devcontainer/devcontainer.json` file.
ARG GEOS_TPL_TAG
FROM docker.io/geosx/ubuntu20.04-gcc10:${GEOS_TPL_TAG}

RUN apt-get update
RUN apt-get install -y \
    ca-certificates \
    curl \
    git-lfs \
    vim \
    ninja-build \
    gdb \
# for paraview headless
    libcurl4 \
    libglapi-mesa
