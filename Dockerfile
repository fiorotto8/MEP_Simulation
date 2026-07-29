# syntax=docker/dockerfile:1.7

FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Current stable versions as of 2026-07-28.
# They can be overridden with --build-arg.
ARG PYTHON_VER=3.14.6
ARG ROOT_VER=6.40.02
ARG GEANT4_VER=11.4.2

# Leave unset to use all available CPU cores.
ARG BUILD_JOBS

ENV PYTHON_PREFIX=/opt/python \
    ROOTSYS=/opt/root \
    GEANT4_PREFIX=/opt/geant4 \
    PATH=/opt/python/bin:/opt/root/bin:/opt/geant4/bin:${PATH} \
    LD_LIBRARY_PATH=/opt/python/lib:/opt/root/lib:/opt/geant4/lib \
    CMAKE_PREFIX_PATH=/opt/root:/opt/geant4 \
    PYTHONPATH=/opt/root/lib \
    PYTHONUNBUFFERED=1

# ------------------------------------------------
# System dependencies
# ------------------------------------------------
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    ca-certificates \
    binutils \
    build-essential \
    cmake \
    dpkg-dev \
    git \
    ninja-build \
    wget \
    curl \
    pkg-config \
    xz-utils \
    \
    gcc \
    g++ \
    gfortran \
    python3 \
    \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev \
    libsqlite3-dev \
    libffi-dev \
    libgdbm-dev \
    libgdbm-compat-dev \
    libncurses-dev \
    tk-dev \
    uuid-dev \
    \
    libx11-dev \
    libxext-dev \
    libxft-dev \
    libxpm-dev \
    libxmu-dev \
    libxrender-dev \
    libxrandr-dev \
    libxinerama-dev \
    libxcursor-dev \
    libxi-dev \
    libxfixes-dev \
    \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libgl1-mesa-dri \
    freeglut3-dev \
    libglew-dev \
    libgl2ps-dev \
    mesa-utils \
    \
    libpcre3-dev \
    libtbb-dev \
    libvdt-dev \
    libfftw3-dev \
    libgsl-dev \
    libxml2-dev \
    libexpat1-dev \
    libxerces-c-dev \
    \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    libgif-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libftgl-dev \
    \
    libboost-all-dev \
    libcfitsio-dev \
    nlohmann-json3-dev \
    \
    libcurl4-openssl-dev \
    liblz4-dev \
    libzstd-dev \
    libxxhash-dev \
    \
    qt6-base-dev \
    qt6-base-dev-tools \
    libqt6opengl6-dev \
    qt6-qpa-plugins \
    \
    libxkbcommon-x11-0 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-xinerama0 \
    \
    fonts-dejavu-core \
    xfonts-base \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# ------------------------------------------------
# Python 3.14.6
#
# Installed independently from Ubuntu's system
# Python so distribution tools continue to use
# /usr/bin/python3.
# ------------------------------------------------
RUN set -eux; \
    JOBS="${BUILD_JOBS:-$(nproc)}"; \
    cd /tmp; \
    wget -q "https://www.python.org/ftp/python/${PYTHON_VER}/Python-${PYTHON_VER}.tgz"; \
    tar -xzf "Python-${PYTHON_VER}.tgz"; \
    cd "Python-${PYTHON_VER}"; \
    ./configure \
        --prefix="${PYTHON_PREFIX}" \
        --enable-shared \
        --with-ensurepip=install; \
    make -j "${JOBS}"; \
    make install; \
    echo "${PYTHON_PREFIX}/lib" > /etc/ld.so.conf.d/python.conf; \
    ldconfig; \
    "${PYTHON_PREFIX}/bin/python3" -m pip install \
        --no-cache-dir \
        --upgrade \
        pip \
        setuptools \
        wheel \
        numpy \
        matplotlib; \
    "${PYTHON_PREFIX}/bin/python3" --version; \
    rm -rf \
        "/tmp/Python-${PYTHON_VER}" \
        "/tmp/Python-${PYTHON_VER}.tgz"

# ------------------------------------------------
# ROOT 6.40.02
# ------------------------------------------------
RUN set -eux; \
    JOBS="${BUILD_JOBS:-$(nproc)}"; \
    cd /tmp; \
    wget -q \
        "https://root.cern/download/root_v${ROOT_VER}.source.tar.gz"; \
    tar -xzf "root_v${ROOT_VER}.source.tar.gz"; \
    \
    cmake \
        -S "/tmp/root-${ROOT_VER}" \
        -B /tmp/root-build \
        -G Ninja \
        -DCMAKE_INSTALL_PREFIX="${ROOTSYS}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -Dgnuinstall=ON \
        -Drpath=ON \
        -Dpyroot=ON \
        -Dtmva=ON \
        -Dtmva-pymva=OFF \
        -Droofit=ON \
        -Dx11=ON \
        -Dopengl=ON \
        -DPython3_ROOT_DIR="${PYTHON_PREFIX}" \
        -DPython3_EXECUTABLE="${PYTHON_PREFIX}/bin/python3"; \
    \
    cmake --build /tmp/root-build --parallel "${JOBS}"; \
    cmake --install /tmp/root-build; \
    test -f "${ROOTSYS}/bin/thisroot.sh"; \
    \
    rm -rf \
        "/tmp/root-${ROOT_VER}" \
        /tmp/root-build \
        "/tmp/root_v${ROOT_VER}.source.tar.gz"
        
# ------------------------------------------------
# Geant4 11.4.2
# ------------------------------------------------
RUN set -eux; \
    JOBS="${BUILD_JOBS:-$(nproc)}"; \
    cd /tmp; \
    wget -q \
        "https://gitlab.cern.ch/geant4/geant4/-/archive/v${GEANT4_VER}/geant4-v${GEANT4_VER}.tar.gz"; \
    tar -xzf "geant4-v${GEANT4_VER}.tar.gz"; \
    \
    cmake \
        -S "/tmp/geant4-v${GEANT4_VER}" \
        -B /tmp/geant4-build \
        -G Ninja \
        -DCMAKE_INSTALL_PREFIX="${GEANT4_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DGEANT4_BUILD_MULTITHREADED=ON \
        -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_USE_GDML=ON \
        -DGEANT4_USE_OPENGL_X11=ON \
        -DGEANT4_USE_QT=ON; \
    \
    cmake --build /tmp/geant4-build --parallel "${JOBS}"; \
    cmake --install /tmp/geant4-build; \
    test -f "${GEANT4_PREFIX}/bin/geant4.sh"; \
    \
    rm -rf \
        "/tmp/geant4-v${GEANT4_VER}" \
        /tmp/geant4-build \
        "/tmp/geant4-v${GEANT4_VER}.tar.gz"

# ------------------------------------------------
# Environment configuration
# ------------------------------------------------
RUN printf '%s\n' \
    '#!/usr/bin/env bash' \
    '' \
    'if [[ -f /opt/root/bin/thisroot.sh ]]; then' \
    '    source /opt/root/bin/thisroot.sh' \
    'fi' \
    '' \
    'if [[ -f /opt/geant4/bin/geant4.sh ]]; then' \
    '    source /opt/geant4/bin/geant4.sh' \
    'fi' \
    > /etc/profile.d/hep.sh \
 && chmod 0644 /etc/profile.d/hep.sh

# Source the environments even when the container is
# invoked with a direct command rather than a login shell.
RUN printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -e' \
    '' \
    'source /opt/root/bin/thisroot.sh' \
    'source /opt/geant4/bin/geant4.sh' \
    '' \
    'exec "$@"' \
    > /usr/local/bin/hep-entrypoint \
 && chmod +x /usr/local/bin/hep-entrypoint

# ------------------------------------------------
# Installation checks
# ------------------------------------------------
RUN source "${ROOTSYS}/bin/thisroot.sh" \
 && source "${GEANT4_PREFIX}/bin/geant4.sh" \
 && python3 --version \
 && root-config --version \
 && geant4-config --version \
 && python3 -c \
    'import ROOT; print("PyROOT:", ROOT.gROOT.GetVersion())'

WORKDIR /workspace

ENTRYPOINT ["/usr/local/bin/hep-entrypoint"]
CMD ["/bin/bash", "-l"]