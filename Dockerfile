FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-lc"]

# ------------------------------------------------
# System dependencies
# ------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git wget curl pkg-config \
    python3 python3-dev python3-pip \
    gcc g++ gfortran \
    libssl-dev \
    libx11-dev libxext-dev libxft-dev libxpm-dev libxmu-dev libxrender-dev \
    libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev \
    libz-dev libbz2-dev liblzma-dev \
    libpcre3-dev \
    libtbb-dev \
    libfftw3-dev \
    libgsl-dev \
    libxml2-dev \
    libjpeg-dev libpng-dev libtiff-dev \
    libxerces-c-dev \
    libreadline-dev libncurses5-dev \
    libboost-all-dev libcfitsio-dev \
    python3 python3-dev python3-pip \
    python3-matplotlib python3-numpy python3-tk fonts-dejavu-core \
    qtbase5-dev \
    qtchooser \
    qt5-qmake \
    qtbase5-dev-tools \
    libqt5opengl5-dev \
    libxkbcommon-x11-0 \
    libqt5widgets5 libqt5gui5 libqt5core5a libqt5opengl5 \  
    libgl1-mesa-dri mesa-utils xfonts-base \
    libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-xinerama0 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt

# ------------------------------------------------
# ROOT 6.20/04
# ------------------------------------------------
ARG ROOT_VER=6.20.04

RUN wget -q https://root.cern/download/root_v${ROOT_VER}.source.tar.gz \
 && tar -xzf root_v${ROOT_VER}.source.tar.gz \
 && rm root_v${ROOT_VER}.source.tar.gz

RUN mkdir -p root-build \
 && cd root-build \
 && cmake \
    -DCMAKE_INSTALL_PREFIX=/opt/root \
    -DCMAKE_BUILD_TYPE=Release \
    -Dgnuinstall=ON \
    -Dpython3=ON \
    -Dminuit2=ON \
    -Droofit=ON \
    -Dx11=ON \
    -Dopengl=ON \
    ../root-${ROOT_VER} \
 && cmake --build . -- -j$(nproc) \
 && cmake --build . --target install -- -j$(nproc) \
 && test -f /opt/root/bin/thisroot.sh

RUN echo 'if [ -f /opt/root/bin/thisroot.sh ]; then source /opt/root/bin/thisroot.sh; fi' > /etc/profile.d/root.sh

# ------------------------------------------------
# Geant4 10.5.1
# ------------------------------------------------
ARG G4_VER=10.5.1
ARG G4_BASE=https://gitlab.cern.ch/geant4/geant4/-/archive/v${G4_VER}

RUN wget -q ${G4_BASE}/geant4-v${G4_VER}.tar.gz \
 && tar -xzf geant4-v${G4_VER}.tar.gz \
 && rm geant4-v${G4_VER}.tar.gz

RUN mkdir -p geant4-build \
 && cd geant4-build \
 && cmake \
    -DCMAKE_INSTALL_PREFIX=/opt/geant4 \
    -DCMAKE_BUILD_TYPE=Release \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_USE_OPENGL_X11=ON \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_USE_GDML=ON \
    ../geant4-v${G4_VER} \
 && cmake --build . -- -j$(nproc) \
 && cmake --build . --target install -- -j10 \
 && test -f /opt/geant4/bin/geant4.sh

RUN echo 'if [ -f /opt/geant4/bin/geant4.sh ]; then source /opt/geant4/bin/geant4.sh; fi' > /etc/profile.d/geant4.sh

# Start an interactive login shell so /etc/profile.d/*.sh is sourced
CMD ["/bin/bash", "-l"]
