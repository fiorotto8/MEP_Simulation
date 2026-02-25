cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(gpd3d)

# ------------------------------------------------
# Options
# ------------------------------------------------
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
add_compile_options(-O3)

# ------------------------------------------------
# Geant4
# ------------------------------------------------
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

# ------------------------------------------------
# ROOT
# ------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
include(${ROOT_USE_FILE})

# ------------------------------------------------
# Boost
# ------------------------------------------------
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# ------------------------------------------------
# CFITSIO (via pkg-config)
# ------------------------------------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED cfitsio)

# Include directories
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${CFITSIO_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIRS}
)

# Link directories (if any)
link_directories(${CFITSIO_LIBRARY_DIRS})

# ------------------------------------------------
# Sources / headers / macros
# ------------------------------------------------
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB_RECURSE macros  RELATIVE ${PROJECT_SOURCE_DIR} macros/*)

# Copy macros to build dir (out-of-source build friendly)
foreach(_file ${macros})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
  )
endforeach()

# ------------------------------------------------
# Build
# ------------------------------------------------
add_executable(gpd3d main.cc ${sources} ${headers})

target_include_directories(gpd3d PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${CFITSIO_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(gpd3d
  ${Geant4_LIBRARIES}
  ${CFITSIO_LIBRARIES}
  ${ROOT_LIBRARIES}
)

install(TARGETS gpd3d DESTINATION bin)
