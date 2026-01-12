cmake_minimum_required(VERSION 3.16)

project(cf_maker_new LANGUAGES CXX)

# ----------------------------
# Global settings
# ----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------
# ROOT
# ----------------------------
find_package(ROOT REQUIRED COMPONENTS
    Core RIO Hist Graf Graf3d Gpad Tree TreePlayer Physics MathCore
)
include(${ROOT_USE_FILE})

# ----------------------------
# OpenSSL (for SHA256)
# ----------------------------
find_package(OpenSSL REQUIRED)

# ----------------------------
# Paths
# ----------------------------
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(SRC_DIR     ${CMAKE_SOURCE_DIR}/src)
set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/external)

# ----------------------------
# Sources
# ----------------------------
file(GLOB_RECURSE PROJECT_SOURCES
    ${SRC_DIR}/*.cpp
)

set(MAIN_SOURCE ${CMAKE_SOURCE_DIR}/main.cpp)

# ----------------------------
# Executable
# ----------------------------
add_executable(main
    ${MAIN_SOURCE}
    ${PROJECT_SOURCES}
)

# ----------------------------
# Include paths
# ----------------------------
target_include_directories(main PRIVATE
    ${INCLUDE_DIR}
    ${EXTERNAL_DIR}
)

# ----------------------------
# Compile options
# ----------------------------
target_compile_options(main PRIVATE
    -Wall
    -Wextra
    -O2
    -DROOT_STRICT_OWNERSHIP
)

# ----------------------------
# Link libraries
# ----------------------------
target_link_libraries(main PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::Physics
    ROOT::MathCore
    OpenSSL::Crypto
)

# ----------------------------
# Install
# ----------------------------
install(TARGETS main DESTINATION bin)
