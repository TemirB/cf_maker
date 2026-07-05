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
# Paths
# ----------------------------
set(INCLUDE_DIR    ${CMAKE_SOURCE_DIR}/include)
set(SRC_DIR        ${CMAKE_SOURCE_DIR}/src)
set(EXTERNAL_DIR   ${CMAKE_SOURCE_DIR}/external)

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
# nlohmann/json
# ----------------------------
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install OFF CACHE INTERNAL "")

add_subdirectory(
    ${EXTERNAL_DIR}/nlohmann_json
    ${CMAKE_BINARY_DIR}/external/nlohmann_json
    EXCLUDE_FROM_ALL
)

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

target_compile_features(main PRIVATE cxx_std_17)

# ----------------------------
# Include paths
# ----------------------------
target_include_directories(main PRIVATE
    ${INCLUDE_DIR}
)

# ----------------------------
# Compile options
# ----------------------------
target_compile_options(main PRIVATE
    -Wall
    -Wextra
    -O2
)

target_compile_definitions(main PRIVATE
    ROOT_STRICT_OWNERSHIP
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
    nlohmann_json::nlohmann_json
)

# ----------------------------
# Install
# ----------------------------
install(TARGETS main DESTINATION bin)