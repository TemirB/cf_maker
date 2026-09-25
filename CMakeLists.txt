cmake_minimum_required(VERSION 3.20)

project(cf_maker_new LANGUAGES CXX)

# ----------------------------
# Global settings
# ----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(CF_MAKER_SANITIZE "Build with Address/Undefined sanitizers" OFF)

# ----------------------------
# ccache
# ----------------------------
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}" CACHE STRING "Compiler launcher" FORCE)
    message(STATUS "ccache enabled: ${CCACHE_PROGRAM}")
endif()

# ----------------------------
# Paths
# ----------------------------
set(SRC_DIR        ${CMAKE_SOURCE_DIR}/src)
set(EXTERNAL_DIR   ${CMAKE_SOURCE_DIR}/external)

# ----------------------------
# ROOT
# ----------------------------
find_package(ROOT REQUIRED COMPONENTS
    Core RIO Hist Graf Graf3d Gpad Tree TreePlayer Physics MathCore
)

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
    CONFIGURE_DEPENDS
    ${SRC_DIR}/*.cpp
)

add_executable(main ${PROJECT_SOURCES})
target_compile_features(main PRIVATE cxx_std_17)

# ----------------------------
# Warnings
# ----------------------------
add_library(project_warnings INTERFACE)
target_compile_options(project_warnings INTERFACE
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wno-sign-conversion
    -Wnon-virtual-dtor
    -Wnull-dereference
    -Wformat=2
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(project_warnings INTERFACE -Wold-style-cast)
endif()

# ----------------------------
# Sanitizers
# ----------------------------
if(CF_MAKER_SANITIZE)
    target_compile_options(main PRIVATE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
    )
    target_link_options(main PRIVATE -fsanitize=address,undefined)
endif()

# ----------------------------
# Include paths
# ----------------------------
target_include_directories(main PRIVATE ${SRC_DIR})
target_include_directories(main SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS})

get_target_property(JSON_INCLUDE_DIRS nlohmann_json::nlohmann_json INTERFACE_INCLUDE_DIRECTORIES)
if(JSON_INCLUDE_DIRS)
    target_include_directories(main SYSTEM PRIVATE ${JSON_INCLUDE_DIRS})
endif()

# ----------------------------
# Compile options and definitions
# ----------------------------
target_compile_definitions(main PRIVATE
    ROOT_STRICT_OWNERSHIP
)

# ----------------------------
# Link libraries
# ----------------------------
target_link_libraries(main PRIVATE
    project_warnings
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

add_executable(generate_gaussian tools/generate_gaussian.cpp)
target_include_directories(generate_gaussian PRIVATE ${SRC_DIR})
target_include_directories(generate_gaussian SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(generate_gaussian PRIVATE
    project_warnings ROOT::Core ROOT::RIO ROOT::Hist nlohmann_json::nlohmann_json
)

# ----------------------------
# Tests
# ----------------------------
enable_testing()

add_executable(test_core
    tests/test_core.cpp
    src/fit/result.cpp
)
target_compile_features(test_core PRIVATE cxx_std_17)
target_compile_options(test_core PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)
target_include_directories(test_core PRIVATE ${SRC_DIR})

add_test(NAME core COMMAND test_core)

add_executable(test_config
    tests/test_config.cpp
    src/config/config.cpp
    src/core/fs.cpp
    src/core/log.cpp
)
target_compile_features(test_config PRIVATE cxx_std_17)
target_compile_options(test_config PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)
target_include_directories(test_config PRIVATE ${SRC_DIR})
target_include_directories(test_config SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(test_config PRIVATE
    ROOT::Core
    ROOT::Gpad
    nlohmann_json::nlohmann_json
)

add_test(NAME config COMMAND test_config)

add_executable(test_fit
    tests/test_fit.cpp
    src/fit/model.cpp
    src/fit/initial_parameters.cpp
    src/fit/result.cpp
    src/core/log.cpp
)
target_compile_features(test_fit PRIVATE cxx_std_17)
target_compile_options(test_fit PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)
target_include_directories(test_fit PRIVATE ${SRC_DIR})
target_include_directories(test_fit SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(test_fit PRIVATE
    ROOT::Core
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::MathCore
)

add_test(NAME fit COMMAND test_fit)

add_executable(test_projections
    tests/test_projections.cpp
    src/analysis/projections1d.cpp
    src/analysis/projections2d.cpp
    src/core/lcms.cpp
    src/core/fs.cpp
    src/core/log.cpp
    src/draw/draw.cpp
    src/io/input.cpp
    src/fit/model.cpp
    src/fit/initial_parameters.cpp
    src/fit/result.cpp
)
target_include_directories(test_projections PRIVATE ${SRC_DIR})
target_include_directories(test_projections SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(test_projections PRIVATE
    project_warnings
    ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::MathCore
)
add_test(NAME projections COMMAND test_projections)

add_executable(test_generated
    tests/test_generated.cpp
    src/io/input.cpp
    src/fit/model.cpp
    src/fit/initial_parameters.cpp
    src/fit/result.cpp
    src/core/log.cpp
)
target_include_directories(test_generated PRIVATE ${SRC_DIR})
target_include_directories(test_generated SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(test_generated PRIVATE
    project_warnings ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::MathCore
    nlohmann_json::nlohmann_json
)
add_test(NAME generate_fixture COMMAND ${CMAKE_COMMAND}
    -DGENERATOR=$<TARGET_FILE:generate_gaussian>
    -DOUTPUT_DIR=${CMAKE_BINARY_DIR}/gaussian_test_data
    -P ${CMAKE_SOURCE_DIR}/tests/generate_fixture.cmake
)
set_tests_properties(generate_fixture PROPERTIES FIXTURES_SETUP gaussian_data)
add_test(NAME generated_fit COMMAND test_generated ${CMAKE_BINARY_DIR}/gaussian_test_data/input.root)
add_test(NAME generated_pipeline COMMAND main ${CMAKE_BINARY_DIR}/gaussian_test_data/config.json)
set_tests_properties(generated_fit generated_pipeline PROPERTIES FIXTURES_REQUIRED gaussian_data)

# ----------------------------
# Install
# ----------------------------
install(TARGETS main DESTINATION bin)
