cmake_minimum_required(VERSION 3.16)

project(cf_maker_new LANGUAGES CXX)

# ===============================
# Общие настройки
# ===============================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ===============================
# ROOT
# ===============================
find_package(ROOT REQUIRED COMPONENTS
    Core
    RIO
    Hist
    Graf
    Graf3d
    Gpad
    Tree
    TreePlayer
    Physics
    MathCore
)

include(${ROOT_USE_FILE})

# ===============================
# Пути
# ===============================
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(SRC_DIR     ${CMAKE_SOURCE_DIR}/src)

include_directories(${INCLUDE_DIR})

# ===============================
# Исходники
# ===============================

# main.cpp отдельно
set(MAIN_SOURCE ${CMAKE_SOURCE_DIR}/main.cpp)

# Все .cpp из src/
file(GLOB SRC_SOURCES CONFIGURE_DEPENDS
    ${SRC_DIR}/*.cpp
)

# ===============================
# Исполняемый файл
# ===============================
add_executable(main
    ${MAIN_SOURCE}
    ${SRC_SOURCES}
)

target_include_directories(main PRIVATE ${INCLUDE_DIR})

target_compile_options(main PRIVATE
    -Wall
    -O2
    -DROOT_STRICT_OWNERSHIP
)

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
)

install(TARGETS main DESTINATION bin)
