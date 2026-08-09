cmake_minimum_required(VERSION 3.16)
project(vurp_m_bpc VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../external/cmake_modules")
find_package(GUROBI REQUIRED)

add_executable(vurp_m_bpc
    src/main.cpp
    src/vurpm.cpp
)

target_include_directories(vurp_m_bpc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${GUROBI_INCLUDE_DIRS}
)

target_link_libraries(vurp_m_bpc PRIVATE ${GUROBI_LIBRARIES})

if (MSVC)
    target_compile_options(vurp_m_bpc PRIVATE /W4 /permissive-)
else ()
    target_compile_options(vurp_m_bpc PRIVATE -Wall -Wextra -Wpedantic)
endif ()
