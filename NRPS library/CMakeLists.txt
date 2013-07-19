project(NRPSDesigner)
cmake_minimum_required(VERSION 2.8)

find_package(LibXml2 REQUIRED)

set(CMAKE_CXX_FLAGS
"${CMAKE_CXX_FLAGS} -std=c++11"
)
include(GenerateExportHeader)
add_compiler_export_flags()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/nrpsdesigner
    ${CMAKE_CURRENT_BINARY_DIR}/lib/nrpsdesigner
    ${LIBXML2_INCLUDE_DIR}
)

add_subdirectory(lib)
add_subdirectory(ui)
