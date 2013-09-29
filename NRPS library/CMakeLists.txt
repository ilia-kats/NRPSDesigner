project(NRPSDesigner)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/kde
)

find_package(LibXml2 REQUIRED)
find_package(CURL REQUIRED)
find_package(MySQL REQUIRED)
find_package(MySQLcpp REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

set(CMAKE_CXX_FLAGS
"${CMAKE_CXX_FLAGS} -std=c++11"
)
include(GenerateExportHeader)
add_compiler_export_flags()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/nrpsdesigner
    ${CMAKE_CURRENT_BINARY_DIR}/lib/nrpsdesigner
    ${LIBXML2_INCLUDE_DIR}
    ${MYSQL_INCLUDE_DIR}
    ${MYSQLCPP_INCLUDE_DIR}
    ${CURL_INCLUDE_DIRS}
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   add_definitions("-D_CXX_CLANG")
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

if (CMAKE_LIBRARY_ARCHITECTURE)
    set (LIB_INSTALL_DIR lib${LIB_SUFFIX}/${CMAKE_LIBRARY_ARCHITECTURE})
else (CMAKE_LIBRARY_ARCHITECTURE)
    set (LIB_INSTALL_DIR lib${LIB_SUFFIX})
endif (CMAKE_LIBRARY_ARCHITECTURE)

add_subdirectory(lib)
add_subdirectory(ui)
