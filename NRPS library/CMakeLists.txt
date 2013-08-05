project(NRPSDesigner)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/kde
)

option(BUILD_STATIC "Link statically." OFF)
if (BUILD_STATIC)

    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -pthread"
    )
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++"
        )
    endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    find_package(Libdl REQUIRED) # dependency of mysqlcpp

    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

    find_package(MySQL REQUIRED) # dependency of mysqlcpp
    find_package(LibIDN REQUIRED) # dependency of curl
    find_package(GnuTLS REQUIRED) # dependency of curl
    find_package(LibSSH2 REQUIRED) # dependency of curl
    find_package(LibLDAP REQUIRED) # dependency of curl
    find_package(LibGCrypt REQUIRED) # dependency of curl
    find_package(GSSAPI REQUIRED) # dependency of curl
    find_package(LibRTMP REQUIRED) # dependency of curl
    find_package(LibGPGError REQUIRED) # dependency of gcrypt
    find_package(Libtasn1 REQUIRED) # dependency of gnutls
    find_package(P11Kit REQUIRED) # dependency of gnutls
    find_package(ZLIB REQUIRED) # dependency of libxml2
    find_package(LibLZMA REQUIRED) # dependency of libxml2
    find_package(Sasl2 REQUIRED) # dependency of libldap
    find_package(BerkeleyDB REQUIRED) # dependency of libsasl2
    find_package(OpenSSL REQUIRED) # dependency of libsasl2
    find_package(Sqlite REQUIRED) # dependency of libsasl2
    find_package(PostgreSQL REQUIRED) # dependency of libsasl2
endif(BUILD_STATIC)

find_package(LibXml2 REQUIRED)
find_package(MySQLcpp REQUIRED)

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
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   add_definitions("-D_CXX_CLANG")
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

add_subdirectory(lib)
add_subdirectory(ui)
