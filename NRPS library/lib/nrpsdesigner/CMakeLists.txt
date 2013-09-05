project(libnrpsdesigner)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/domains
)

set(SRCS
    monomer.cpp

    domains/domaintypea.cpp
    domains/domaintypec.cpp
    domains/domaintypet.cpp
    domains/domaintypee.cpp
    domains/domaintypete.cpp

    global_enums.cpp

    domain.cpp

    origin.cpp
    product.cpp
    abstractdatabaseconnector.cpp
    mysqldatabaseconnector.cpp
    taxon.cpp
    nrps.cpp
    nrpsbuilder.cpp
)

if (UNIX)
    add_definitions("-DUSE_STRPTIME") # GCC does not yet implement std::get_time,
                                      # see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=54354
endif(UNIX)

add_library(nrpsdesigner ${SRCS})
generate_export_header(nrpsdesigner)
target_link_libraries(nrpsdesigner ${LIBXML2_LIBRARIES} ${CURL_LIBRARIES} ${MYSQL_LIBRARIES} ${MYSQLCPP_LIBRARIES})
