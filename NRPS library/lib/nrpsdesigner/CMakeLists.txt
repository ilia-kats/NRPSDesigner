project(libnrpsdesigner)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/domains
)

set(SRCS
    monomer.cpp
    nrp.cpp

    abstractdomaintype.cpp
    domains/domaintypea.cpp
    domains/domaintypeac.cpp
    domains/domaintypet.cpp
    domains/domaintypee.cpp
    domains/domaintypete.cpp

    domain.cpp

    generalpathway.cpp
    pathway.cpp
    abstractdatabaseconnector.cpp
    mysqldatabaseconnector.cpp
)

add_library(nrpsdesigner ${SRCS})
generate_export_header(nrpsdesigner)
target_link_libraries(nrpsdesigner ${LIBXML2_LIBRARIES} ${MYSQLCPP_LIBRARIES})
