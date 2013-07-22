project(libnrpsdesigner)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/domaintypes
)

set(SRCS
    monomer.cpp
    nrp.cpp

    abstractdomaintype.cpp
    domaintypes/domaintypea.cpp
    domaintypes/domaintypeac.cpp
    domaintypes/domaintypet.cpp
    domaintypes/domaintypee.cpp
    domaintypes/domaintypete.cpp

    generalpathway.cpp
)

add_library(nrpsdesigner ${SRCS})
generate_export_header(nrpsdesigner)
target_link_libraries(nrpsdesigner ${LIBXML2_LIBRARIES})
