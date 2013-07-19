project(libnrpsdesigner)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(SRCS
    monomer.cpp
    nrp.cpp
)

add_library(nrpsdesigner ${SRCS})
generate_export_header(nrpsdesigner)
target_link_libraries(nrpsdesigner ${LIBXML2_LIBRARIES})
