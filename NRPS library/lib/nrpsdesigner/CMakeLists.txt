project(libnrpsdesigner)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(SRCS
    monomer.cpp
    nrp.cpp
)

add_library(nrpsdesign ${SRCS})
generate_export_header(nrpsdesign)
target_link_libraries(nrpsdesign ${LIBXML2_LIBRARIES})
