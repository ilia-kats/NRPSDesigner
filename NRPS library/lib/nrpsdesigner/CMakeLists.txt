project(libnrpsdesigner)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/domains
    ${CMAKE_CURRENT_BINARY_DIR}
)

configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

set(SRCS
    domains/domaintypea.cpp
    domains/domaintypec.cpp
    domains/domaintypet.cpp
    domains/domaintypee.cpp
    domains/domaintypete.cpp
    domains/domaintypeaoxa.cpp
    domain.cpp
    monomer.cpp
    origin.cpp
    product.cpp
    abstractdatabaseconnector.cpp
    mysqldatabaseconnector.cpp
    networkoptions.cpp
    taxon.cpp
    taxonbuilder.cpp
    nrps.cpp
    nrpsbuilder.cpp
    nrpslibrary.cpp

    global_enums.cpp

    globals_internal.cpp
)
set(HEADERS
    monomer.h
    domains/domaintypea.h
    domains/domaintypec.h
    domains/domaintypet.h
    domains/domaintypee.h
    domains/domaintypete.h
    domains/domaintypeaoxa.h
    global_enums.h
    domain.h
    origin.h
    product.h
    abstractdatabaseconnector.h
    networkoptions.h
    taxon.h
    taxonbuilder.h
    nrps.h
    nrpsbuilder.h
    nrpslibrary.h
    exceptions.h
    ${CMAKE_CURRENT_BINARY_DIR}/nrpsdesigner_export.h
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

if (UNIX)
    add_definitions("-DUSE_STRPTIME") # GCC does not yet implement std::get_time,
                                      # see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=54354
endif(UNIX)

add_library(nrpsdesigner ${SRCS})
generate_export_header(nrpsdesigner)
target_link_libraries(nrpsdesigner ${LIBXML2_LIBRARIES} ${CURL_LIBRARIES} ${MYSQL_LIBRARIES} ${MYSQLCPP_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (WITH_SBOL)
    target_link_libraries(nrpsdesigner ${SBOL_LIBRARY})
endif (WITH_SBOL)

set_target_properties(nrpsdesigner PROPERTIES PUBLIC_HEADER "${HEADERS}")
install(TARGETS nrpsdesigner LIBRARY DESTINATION ${LIB_SUFFIX} RUNTIME DESTINATION bin ARCHIVE DESTINATION lib PUBLIC_HEADER DESTINATION include)
