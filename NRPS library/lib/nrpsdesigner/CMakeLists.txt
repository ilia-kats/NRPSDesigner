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
    networkoptions.cpp
    taxon.cpp
    taxonbuilder.cpp
    nrps.cpp
    nrpsbuilder.cpp
)
set(HEADERS
    monomer.h
    domains/domaintypea.h
    domains/domaintypec.h
    domains/domaintypet.h
    domains/domaintypee.h
    domains/domaintypete.h
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
    exceptions.h
    ${CMAKE_CURRENT_BINARY_DIR}/nrpsdesigner_export.h
)

if (UNIX)
    add_definitions("-DUSE_STRPTIME") # GCC does not yet implement std::get_time,
                                      # see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=54354
endif(UNIX)

add_library(nrpsdesigner ${SRCS})
generate_export_header(nrpsdesigner)
target_link_libraries(nrpsdesigner ${LIBXML2_LIBRARIES} ${CURL_LIBRARIES} ${MYSQL_LIBRARIES} ${MYSQLCPP_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
set_target_properties(nrpsdesigner PROPERTIES PUBLIC_HEADER "${HEADERS}")
install(TARGETS nrpsdesigner LIBRARY DESTINATION ${LIB_SUFFIX} RUNTIME DESTINATION bin ARCHIVE DESTINATION lib PUBLIC_HEADER DESTINATION include)
