set(SRCS
    main.cpp
)

find_package(CURL REQUIRED)
include_directories(${CURL_INCLUDE_DIRS})


add_executable(nrpsdesigner-bin ${SRCS})
target_link_libraries(nrpsdesigner-bin nrpsdesigner ${CURL_LIBRARIES} ${MYSQLCPP_LIBRARIES})
set_target_properties(nrpsdesigner-bin PROPERTIES OUTPUT_NAME nrpsdesigner)
