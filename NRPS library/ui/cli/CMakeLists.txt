set(SRCS
    main.cpp
)

add_executable(nrpsdesigner-bin ${SRCS})
target_link_libraries(nrpsdesigner-bin nrpsdesigner)
set_target_properties(nrpsdesigner-bin PROPERTIES OUTPUT_NAME nrpsdesigner)
install(TARGETS nrpsdesigner-bin RUNTIME DESTINATION bin)
