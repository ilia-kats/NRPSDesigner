set(SRCS
    main.cpp
)

add_executable(nrpsdesigner ${SRCS})
target_link_libraries(nrpsdesigner nrpsdesign)
