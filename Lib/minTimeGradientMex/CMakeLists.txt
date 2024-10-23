cmake_minimum_required(VERSION 3.20)

project(minTimeGradient)


find_package(Matlab COMPONENTS MX_LIBRARY REQUIRED)

set(MEX_SRC ${PROJECT_SOURCE_DIR}/Sources/minTimeGradient.c)

matlab_add_mex(
    NAME minTimeGradient
    SRC ${MEX_SRC}
)

set_target_properties(minTimeGradient
    PROPERTIES
    COMPILE_FLAGS "-fvisibility=default"
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)

add_custom_command(TARGET minTimeGradient POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/Sources/vdSpiralDesign.m" "${PROJECT_SOURCE_DIR}/bin"
)
