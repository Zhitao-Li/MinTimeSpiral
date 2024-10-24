cmake_minimum_required(VERSION 3.20)

project(SDC3)


find_package(Matlab COMPONENTS MX_LIBRARY REQUIRED)

set(MEX_SRC
    ${PROJECT_SOURCE_DIR}/Sources/sdc3_MAT.c
)

matlab_add_mex(
    NAME sdc3_MAT
    SRC ${MEX_SRC}
)

set_target_properties(sdc3_MAT
    PROPERTIES
    COMPILE_FLAGS "-fvisibility=default"
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
