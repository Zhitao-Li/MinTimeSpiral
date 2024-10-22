cmake_minimum_required(VERSION 3.20)

project(SOCP)


find_package(Matlab COMPONENTS MX_LIBRARY REQUIRED)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(MEX_SRC ${PROJECT_SOURCE_DIR}/Sources/socp_mex.c)

matlab_add_mex(
    NAME socp_mex
    SRC ${MEX_SRC}
    LINK_TO ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

set_target_properties(socp_mex
    PROPERTIES
    COMPILE_FLAGS "-fvisibility=default"
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)

add_custom_command(TARGET socp_mex POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/Sources/socp.m" "${PROJECT_SOURCE_DIR}/bin"
)
