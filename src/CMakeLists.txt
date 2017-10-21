########################################
# Set up how to compile the source files
########################################

# Define the executable name
SET(2DEULEREXE 2deuler.exe)
SET(MESHGENEXE meshgen.exe)

# Add the source files
SET(2DEULER_Sources ${SRC}/main.f90
                    ${SRC}/readinput.f90
                    ${SRC}/mesh.f90
                    ${SRC}/solution.f90
                    ${SRC}/hllc.f90
                    ${SRC}/movingmesh.f90
                    ${SRC}/library.f90
)

SET(MESHGEN_Sources ${SRC}/meshgen.f90)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${2DEULEREXE} ${2DEULER_Sources})
ADD_EXECUTABLE(${MESHGENEXE} ${MESHGEN_Sources})

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${2DEULEREXE} ${SOME_LIBRARY})

#####################################
# Tell how to install this executable
#####################################

set_target_properties( ${2DEULEREXE}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

set_target_properties( ${MESHGENEXE}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
