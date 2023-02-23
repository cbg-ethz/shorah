# file adapted from: https://gitlab.inria.fr/cado-nfs/cado-nfs/
# original License: Gnu Lesser General Public License (LGPL) version 2.1 (or any later version).

# popcnt
message(STATUS "Testing whether popcnt code can be used")

try_run(popcnt_runs popcnt_compiles
    ${PROJECT_BINARY_DIR}/cmake
    ${PROJECT_SOURCE_DIR}/cmake/popcnt.c
    )
if(popcnt_compiles)
    if (popcnt_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether popcnt code can be used -- No")
        set (HAVE_POPCNT 0)
    else()
        message(STATUS "Testing whether popcnt code can be used -- Yes")
        set (HAVE_POPCNT 1)
    endif()
else()
    try_run(popcnt_runs popcnt_compiles
        ${PROJECT_BINARY_DIR}/cmake
        ${PROJECT_SOURCE_DIR}/cmake/popcnt.c
        COMPILE_DEFINITIONS -mpopcnt
        )
    if(popcnt_compiles)
        if (popcnt_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether popcnt code can be used -- No")
            set (HAVE_POPCNT 0)
        else()
            message(STATUS "Testing whether popcnt code can be used -- Yes, with -mpopcnt")
            set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mpopcnt")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
            set (HAVE_POPCNT 1)
        endif()
    else()
        message(STATUS "Testing whether popcnt code can be used -- No")
        set (HAVE_POPCNT 0)
    endif()
endif()
