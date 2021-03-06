include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(myproject_tests_SRCS
        #tests
        ${CMAKE_CURRENT_SOURCE_DIR}/test_hello_world.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testVerlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testLenardJonesDirectionSummation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testThermostat.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testNeighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testEmbeddedAtomPotential.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testHelperfunctions.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testMilstone.cpp
        #additional files
        ${CMAKE_SOURCE_DIR}/Sourcefiles/verlet.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/lenardJonesDirectionSummation.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/berendsenThermostat.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/helperfunctions.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/xyz.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/neighbors.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/embeddedAtomPotential.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/milestones.cpp
        ${CMAKE_SOURCE_DIR}/Sourcefiles/clusterGenerator.cpp
        )

# Tests
include_directories(yamd_tests ${CMAKE_SOURCE_DIR})
add_executable(myproject_tests ${myproject_SRCS} ${myproject_HDRS} ${myproject_tests_SRCS})
target_link_libraries(myproject_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(myproject_tests)

add_test(NAME myproject_tests COMMAND myproject_tests)
