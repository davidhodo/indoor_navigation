project(cas_line_extraction)

# Add the include folder to the include path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add default source files
set(CAS_SRCS src/cas_line_extraction.cpp)

# Add header files
set(CAS_HEADERS ${PROJECT_SOURCE_DIR}/include/cas_line_extraction/cas_line_extraction.h)

## Build the library
#set(CAS_LINK_LIB  eigen)

# Compile the libmdc2250 Library
add_library(cas_line_extraction ${CAS_SRCS} ${CAS_HEADERS})
#target_link_libraries(cas_line_extraction )

## Build Examples
option(CAS_BUILD_EXAMPLES "Build all of the CAS examples." OFF)
IF(CAS_BUILD_EXAMPLES)
  # Compile the Example program
  add_executable(cas_line_extraction_example examples/example.cpp)
  # Link the Test program to the cas_line_extraction library
  target_link_libraries(cas_line_extraction_example cas_line_extraction ${CAS_LINK_LIB})
ENDIF(CAS_BUILD_EXAMPLES)

## Build tests
option(CAS_BUILD_TESTS "Build all of the CAS tests." OFF)
IF(CAS_BUILD_TESTS)
    # Find Google Test
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})

    # Compile the tests program
    add_executable(cas_line_extraction_tests tests/line_extraction_tests.cpp)
    # Link the Test program to the cas_line_extraction library
    target_link_libraries(cas_line_extraction_tests cas_line_extraction ${GTEST_BOTH_LIBRARIES} pthread)
ENDIF(CAS_BUILD_TESTS)
