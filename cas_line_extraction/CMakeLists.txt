project(cas_line_extraction)

# Add the include folder to the include path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add default source files
set(CAS_SRCS src/cas_line_extraction.cpp)

# Add header files
set(CAS_HEADERS ${PROJECT_SOURCE_DIR}/include/GPSNetworkInterface.h ${PROJECT_SOURCE_DIR}/include/MockGPS.h)

## Build the library
#set(CAS_LINK_LIBS  eigen)

# Compile the libmdc2250 Library
add_library(cas_line_extraction ${CAS_SRCS})
#target_link_libraries(cas_line_extraction )

## Build Examples
# If specified
IF(GPSINTERFACE_BUILD_EXAMPLES)
  # Compile the Example program
  #add_executable(gps_example examples/atrv_example.cpp)
  # Link the Test program to the atrv library
  #target_link_libraries(gps_example gps_interfaces ${GPS_LINK_LIBS})
ENDIF(GPSINTERFACE_BUILD_EXAMPLES)
