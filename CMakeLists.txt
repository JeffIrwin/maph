cmake_minimum_required(VERSION 3.7)

set(PROJECT maph)

if (LINUX)
	set(CMAKE_CXX_FLAGS "-Wall -Wextra")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
elseif (APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	# Windows
	# TODO:  set debug/release flags
endif()

project(${PROJECT})

set(MAPH_DIR src)

set(SUBMODULE_DIR submodules)

set(PNG_DIR ${SUBMODULE_DIR}/lodepng)
set(COLORMAPPER_DIR ${SUBMODULE_DIR}/colormapper)

add_subdirectory(${COLORMAPPER_DIR})

include_directories(${MAPH_DIR})

include_directories(
	${PNG_DIR}
	${COLORMAPPER_DIR}
	)

add_executable(${PROJECT} src/main.cpp ${PNG_DIR}/lodepng.cpp)

target_link_libraries(${PROJECT}
	colormapper
	)

