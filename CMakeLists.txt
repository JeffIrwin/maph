cmake_minimum_required(VERSION 3.7)

set(PROJECT maph)

if (LINUX)
	set(CMAKE_CXX_FLAGS "-Wall -Wextra")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
else()
	# Windows
	# TODO:  set debug/release flags
endif()

project(${PROJECT})

set(MAPH_DIR src)

set(SUBMODULE_DIR submodules)

set(PNG_DIR ${SUBMODULE_DIR}/lodepng)
set(JSON_DIR ${SUBMODULE_DIR}/json/single_include/nlohmann)
set(XML_DIR ${SUBMODULE_DIR}/pugixml/src)

include_directories(${MAPH_DIR})

include_directories(${PNG_DIR})
include_directories(${JSON_DIR})
include_directories(${XML_DIR})

add_executable(${PROJECT} src/main.cpp ${PNG_DIR}/lodepng.cpp ${XML_DIR}/pugixml.cpp)

