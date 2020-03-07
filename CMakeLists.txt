cmake_minimum_required(VERSION 3.0)
project(RayTracing VERSION 1.0 LANGUAGES CXX)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 11)

if (NOT UNIX)
	message(FATAL_ERROR "Not Unix!")
endif()

include_directories(include)

file(GLOB SOURCES src/*.cpp)

add_executable(rt ${SOURCES})

