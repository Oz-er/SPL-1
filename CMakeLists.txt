cmake_minimum_required(VERSION 3.10)

project(DomainAdaption)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(DomainAdaption
    Main.cpp
    TCA.cpp
    BDA.cpp
    CORAL.cpp
)