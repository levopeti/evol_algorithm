cmake_minimum_required(VERSION 3.10)
project(Evolutionary_algorithm)

set(CMAKE_CXX_STANDARD 11)

set(APP_HEADERS ./Headers)
set(APP_SRC ./)
include_directories(${APP_HEADERS})
include_directories(${APP_SRC})



set(
        SOURCE_FILES
        main.cpp
        ${APP_HEADERS}/Gene.h
        ${APP_HEADERS}/Iteration.h
        ${APP_HEADERS}/Swarm.h
        ${APP_HEADERS}/Population.h
        ${APP_HEADERS}/Chromosome.h
        ${APP_HEADERS}/report.h
        ${APP_HEADERS}/python_test.h
        ${APP_SRC}/Gene.cpp
        ${APP_SRC}/Swarm.cpp
        ${APP_SRC}/Iteration.cpp
        ${APP_SRC}/Population.cpp
        ${APP_SRC}/Chromosome.cpp
        ${APP_SRC}/python_test.cpp
        Headers/Chromosome.h Headers/Population.h Headers/Grid_Search.h Headers/python_test.h)

add_executable(untitled ${SOURCE_FILES})
