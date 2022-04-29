cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

project(Practica3)

add_executable(Practica3
    src/Fluid2.h
    src/Fluid2Exercise.cpp
    src/FluidVisualizer2.cpp
    src/FluidVisualizer2.h
    src/Grid2.h
    src/main.cpp
    src/Particles2.h
    src/Scene.cpp
    src/Scene.h
    src/Containers/Array2.h
    src/Containers/Index2.cpp
    src/Containers/Index2.h
    src/Math/AABox2.h
    src/Math/Vector2.cpp
    src/Math/Vector2.h
    src/Math/Vector3.cpp
    src/Math/Vector3.h
    src/Numeric/SparseMatrix.h
    src/Numeric/PCGSolver.h
    src/Std/CMath.h
    src/Std/CStdInt.h
    )

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

target_include_directories(Practica3
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/src
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    )

target_link_libraries(Practica3
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARY}
    )

if (WIN32)
    add_custom_command(
        TARGET Practica3
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} ARGS -E copy_if_different
        "${PROJECT_SOURCE_DIR}/External/freeglut-3.2.1/build-vs2017/bin/$<CONFIG>/freeglut$<$<CONFIG:Debug>:d>.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
        )
endif ()
