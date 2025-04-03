cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(FemGauntlet LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(BUILD_SHARED_LIBS ON)

include_directories(.)

add_subdirectory(project)

set(EXEC_NAME solver)

# Set GMSH Version : example gmsh-4.12.2 -> GMSH_VERSION=gmsh-4.12, GMSH_MINOR=2
set(GMSH_VERSION gmsh-4.13)
set(GMSH_MINOR 1)

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

# GLFW
add_subdirectory(genmesh/glfw)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/glfw/include)

# GMSH
if(MSVC OR MSYS OR MINGW)
    set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Windows64-sdk/lib/gmsh.lib)
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Windows64-sdk/include)
    #Copy .dll into Debug
    set(GMSH_DLL ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Windows64-sdk/lib/${GMSH_VERSION}.dll)
    if(MSVC)
        file(COPY ${GMSH_DLL} DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE})
    else()
        file(COPY ${GMSH_DLL} DESTINATION ${CMAKE_BINARY_DIR})
    endif()
endif()

if(APPLE)
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
        set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSARM-sdk/lib/libgmsh.dylib)
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSARM-sdk/include)
    else()
        set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSX-sdk/lib/libgmsh.dylib)
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-MacOSX-sdk/include)
    endif()
endif()

if(UNIX AND NOT APPLE)
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Linux64-sdk/include)
    set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/genmesh/gmsh/${GMSH_VERSION}.${GMSH_MINOR}-Linux64-sdk/lib/libgmsh.so)
endif()

message("GMSH_INCLUDE_DIR : ${GMSH_INCLUDE_DIR}")
message("GMSH_LIBRARY : ${GMSH_LIBRARY}")
include_directories(${GMSH_INCLUDE_DIR})

# solver
add_executable(${EXEC_NAME}
    genmesh/src/glfem.c
    main.c
)


target_link_libraries(
    ${EXEC_NAME}
    PRIVATE
        femsolver
        glfw
        ${OPENGL_gl_LIBRARY}
        ${GMSH_LIBRARY}
)
target_compile_options(${EXEC_NAME} PUBLIC "-Wno-stringop-truncation") 
if(APPLE)
    target_compile_options(${EXEC_NAME} PRIVATE "-Wdeprecated-declarations")
    target_compile_definitions(${EXEC_NAME} PRIVATE -D GL_SILENCE_DEPRECATION)
endif()

if (MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXEC_NAME})
    target_compile_options(${EXEC_NAME} PRIVATE "/W2")
endif()