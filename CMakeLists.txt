cmake_minimum_required (VERSION 3.21)
project (vole LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# EMP toolkit dependency: emp-tool and emp-ot 1.0 (main line), C++20.
#
# Default: build both from the git submodules in thirdparty/ (run
#   git submodule update --init
# after cloning). -DVOLE_USE_SYSTEM_EMP=ON uses installed packages instead
# (find_package(emp-tool 1.0) / find_package(emp-ot 1.0)).
# ---------------------------------------------------------------------------
option(VOLE_USE_SYSTEM_EMP "Use installed emp-tool/emp-ot 1.0 instead of the submodules" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(VOLE_USE_SYSTEM_EMP)
  find_package(emp-tool 1.0 REQUIRED)
  find_package(emp-ot 1.0 REQUIRED)
else()
  set(EMP_TOOL_DIR ${CMAKE_SOURCE_DIR}/thirdparty/emp-tool)
  set(EMP_OT_DIR   ${CMAKE_SOURCE_DIR}/thirdparty/emp-ot)
  if(NOT EXISTS ${EMP_TOOL_DIR}/CMakeLists.txt OR NOT EXISTS ${EMP_OT_DIR}/CMakeLists.txt)
    message(FATAL_ERROR "EMP submodules are missing. Run:  git submodule update --init")
  endif()
  set(EMP_TOOL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(EMP_TOOL_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
  set(EMP_TOOL_INSTALL OFF CACHE BOOL "" FORCE)
  set(EMP_OT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(EMP_OT_INSTALL OFF CACHE BOOL "" FORCE)
  add_subdirectory(${EMP_TOOL_DIR} ${CMAKE_BINARY_DIR}/emp-tool)
  # emp-ot's CMakeLists calls find_package(emp-tool 1.0 REQUIRED). Point it at
  # a shim config that just confirms the in-tree target exists.
  set(emp-tool_DIR ${CMAKE_SOURCE_DIR}/cmake/emp-tool-shim CACHE PATH "" FORCE)
  add_subdirectory(${EMP_OT_DIR} ${CMAKE_BINARY_DIR}/emp-ot)
endif()

# Flags for this repo's targets (the EMP targets carry their own via
# target_compile_features / options).
add_library(vole_flags INTERFACE)
target_compile_features(vole_flags INTERFACE cxx_std_20)
target_compile_options(vole_flags INTERFACE
  -O3 -Wall -Wno-ignored-attributes -Wno-unused-result -funroll-loops
  -march=native -maes -mrdseed -pthread)
target_include_directories(vole_flags INTERFACE ${CMAKE_SOURCE_DIR})
target_link_libraries(vole_flags INTERFACE emp-ot::emp-ot emp-tool::emp-tool pthread)
# Name used by test/ and bench/ CMakeLists.
set(EMP-OT_LIBRARIES vole_flags)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)   # correctness tests
ADD_SUBDIRECTORY(bench)  # runtime benchmarks
