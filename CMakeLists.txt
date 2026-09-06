# 3.5 is the lowest floor accepted by CMake 4.x (it dropped <3.5 compatibility).
cmake_minimum_required (VERSION 3.5)
project (vole)

# ---------------------------------------------------------------------------
# EMP toolkit dependency: emp-tool + emp-ot at their 0.3.0 releases.
#
# Default: build from the git submodules in thirdparty/ (pinned commits; run
#   git submodule update --init
# after cloning). emp-ot is header-only; emp-tool contributes one library.
#
# -DVOLE_USE_SYSTEM_EMP=ON uses an installed EMP instead (must be 0.3.0: the
# 1.0 rewrite removed BaseCot/ferret and does not build this code).
# ---------------------------------------------------------------------------
option(VOLE_USE_SYSTEM_EMP "Use an installed EMP 0.3.0 instead of the submodules" OFF)

if(VOLE_USE_SYSTEM_EMP)
  find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
  include(${CMAKE_FOLDER}/cmake/emp-base.cmake)
  # emp-ot's config pulls in emp-tool; EMP-OT_LIBRARIES is emp-tool's library.
  find_package(emp-ot REQUIRED)
  include_directories(${EMP-OT_INCLUDE_DIRS})
else()
  set(EMP_TOOL_DIR ${CMAKE_SOURCE_DIR}/thirdparty/emp-tool)
  set(EMP_OT_DIR   ${CMAKE_SOURCE_DIR}/thirdparty/emp-ot)
  if(NOT EXISTS ${EMP_TOOL_DIR}/emp-tool/emp-tool.h OR
     NOT EXISTS ${EMP_OT_DIR}/emp-ot/emp-ot.h)
    message(FATAL_ERROR
      "EMP submodules are missing. Run:  git submodule update --init")
  endif()
  # emp-base.cmake: compiler flags (-march=native -maes -std=c++11 -O3),
  # OpenSSL root on macOS, bin/ output dir, optional THREADING/rdseed flags.
  set(CMAKE_FOLDER ${EMP_TOOL_DIR})
  include(${EMP_TOOL_DIR}/cmake/emp-base.cmake)
  find_package(OpenSSL REQUIRED)
  include_directories(${OPENSSL_INCLUDE_DIR} ${EMP_TOOL_DIR} ${EMP_OT_DIR})
  add_definitions(-DEMP_CIRCUIT_PATH=${EMP_TOOL_DIR}/emp-tool/circuits/files/)
  add_library(emp-tool STATIC ${EMP_TOOL_DIR}/emp-tool/emp-tool.cpp)
  target_link_libraries(emp-tool ${OPENSSL_LIBRARIES})
  set(EMP-OT_LIBRARIES emp-tool)
endif()

include_directories(${CMAKE_SOURCE_DIR})  # so "vole/..." headers resolve

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)   # correctness tests
ADD_SUBDIRECTORY(bench)  # runtime benchmarks
