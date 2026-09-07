# Shim used when emp-tool is built in-tree via add_subdirectory: emp-ot's
# find_package(emp-tool 1.0 REQUIRED) lands here and only has to confirm the
# target already exists.
if(NOT TARGET emp-tool::emp-tool)
  message(FATAL_ERROR "emp-tool-shim: emp-tool::emp-tool target not defined (add_subdirectory(thirdparty/emp-tool) first)")
endif()
set(emp-tool_FOUND TRUE)
