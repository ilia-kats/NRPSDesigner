# - Try to find LIBDL
# Find LIBDL headers, libraries
#
#  LIBDL_FOUND               True if libdl got found
#  LIBDL_INCLUDE_DIR        Location of libdl headers
#  LIBDL_LIBRARIES           List of libaries to use libdl
#

find_path(LIBDL_INCLUDE_DIR dlfcn.h)
find_library(LIBDL_LIBRARIES NAMES dl ltdl)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libdl  DEFAULT_MSG  LIBDL_LIBRARIES LIBDL_INCLUDE_DIR)

mark_as_advanced(LIBDL_INCLUDE_DIR LIBDL_LIBRARIES)
