# - Try to find LIBTASN
# Find LIBTASN headers, libraries
#
#  LIBTASN_FOUND               True if libtasn got found
#  LIBTASN_INCLUDE_DIR        Location of libtasn headers
#  LIBTASN_LIBRARIES           List of libaries to use libtasn
#

find_path(LIBTASN_INCLUDE_DIR libtasn1.h)
find_library(LIBTASN_LIBRARIES tasn1)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libtasn1 DEFAULT_MSG  LIBTASN_LIBRARIES LIBTASN_INCLUDE_DIR)

mark_as_advanced(LIBTASN_INCLUDE_DIR LIBTASN_LIBRARIES)
