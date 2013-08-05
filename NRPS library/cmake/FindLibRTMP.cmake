# - Try to find LIBRTMP
# Find LIBRTMP headers, libraries
#
#  LIBRTMP_FOUND               True if librtmp got found
#  LIBRTMP_INCLUDE_DIR        Location of librtmp headers
#  LIBRTMP_LIBRARIES           List of libaries to use librtmp
#

find_path(LIBRTMP_INCLUDE_DIR librtmp/rtmp.h)
find_library(LIBRTMP_LIBRARIES rtmp)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libRTMP  DEFAULT_MSG  LIBRTMP_LIBRARIES LIBRTMP_INCLUDE_DIR)

mark_as_advanced(LIBRTMP_INCLUDE_DIR LIBRTMP_LIBRARIES)
