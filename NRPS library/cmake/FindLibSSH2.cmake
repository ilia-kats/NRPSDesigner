# - Try to find LIBSSH2
# Find LIBSSH2 headers, libraries
#
#  LIBSSH2_FOUND               True if libssh2 got found
#  LIBSSH2_INCLUDE_DIR        Location of libssh2 headers
#  LIBSSH2_LIBRARIES           List of libaries to use libssh2
#

find_path(LIBSSH2_INCLUDE_DIR libssh2.h)
find_library(LIBSSH2_LIBRARIES ssh2)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libSSH2  DEFAULT_MSG  LIBSSH2_LIBRARIES LIBSSH2_INCLUDE_DIR)

mark_as_advanced(LIBSSH2_INCLUDE_DIR LIBSSH2_LIBRARIES)
