# - Try to find LIBGPGERROR
# Find LIBGPGERROR headers, libraries
#
#  LIBGPGERROR_FOUND               True if libgpgerror got found
#  LIBGPGERROR_INCLUDE_DIR        Location of libgpgerror headers
#  LIBGPGERROR_LIBRARIES           List of libaries to use libgpgerror
#

find_path(LIBGPGERROR_INCLUDE_DIR gpg-error.h)
find_library(LIBGPGERROR_LIBRARIES gpg-error)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libGPGError  DEFAULT_MSG  LIBGPGERROR_LIBRARIES LIBGPGERROR_INCLUDE_DIR)

mark_as_advanced(LIBGPGERROR_INCLUDE_DIR LIBGPGERROR_LIBRARIES)
