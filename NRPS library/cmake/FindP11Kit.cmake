# - Try to find LIBP11KIT
# Find P11KIT headers, libraries
#
#  P11KIT_FOUND               True if libp11kit got found
#  P11KIT_INCLUDE_DIR        Location of libp11kit headers
#  P11KIT_LIBRARIES           List of libaries to use libp11kit
#

find_path(P11KIT_INCLUDE_DIR p11-kit/p11-kit.h PATH_SUFFIXES p11-kit-1)
find_library(P11KIT_LIBRARIES p11-kit)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(p11kit DEFAULT_MSG  P11KIT_LIBRARIES P11KIT_INCLUDE_DIR)

mark_as_advanced(P11KIT_INCLUDE_DIR P11KIT_LIBRARIES)
