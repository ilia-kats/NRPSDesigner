# - Try to find LIBGCRYPT
# Find LIBGCRYPT headers, libraries
#
#  LIBGCRYPT_FOUND               True if libgcrypt got found
#  LIBGCRYPT_INCLUDE_DIR        Location of libgcrypt headers
#  LIBGCRYPT_LIBRARIES           List of libaries to use libgcrypt
#

find_path(LIBGCRYPT_INCLUDE_DIR gcrypt.h)
find_library(LIBGCRYPT_LIBRARIES gcrypt)

if(LIBGCRYPT_LIBRARIES AND LIBGCRYPT_INCLUDE_DIR)
    message(STATUS "Found libgcrypt: ${LIBGCRYPT_LIBRARIES}")
    set(LIBGCRYPT_FOUND TRUE)
else(LIBGCRYPT_LIBRARIES AND LIBGCRYPT_INCLUDE_DIR)
    message(STATUS "Could NOT find libgcrypt")
    set(LIBGCRYPT_FOUND FALSE)
endif(LIBGCRYPT_LIBRARIES AND LIBGCRYPT_INCLUDE_DIR)

mark_as_advanced(LIBGCRYPT_INCLUDE_DIR LIBGCRYPT_LIBRARIES)
