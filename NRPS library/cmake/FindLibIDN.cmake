# - Try to find LIBIDN
# Find LIBIDN headers, libraries
#
#  LIBIDN_FOUND               True if libidn got found
#  LIBIDN_INCLUDE_DIR        Location of libidn headers
#  LIBIDN_LIBRARIES           List of libaries to use libidn
#

find_path(LIBIDN_INCLUDE_DIR stringprep.h)
find_library(LIBIDN_LIBRARIES idn)

if(LIBIDN_LIBRARIES AND LIBIDN_INCLUDE_DIR)
    message(STATUS "Found libidn: ${LIBIDN_LIBRARIES}")
    set(LIBIDN_FOUND TRUE)
else(LIBIDN_LIBRARIES AND LIBIDN_INCLUDE_DIR)
    message(STATUS "Could NOT find libidn")
    set(LIBIDN_FOUND FALSE)
endif(LIBIDN_LIBRARIES AND LIBIDN_INCLUDE_DIR)

mark_as_advanced(LIBIDN_INCLUDE_DIR LIBIDN_LIBRARIES)
