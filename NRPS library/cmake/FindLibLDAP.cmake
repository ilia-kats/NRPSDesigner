# - Try to find LIBLDAP
# Find LIBLDAP headers, libraries
#
#  LIBLDAP_FOUND               True if libldap got found
#  LIBLDAP_INCLUDE_DIR        Location of libldap headers
#  LIBLDAP_LIBRARIES           List of libraries to use libldap
#  LIBLDAP_LIBRARY            The libldap library
#  LIBLDAPR_LIBRARY           The libldap_r library
#  LIBLBER_LIBRARY            The liblber library
#

find_path(LIBLDAP_INCLUDE_DIR ldap.h)
find_library(LIBLDAP_LIBRARY ldap)
find_library(LIBLDAPR_LIBRARY ldap_r)
find_library(LIBLBER_LIBRARY lber)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libLDAP DEFAULT_MSG  LIBLDAP_LIBRARY LIBLDAPR_LIBRARY LIBLBER_LIBRARY LIBLDAP_INCLUDE_DIR)
if (LIBLDAP_FOUND)
    set(LIBLDAP_LIBRARIES ${LIBLDAP_LIBRARY} ${LIBLDAPR_LIBRARY} ${LIBLBER_LIBRARY})
endif (LIBLDAP_FOUND)

mark_as_advanced(LIBLDAP_INCLUDE_DIR LIBLDAP_LIBRARIES LIBLDAP_LIBRARY LIBLDAPR_LIBRARY LIBLBER_LIBRARY)
