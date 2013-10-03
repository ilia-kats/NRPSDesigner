# Search for libSBOL headers
INCLUDE(FindPackageMessage)

FIND_PATH(SBOL_INCLUDE_DIR sbol.h)
FIND_LIBRARY(SBOL_LIBRARY NAMES sbol)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SBOL DEFAULT_MSG SBOL_LIBRARY SBOL_INCLUDE_DIR)

MARK_AS_ADVANCED(
    SBOL_INCLUDE_DIR
    SBOL_LIBRARY
    SBOL_FOUND
)

