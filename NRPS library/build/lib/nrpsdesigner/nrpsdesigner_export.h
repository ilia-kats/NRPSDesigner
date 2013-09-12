
#ifndef NRPSDESIGNER_EXPORT_H
#define NRPSDESIGNER_EXPORT_H

#ifdef NRPSDESIGNER_STATIC_DEFINE
#  define NRPSDESIGNER_EXPORT
#  define NRPSDESIGNER_NO_EXPORT
#else
#  ifndef NRPSDESIGNER_EXPORT
#    ifdef nrpsdesigner_EXPORTS
        /* We are building this library */
#      define NRPSDESIGNER_EXPORT 
#    else
        /* We are using this library */
#      define NRPSDESIGNER_EXPORT 
#    endif
#  endif

#  ifndef NRPSDESIGNER_NO_EXPORT
#    define NRPSDESIGNER_NO_EXPORT 
#  endif
#endif

#ifndef NRPSDESIGNER_DEPRECATED
#  define NRPSDESIGNER_DEPRECATED __attribute__ ((__deprecated__))
#  define NRPSDESIGNER_DEPRECATED_EXPORT NRPSDESIGNER_EXPORT __attribute__ ((__deprecated__))
#  define NRPSDESIGNER_DEPRECATED_NO_EXPORT NRPSDESIGNER_NO_EXPORT __attribute__ ((__deprecated__))
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define NRPSDESIGNER_NO_DEPRECATED
#endif

#endif
