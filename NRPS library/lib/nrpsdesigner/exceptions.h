#ifndef NRPSDESIGNER_EXCEPTIONS_H
#define NRPSDESIGNER_EXCEPTIONS_H

#include "nrpsdesigner_export.h"

#include <exception>
#include <stdexcept>

namespace nrps
{
class NRPSDESIGNER_EXPORT NetworkError : public std::runtime_error
{using std::runtime_error::runtime_error;};

class NRPSDESIGNER_EXPORT NCBITaxonomyError : public NetworkError
{using NetworkError::NetworkError;};

class NRPSDESIGNER_EXPORT TaxonomyDumpError : public std::runtime_error
{using std::runtime_error::runtime_error;};

class NRPSDESIGNER_EXPORT DatabaseError : public std::runtime_error
{using std::runtime_error::runtime_error;};

}

#endif
