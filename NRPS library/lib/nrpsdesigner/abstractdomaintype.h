#ifndef NRPSDESIGNER_ABSTRACTDOMAINTYPE_H
#define NRPSDESIGNER_ABSTRACTDOMAINTYPE_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"

#include <cstddef>
#include <functional>

namespace nrps
{
class NRPSDESIGNER_EXPORT AbstractDomainType
{
public:
    virtual ~AbstractDomainType();

    DomainType type() const;
    virtual std::size_t hash() const;
    virtual bool full() const;

protected:
    AbstractDomainType(DomainType type);

    DomainType m_type;
};
}

namespace std
{
template<>
struct NRPSDESIGNER_EXPORT hash<nrps::AbstractDomainType>
{
public:
    std::size_t operator()(const nrps::AbstractDomainType &d) const;
};
}

#endif
