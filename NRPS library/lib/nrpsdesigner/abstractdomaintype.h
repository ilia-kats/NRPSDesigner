#ifndef NRPSDESIGNER_ABSTRACTDOMAINTYPE_H
#define NRPSDESIGNER_ABSTRACTDOMAINTYPE_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"

namespace nrps
{
class NRPSDESIGNER_EXPORT AbstractDomainType
{
public:
    virtual ~AbstractDomainType();

    DomainType type() const;

protected:
    AbstractDomainType(DomainType type);

    DomainType m_type;
};
}

#endif
