#ifndef NRPSDESIGNER_DOMAINTYPETE_H
#define NRPSDESIGNER_DOMAINTYPETE_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"

namespace nrps
{
class NRPSDESIGNER_EXPORT DomainTypeTe : public AbstractDomainType
{
public:
    DomainTypeTe(bool);
    virtual ~DomainTypeTe();

    bool circularizing() const;

private:
    bool m_circularizing;
};
}

#endif
