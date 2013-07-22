#ifndef NRPSDESIGNER_DOMAINTYPET_H
#define NRPSDESIGNER_DOMAINTYPET_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"
#include "abstractdomaintype.h"

namespace nrps
{
class NRPSDESIGNER_EXPORT DomainTypeT : public AbstractDomainType
{
public:
    DomainTypeT(DomainTPosition);
    virtual ~DomainTypeT();
    
    DomainTPosition position() const;
    
private:
    DomainTPosition m_position;
};
}

#endif
