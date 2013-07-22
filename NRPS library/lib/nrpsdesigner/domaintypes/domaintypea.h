#ifndef NRPSDESIGNER_DOMAINTYPEA_H
#define NRPSDESIGNER_DOMAINTYPEA_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"

#include <cstdint>

namespace nrps
{
class NRPSDESIGNER_EXPORT DomainTypeA : public AbstractDomainType
{
public:
    DomainTypeA(uint32_t);
    virtual ~DomainTypeA();
    
    uint32_t substrate() const;
    
protected:
    DomainTypeA(uint32_t, DomainType);
    
private:
    uint32_t m_substrate;
};
}

#endif
