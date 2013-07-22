#ifndef NRPSDESIGNER_DOMAINTYPEAC_H
#define NRPSDESIGNER_DOMAINTYPEAC_H

#include "nrpsdesigner_export.h"
#include "domaintypea.h"
#include "global_enums.h"

namespace nrps
{
class NRPSDESIGNER_EXPORT DomainTypeAC : public DomainTypeA
{
public:
    DomainTypeAC(uint32_t, Configuration);
    virtual ~DomainTypeAC();
    
    Configuration chirality() const;
    
private:
    Configuration m_chirality;
};
}

#endif

