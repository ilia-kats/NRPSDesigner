#ifndef NRPSDESIGNER_DOMAINTYPEAOXA_H
#define NRPSDESIGNER_DOMAINTYPEAOXA_H

#include "nrpsdesigner_export.h"
#include "domaintypea.h"

namespace nrps
{
class DomainTypeAOxA : public DomainTypeA
{
public:
    DomainTypeAOxA(uint32_t);
    DomainTypeAOxA(uint32_t, uint32_t);
    virtual ~DomainTypeAOxA();
};
}

#endif
