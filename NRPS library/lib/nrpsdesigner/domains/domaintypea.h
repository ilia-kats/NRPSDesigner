#ifndef NRPSDESIGNER_DOMAINTYPEA_H
#define NRPSDESIGNER_DOMAINTYPEA_H

#include "nrpsdesigner_export.h"
#include "domain.h"

namespace nrps
{
class DomainTypeA : public Domain
{
public:
    DomainTypeA(uint32_t);
    DomainTypeA(uint32_t, uint32_t);
    virtual ~DomainTypeA();

    uint32_t substrate() const;
    void setSubstrate(uint32_t);

protected:
    DomainTypeA(DomainType, uint32_t);
    DomainTypeA(DomainType, uint32_t, uint32_t);
    virtual void writeXml(xmlTextWriterPtr) const;

private:
    uint32_t m_substrate;
};
}

#endif
