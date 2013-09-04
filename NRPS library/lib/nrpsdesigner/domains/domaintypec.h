#ifndef NRPSDESIGNER_DOMAINTYPEC_H
#define NRPSDESIGNER_DOMAINTYPEC_H

#include "nrpsdesigner_export.h"
#include "domain.h"

namespace nrps
{
class DomainTypeC : public Domain
{
public:
    DomainTypeC(uint32_t);
    DomainTypeC(uint32_t, uint32_t);
    DomainTypeC(uint32_t, uint32_t, Configuration);
    virtual ~DomainTypeC();

    uint32_t substrate() const;
    void setSubstrate(uint32_t);

    Configuration chirality() const;
    void setChirality(Configuration);

protected:
    virtual void writeXml(xmlTextWriterPtr) const;

private:
    uint32_t m_substrate;
    Configuration m_chirality;
};
}

#endif
