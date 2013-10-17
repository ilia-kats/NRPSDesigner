#ifndef NRPSDESIGNER_DOMAINTYPETE_H
#define NRPSDESIGNER_DOMAINTYPETE_H

#include "nrpsdesigner_export.h"
#include "domain.h"

namespace nrps
{
class DomainTypeTe : public Domain
{
public:
    DomainTypeTe(uint32_t);
    DomainTypeTe(uint32_t, bool);
    virtual ~DomainTypeTe();

    bool circularizing() const;
    void setCircularizing(bool);

#ifdef WITH_INTERNAL_XML
// protected:
//     virtual void writeXml(xmlTextWriterPtr) const;
#endif

private:
    bool m_circularizing;
};
}

#endif
