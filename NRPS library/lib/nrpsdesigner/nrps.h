#ifndef NRPSDESIGNER_NRPS_H
#define NRPSDESIGNER_NRPS_H

#include "config.h"
#include "nrpsdesigner_export.h"
#include "monomer.h"
#include "domain.h"

#include <vector>
#include <string>
#include <ostream>

#ifdef WITH_INTERNAL_XML
#include <libxml/xmlwriter.h>
#endif
#ifdef WITH_SBOL
extern "C" {
#include <sbol.h>
}
#endif

namespace nrps
{
class NRPSDESIGNER_EXPORT Nrps : public std::vector<std::shared_ptr<Domain>>
{
public:
    Nrps();
    bool isIndigoidineTagged() const;
    void setIndigoidineTagged(bool);

#ifdef WITH_INTERNAL_XML
    std::string toXml() const;
    void toXml(std::ostream&) const;
    void toXml(const std::string&) const;
    void toXml(const char*) const;
    void toXml(int) const;
    void toXml(xmlTextWriterPtr) const;
#endif

#ifdef WITH_SBOL
    std::string toSbol() const;
    void toSbol(std::ostream&) const;
    void toSbol(const std::string&) const;
    void toSbol(const char*) const;
    void toSbol(int) const;
    DNAComponent* toSbol(Document*) const;
#endif

private:
    bool m_indigoidineTagged;
};
}

#endif
