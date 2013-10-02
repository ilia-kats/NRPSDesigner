#ifndef NRPSDESIGNER_NRPS_H
#define NRPSDESIGNER_NRPS_H

#include "nrpsdesigner_export.h"
#include "monomer.h"
#include "domain.h"

#include <vector>
#include <string>
#include <ostream>

#include <libxml/xmlwriter.h>

namespace nrps
{
class NRPSDESIGNER_EXPORT Nrps : public std::vector<std::shared_ptr<Domain>>
{
public:
    Nrps(const std::vector<Monomer>&);
    bool isIndigoidineTagged() const;
    void setIndigoidineTagged(bool);
    std::string toXml() const;
    void toXml(std::ostream&) const;
    void toXml(const std::string&) const;
    void toXml(const char*) const;
    void toXml(int) const;

#ifdef WITH_SBOL
    std::string toSbol() const;
    void toSbol(std::ostream&) const;
    void toSbol(const std::string&) const;
    void toSbol(const char*) const;
    void toSbol(int) const;
#endif

private:
    const std::vector<Monomer>& m_nrp;
    bool m_indigoidineTagged;
    void toXml(xmlTextWriterPtr) const;
};
}

#endif
