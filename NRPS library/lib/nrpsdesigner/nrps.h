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
    std::string toXml() const;
    void toXml(std::ostream&) const;
    void toXml(const std::string&) const;
    void toXml(const char*) const;
    void toXml(int) const;

private:
    const std::vector<Monomer>& m_nrp;
    void toXml(xmlTextWriterPtr) const;
};
}

#endif
