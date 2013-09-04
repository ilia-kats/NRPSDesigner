#ifndef NRPSDESIGNER_NRPS_H

#include "nrpsdesigner_export.h"
//#include "global_enums.h"
#include "monomer.h"
#include "domain.h"
#include "taxon.h"

#include <vector>
#include <string>
#include <ostream>

#include <libxml/xmlwriter.h>

namespace nrps
{
class Node;
class AbstractDatabaseConnector;
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
    float makeWeight(const Taxon&, const Taxon&);
    const std::vector<Monomer>& m_nrp;
    void toXml(xmlTextWriterPtr) const;
    void buildModule(int, AbstractDatabaseConnector*, std::vector<std::shared_ptr<std::vector<Node*>>>&);
};
}

#endif
