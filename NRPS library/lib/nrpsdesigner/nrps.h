#ifndef NRPSDESIGNER_NRPS_H

#include "nrpsdesigner_export.h"
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
    Nrps(const std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>>&);
    std::string toXml() const;
    void toXml(std::ostream&) const;
    void toXml(const std::string&) const;
    void toXml(const char*) const;
    void toXml(int) const;

private:
    void toXml(xmlTextWriterPtr) const;
};
}

#endif
