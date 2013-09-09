#ifndef NRPSDESIGNER_ORIGIN_H
#define NRPSDESIGNER_ORIGIN_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"

#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <libxml/xmlwriter.h>

namespace nrps
{
class NRPSDESIGNER_EXPORT Origin
{
public:
    ~Origin();
    uint32_t id() const;
    OriginSourceType sourceType() const;
    const std::string& source() const;
    const std::string& species() const;
    const std::string& description() const;
    uint32_t taxId() const;
    Origin* parent() const;

    void setId(uint32_t);
    void setSourceType(OriginSourceType);
    void setSource(const std::string&);
    void setSource(std::string&&);
    void setSpecies(const std::string&);
    void setSpecies(std::string&&);
    void setDescription(const std::string&);
    void setDescription(std::string&&);
    void setParent(Origin*);

    void toXml(xmlTextWriterPtr) const;

    static Origin* makeOrigin(uint32_t);

private:
    Origin(uint32_t);
    void setTaxId(uint32_t);

    static std::unordered_map<uint32_t, Origin*> s_origins;

    uint32_t m_id;
    uint32_t m_taxid;
    Origin* m_parent;
    OriginSourceType m_sourceType;
    std::string m_source;
    std::string m_species;
    std::string m_product;
    std::string m_description;
    std::unordered_set<Origin*> m_children;

};
}

#endif
