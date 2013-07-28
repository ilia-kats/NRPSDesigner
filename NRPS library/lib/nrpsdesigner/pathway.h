#ifndef NRPSDESIGNER_PATHWAY_H
#define NRPSDESIGNER_PATHWAY_H

#include "nrpsdesigner_export.h"

#include <cstdint>
#include <string>

namespace nrps
{
class NRPSDESIGNER_EXPORT Pathway
{
public:
    // move constructor will be implicitly generated if no user-defined destructor is present
    Pathway() = default;
    Pathway(uint32_t id);
    Pathway(uint32_t pathwayId, uint32_t taxId, const std::string &pathway, const std::string &linkout, const std::string &uniprotId, const std::string &norineId, const std::string &description);

    uint32_t pathwayId() const;
    uint32_t taxId() const;
    const std::string& pathway() const;
    const std::string& linkout() const;
    const std::string& uniProtId() const;
    const std::string& norineId() const;
    const std::string& description() const;

    void setPathwayId(uint32_t);
    void setTaxId(uint32_t);
    void setPathway(const std::string&);
    void setPathway(std::string&&);
    void setLinkout(const std::string&);
    void setLinkout(std::string&&);
    void setUniProtId(const std::string&);
    void setUniProtId(std::string&&);
    void setNorineId(const std::string&);
    void setNorineId(std::string&&);
    void setDescription(const std::string&);
    void setDescription(std::string&&);

private:
    uint32_t m_pathwayId;
    uint32_t m_taxId;
    std::string m_pathway;
    std::string m_linkout;
    std::string m_uniProtId;
    std::string m_norineId;
    std::string m_description;
};
}

#endif
