#ifndef NRPSDESIGNER_MONOMER_H
#define NRPSDESIGNER_MONOMER_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"

#include <string>
#include <unordered_set>
#include <cstdint>

namespace nrps
{
class NRPSDESIGNER_EXPORT Monomer
{
public:
    Monomer();
    Monomer(uint32_t);

    uint32_t id() const;
    void setId(uint32_t);

    uint32_t parentId() const;
    void setParentId(uint32_t);

    uint32_t enantiomerId() const;
    void setEnantiomerId(uint32_t);

    const std::string& name() const;
    void setName(const std::string&);
    void setName(std::string&&);

    Configuration configuration() const;
    void setConfiguration(Configuration);

    const std::unordered_set<uint32_t>& modifications() const;
    void setModifications(const std::unordered_set<uint32_t>&);
    void setModifications(std::unordered_set<uint32_t>&&);
    void addModification(uint32_t);
    bool removeModification(uint32_t);
    bool hasModification(uint32_t) const;

private:
    uint32_t m_id;
    uint32_t m_parentId;
    uint32_t m_enantiomerId;
    std::string m_name;
    Configuration m_configuration;
    std::unordered_set<uint32_t> m_modifications;
};
}

#endif
