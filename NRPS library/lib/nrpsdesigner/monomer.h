#ifndef NRPSDESIGNER_MONOMER_H
#define NRPSDESIGNER_MONOMER_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"

#include <string>
#include <cstdint>

namespace nrps
{
class NRPSDESIGNER_EXPORT Monomer
{
public:
    typedef uint16_t modification_type;
    enum class Modification : modification_type {None = 0x0, Nmethyl = 0x1};

    Monomer();
    Monomer(Monomer&&);
    Monomer(uint32_t);
    Monomer(uint32_t, const std::string&, Configuration = Configuration::L, modification_type mod = static_cast<modification_type>(Modification::None));
    Monomer(uint32_t, std::string&&, Configuration = Configuration::L, modification_type mod = static_cast<modification_type>(Modification::None));
    ~Monomer();

    uint32_t id() const;
    void setId(uint32_t);

    const std::string& name() const;
    void setName(const std::string&);
    void setName(std::string&&);

    Configuration configuration() const;
    void setConfiguration(Configuration);

    modification_type modifications() const;
    void setModifications(modification_type);
    void addModification(Modification);
    bool removeModification(Modification);
    bool hasModification(Modification) const;

private:
    uint32_t m_id;
    std::string m_name;
    Configuration m_configuration;
    modification_type m_modifications;
};
}

#endif
