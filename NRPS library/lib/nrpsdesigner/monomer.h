#ifndef NRPSDESIGNER_MONOMER_H
#define NRPSDESIGNER_MONOMER_H

#include "nrpsdesigner_export.h"

#include <string>

namespace nrps
{
class NRPSDESIGNER_EXPORT Monomer
{
public:
    typedef uint16_t modification_type;
    enum class Configuration {L, D};
    enum class Modification : modification_type {None = 0x0, Nmethyl = 0x1};
    
    Monomer();
    Monomer(Monomer&&);
    Monomer(int);
    Monomer(int, const std::string&, Configuration = Configuration::L, modification_type mod = static_cast<modification_type>(Modification::None));
    Monomer(int, std::string&&, Configuration = Configuration::L, modification_type mod = static_cast<modification_type>(Modification::None));
    ~Monomer();
    
    int id() const;
    void setId(int);
    
    const std::string& name() const;
    void setName(const std::string&);
    void setName(std::string&&);
    
    Configuration configuration() const;
    void setConfiguration(Configuration);
    
    modification_type modifications() const;
    void setModifications(modification_type);
    void addModification(Modification);
    bool removeModification(Modification);
    
private:
    int m_id;
    std::string m_name;
    Configuration m_configuration;
    modification_type m_modifications;
};
}

#endif
