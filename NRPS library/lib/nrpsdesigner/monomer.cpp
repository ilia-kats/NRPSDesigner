#include "monomer.h"

using namespace nrps;

Monomer::Monomer()
: m_id(-1), m_name(), m_configuration(Configuration::L), m_modifications(static_cast<modification_type>(Modification::None))
{}

Monomer::Monomer(Monomer&& other)
: m_id(other.m_id), m_name(std::move(other.m_name)), m_configuration(other.m_configuration), m_modifications(other.m_modifications)
{}

Monomer::Monomer(int id)
: m_id(id), m_name(), m_configuration(Configuration::L), m_modifications(static_cast<modification_type>(Modification::None))
{}

Monomer::Monomer(int id, const std::string &name, Configuration conf, modification_type mod)
: m_id(id), m_name(name), m_configuration(conf), m_modifications(mod)
{}

Monomer::Monomer(int id, std::string&& name, Configuration conf, modification_type mod)
: m_id(id), m_name(std::move(name)), m_configuration(conf), m_modifications(mod)
{}

Monomer::~Monomer()
{}

int Monomer::id() const
{
    return m_id;
}

void Monomer::setId(int id)
{
    m_id = id;
}

const std::string& Monomer::name() const
{
    return m_name;
}

void Monomer::setName(const std::string &name)
{
    m_name = name;
}

void Monomer::setName(std::string &&name)
{
    m_name = std::move(name);
}

Monomer::Configuration Monomer::configuration() const
{
    return m_configuration;
}

void Monomer::setConfiguration(Monomer::Configuration conf)
{
    m_configuration = conf;
}

Monomer::modification_type Monomer::modifications() const
{
    return m_modifications;
}

void Monomer::setModifications(modification_type mods)
{
    m_modifications = mods;
}

void Monomer::addModification(Monomer::Modification mod)
{
    m_modifications |= static_cast<modification_type>(mod);
}

bool Monomer::removeModification(Monomer::Modification mod)
{
    modification_type modi = static_cast<modification_type>(mod);
    bool ret = m_modifications & modi;
    m_modifications &= ~modi;
    return ret;
}
