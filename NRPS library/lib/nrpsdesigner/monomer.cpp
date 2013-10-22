#include "monomer.h"

using namespace nrps;

Monomer::Monomer()
: m_id(0), m_parentId(0), m_enantiomerId(0), m_name(), m_configuration(Configuration::L)
{}

Monomer::Monomer(uint32_t id)
: m_id(id), m_parentId(0), m_enantiomerId(0), m_name(), m_configuration(Configuration::L)
{}

uint32_t Monomer::id() const
{
    return m_id;
}

void Monomer::setId(uint32_t id)
{
    m_id = id;
}

uint32_t Monomer::parentId() const
{
    return m_parentId;
}

void Monomer::setParentId(uint32_t id)
{
    m_parentId = id;
}

uint32_t Monomer::enantiomerId() const
{
    return m_enantiomerId;
}

void Monomer::setEnantiomerId(uint32_t id)
{
    m_enantiomerId = id;
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

Configuration Monomer::configuration() const
{
    return m_configuration;
}

void Monomer::setConfiguration(Configuration conf)
{
    m_configuration = conf;
}

const std::unordered_set<uint32_t>& Monomer::modifications() const
{
    return m_modifications;
}

void Monomer::setModifications(const std::unordered_set<uint32_t>& mods)
{
    m_modifications = mods;
}

void Monomer::setModifications(std::unordered_set<uint32_t>&& mods)
{
    m_modifications = std::move(mods);
}

void Monomer::addModification(uint32_t mod)
{
    m_modifications.insert(mod);
}

bool Monomer::removeModification(uint32_t mod)
{
    return m_modifications.erase(mod) > 0;
}

bool Monomer::hasModification(uint32_t mod) const
{
    return m_modifications.count(mod) > 0;
}
