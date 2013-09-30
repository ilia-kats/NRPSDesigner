#include "domaintypea.h"

#define SUBSTRATE_NODE "substrate"

using namespace nrps;

DomainTypeA::DomainTypeA(uint32_t id)
: Domain(DomainType::A, id)
{}

DomainTypeA::DomainTypeA(uint32_t id, uint32_t substrate)
: Domain(DomainType::A, id), m_substrate(substrate)
{}

DomainTypeA::DomainTypeA(DomainType t, uint32_t id)
: Domain(t, id)
{}

DomainTypeA::DomainTypeA(DomainType t, uint32_t id, uint32_t substrate)
: Domain(t, id), m_substrate(substrate)
{}

DomainTypeA::~DomainTypeA()
{}

uint32_t DomainTypeA::substrate() const
{
    return m_substrate;
}

void DomainTypeA::setSubstrate(uint32_t s)
{
    m_substrate = s;
}

void DomainTypeA::writeXml(xmlTextWriterPtr writer) const
{
    Domain::writeXml(writer);
    xmlTextWriterWriteElement(writer, BAD_CAST SUBSTRATE_NODE, BAD_CAST std::to_string(substrate()).c_str());
}

