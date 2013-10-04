#include "domaintypec.h"

#define SUBSTRATE_NODE "substrate"
#define CHIRALITY_NODE "chirality"

using namespace nrps;

DomainTypeC::DomainTypeC(uint32_t id)
: Domain(DomainType::C, id), m_substrate(0), m_chirality(Configuration::L)
{}

DomainTypeC::DomainTypeC(uint32_t id, uint32_t substrate)
: Domain(DomainType::C, id), m_substrate(substrate), m_chirality(Configuration::L)
{}

DomainTypeC::DomainTypeC(uint32_t id, uint32_t substrate, Configuration chirality)
: Domain(DomainType::C, id), m_substrate(substrate), m_chirality(chirality)
{}

DomainTypeC::~DomainTypeC()
{}

uint32_t DomainTypeC::substrate() const
{
    return m_substrate;
}

void DomainTypeC::setSubstrate(uint32_t s)
{
    m_substrate = s;
}

Configuration DomainTypeC::chirality() const
{
    return m_chirality;
}

void DomainTypeC::setChirality(Configuration c)
{
    m_chirality = c;
}

#ifdef WITH_INTERNAL_XML
void DomainTypeC::writeXml(xmlTextWriterPtr writer) const
{
    Domain::writeXml(writer);
    xmlTextWriterWriteElement(writer, BAD_CAST SUBSTRATE_NODE, BAD_CAST std::to_string(substrate()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST CHIRALITY_NODE, BAD_CAST nrps::toString(chirality()).c_str());
}
#endif
