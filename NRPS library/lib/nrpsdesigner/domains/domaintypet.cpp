#include "domaintypet.h"

#define POSITION_NODE "position"

using namespace nrps;

DomainTypeT::DomainTypeT(uint32_t id)
: Domain(DomainType::T, id)
{}

DomainTypeT::DomainTypeT(uint32_t id, DomainTPosition pos)
: Domain(DomainType::T, id), m_position(pos)
{}

DomainTypeT::~DomainTypeT()
{}

DomainTPosition DomainTypeT::position() const
{
    return m_position;
}

void DomainTypeT::setPosition(DomainTPosition p)
{
    m_position = p;
}

#ifdef WITH_INTERNAL_XML
void DomainTypeT::writeXml(xmlTextWriterPtr writer) const
{
    Domain::writeXml(writer);
    xmlTextWriterWriteElement(writer, BAD_CAST POSITION_NODE, BAD_CAST nrps::toString(position()).c_str());
}
#endif
