#include "abstractdomaintype.h"

using namespace nrps;

AbstractDomainType::AbstractDomainType(DomainType type)
: m_type(type)
{}

AbstractDomainType::~AbstractDomainType()
{}

DomainType AbstractDomainType::type() const
{
    return m_type;
}
