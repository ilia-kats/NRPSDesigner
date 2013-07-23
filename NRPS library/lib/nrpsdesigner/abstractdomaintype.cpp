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

bool AbstractDomainType::full() const
{
    return false;
}

std::size_t AbstractDomainType::hash() const
{
    return std::hash<int>()(static_cast<int>(type()));
}

std::size_t std::hash<AbstractDomainType>::operator()(const AbstractDomainType &d) const
{
    return d.hash();
}
