#include "domaintypete.h"

using namespace nrps;

DomainTypeTe::DomainTypeTe(bool circularizing)
: AbstractDomainType(DomainType::Te), m_circularizing(circularizing)
{}

DomainTypeTe::~DomainTypeTe()
{}

bool DomainTypeTe::circularizing() const
{
    return m_circularizing;
}
