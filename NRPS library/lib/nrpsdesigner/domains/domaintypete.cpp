#include "domaintypete.h"

using namespace nrps;

DomainTypeTe::DomainTypeTe(uint32_t id)
: Domain(DomainType::Te, id), m_circularizing(false)
{}

DomainTypeTe::DomainTypeTe(uint32_t id, bool circularizing)
: Domain(DomainType::Te, id), m_circularizing(circularizing)
{}

DomainTypeTe::~DomainTypeTe()
{}

bool DomainTypeTe::circularizing() const
{
    return m_circularizing;
}

void DomainTypeTe::setCircularizing(bool c)
{
    m_circularizing = c;
}
