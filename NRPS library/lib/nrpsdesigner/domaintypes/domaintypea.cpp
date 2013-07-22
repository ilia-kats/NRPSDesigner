#include "domaintypea.h"

using namespace nrps;

DomainTypeA::DomainTypeA(uint32_t substrate)
: AbstractDomainType(DomainType::A), m_substrate(substrate)
{}

DomainTypeA::DomainTypeA(uint32_t substrate, DomainType type)
: AbstractDomainType(type), m_substrate(substrate)
{}

DomainTypeA::~DomainTypeA()
{}

uint32_t DomainTypeA::substrate() const
{
    return m_substrate;
}
