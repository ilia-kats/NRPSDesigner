#include "domaintypeac.h"

using namespace nrps;

DomainTypeAC::DomainTypeAC(uint32_t substrate, Configuration chirality)
: DomainTypeA(substrate, DomainType::AC), m_chirality(chirality)
{}

DomainTypeAC::~DomainTypeAC()
{}

Configuration DomainTypeAC::chirality() const
{
    return m_chirality;
}
