#include "domaintypet.h"

using namespace nrps;

DomainTypeT::DomainTypeT(DomainTPosition pos)
: AbstractDomainType(DomainType::T), m_position(pos)
{}

DomainTypeT::~DomainTypeT()
{}

DomainTPosition DomainTypeT::position() const
{
    return m_position;
}
