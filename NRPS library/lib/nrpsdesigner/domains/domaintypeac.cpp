#include "domaintypeac.h"

using namespace nrps;

template<bool full>
DomainTypeAC<full>::DomainTypeAC(uint32_t substrate, Configuration chirality)
: DomainTypeA<full>(substrate, DomainType::AC), m_chirality(chirality)
{}

template<bool full>
DomainTypeAC<full>::~DomainTypeAC()
{}

template<bool full>
Configuration DomainTypeAC<full>::chirality() const
{
    return m_chirality;
}

template<>
std::size_t DomainTypeAC<false>::hash() const
{
    return std::hash<int>()(static_cast<int>(type())) ^ std::hash<uint32_t>()(substrate()) ^ std::hash<int>()(static_cast<int>(chirality()));
}
template<>
std::size_t DomainTypeAC<true>::hash() const
{
    return DomainTypeA<true>::hash();
}

template class NRPSDESIGNER_EXPORT DomainTypeAC<true>;
template class NRPSDESIGNER_EXPORT DomainTypeAC<false>;
