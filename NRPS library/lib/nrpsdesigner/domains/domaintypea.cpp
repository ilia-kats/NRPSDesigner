#include "domaintypea.h"

using namespace nrps;

template<bool full>
DomainTypeA<full>::DomainTypeA(uint32_t substrate)
: DomainBaseType<full>::type(DomainType::A), m_substrate(substrate)
{}

template<bool full>
DomainTypeA<full>::DomainTypeA(uint32_t substrate, DomainType type)
: DomainBaseType<full>::type(type), m_substrate(substrate)
{}

template<bool full>
DomainTypeA<full>::~DomainTypeA()
{}

template<bool full>
uint32_t DomainTypeA<full>::substrate() const
{
    return m_substrate;
}

template<>
std::size_t DomainTypeA<false>::hash() const
{
    return std::hash<int>()(static_cast<int>(this->type())) ^ std::hash<uint32_t>()(this->substrate());
}

template<>
std::size_t DomainTypeA<true>::hash() const
{
    return DomainBaseType<true>::type::hash();
}

template class DomainTypeA<true>;
template class DomainTypeA<false>;
