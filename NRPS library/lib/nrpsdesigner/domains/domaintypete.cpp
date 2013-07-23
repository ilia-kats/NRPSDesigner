#include "domaintypete.h"

using namespace nrps;

template<bool full>
DomainTypeTe<full>::DomainTypeTe(bool circularizing)
: DomainBaseType<full>::type(DomainType::Te), m_circularizing(circularizing)
{}

template<bool full>
DomainTypeTe<full>::~DomainTypeTe()
{}

template<bool full>
bool DomainTypeTe<full>::circularizing() const
{
    return m_circularizing;
}

template<>
std::size_t DomainTypeTe<false>::hash() const
{
    return std::hash<int>()(static_cast<int>(type())) ^ std::hash<bool>()(circularizing());
}
template<>
std::size_t DomainTypeTe<true>::hash() const
{
    return DomainBaseType<true>::type::hash();
}

template class DomainTypeTe<true>;
template class DomainTypeTe<false>;
