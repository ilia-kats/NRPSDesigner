#include "domaintypet.h"

using namespace nrps;

template<bool full>
DomainTypeT<full>::DomainTypeT(DomainTPosition pos)
: DomainBaseType<full>::type(DomainType::T), m_position(pos)
{}

template<bool full>
DomainTypeT<full>::~DomainTypeT()
{}

template<bool full>
DomainTPosition DomainTypeT<full>::position() const
{
    return m_position;
}

template<>
std::size_t DomainTypeT<false>::hash() const
{
    return std::hash<int>()(static_cast<int>(type())) ^ std::hash<int>()(static_cast<int>(position()));
}
template<>
std::size_t DomainTypeT<true>::hash() const
{
    return DomainBaseType<true>::type::hash();
}

template class DomainTypeT<true>;
template class DomainTypeT<false>;
