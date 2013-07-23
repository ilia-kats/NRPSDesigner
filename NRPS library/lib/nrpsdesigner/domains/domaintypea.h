#ifndef NRPSDESIGNER_DOMAINTYPEA_H
#define NRPSDESIGNER_DOMAINTYPEA_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"
#include "domain.h"

#include <cstdint>
#include <type_traits>

namespace nrps
{
template <bool full = false>
class DomainTypeA : public DomainBaseType<full>::type
{
public:
    DomainTypeA(uint32_t);
    virtual ~DomainTypeA();

    virtual size_t hash() const;
    uint32_t substrate() const;

protected:
    DomainTypeA(uint32_t, DomainType);

private:
    uint32_t m_substrate;
};

template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeA<false>::hash() const;
template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeA<true>::hash() const;

extern template class NRPSDESIGNER_EXPORT DomainTypeA<true>;
extern template class NRPSDESIGNER_EXPORT DomainTypeA<false>;
}

namespace std
{
template <> struct NRPSDESIGNER_EXPORT hash<nrps::DomainTypeA<false>>;
}

#endif
