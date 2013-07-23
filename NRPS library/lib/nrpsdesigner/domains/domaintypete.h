#ifndef NRPSDESIGNER_DOMAINTYPETE_H
#define NRPSDESIGNER_DOMAINTYPETE_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"
#include "domain.h"

namespace nrps
{
template <bool full = false>
class DomainTypeTe : public DomainBaseType<full>::type
{
public:
    DomainTypeTe(bool);
    virtual ~DomainTypeTe();

    virtual std::size_t hash() const;
    bool circularizing() const;

private:
    bool m_circularizing;
};

template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeTe<false>::hash() const;
template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeTe<true>::hash() const;

extern template class NRPSDESIGNER_EXPORT DomainTypeTe<true>;
extern template class NRPSDESIGNER_EXPORT DomainTypeTe<false>;
}

#endif
