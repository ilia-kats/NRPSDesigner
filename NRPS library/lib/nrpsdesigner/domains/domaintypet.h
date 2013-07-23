#ifndef NRPSDESIGNER_DOMAINTYPET_H
#define NRPSDESIGNER_DOMAINTYPET_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"
#include "abstractdomaintype.h"
#include "domain.h"

namespace nrps
{
template <bool full = false>
class DomainTypeT : public DomainBaseType<full>::type
{
public:
    DomainTypeT(DomainTPosition);
    virtual ~DomainTypeT();
    
    virtual std::size_t hash() const;
    DomainTPosition position() const;
    
private:
    DomainTPosition m_position;
};

template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeT<false>::hash() const;
template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeT<true>::hash() const;

extern template class NRPSDESIGNER_EXPORT DomainTypeT<true>;
extern template class NRPSDESIGNER_EXPORT DomainTypeT<false>;
}

#endif
