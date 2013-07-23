#ifndef NRPSDESIGNER_DOMAINTYPEAC_H
#define NRPSDESIGNER_DOMAINTYPEAC_H

#include "nrpsdesigner_export.h"
#include "domaintypea.h"
#include "global_enums.h"

namespace nrps
{
template <bool full = false>
class DomainTypeAC : public DomainTypeA<full>
{
public:
    DomainTypeAC(uint32_t, Configuration);
    virtual ~DomainTypeAC();
    
    virtual size_t hash() const;
    Configuration chirality() const;
    
private:
    Configuration m_chirality;
};

template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeAC<false>::hash() const;
template<>
std::size_t NRPSDESIGNER_EXPORT DomainTypeAC<true>::hash() const;

extern template class NRPSDESIGNER_EXPORT DomainTypeAC<true>;
extern template class NRPSDESIGNER_EXPORT DomainTypeAC<false>;
}

#endif

