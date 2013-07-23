#ifndef NRPSDESIGNER_DOMAINTYPEE_H
#define NRPSDESIGNER_DOMAINTYPEE_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"
#include "domain.h"

namespace nrps
{
template <bool full = false>
class DomainTypeE : public DomainBaseType<full>::type
{
public:
    DomainTypeE();
    virtual ~DomainTypeE();
};

extern template class NRPSDESIGNER_EXPORT DomainTypeE<true>;
extern template class NRPSDESIGNER_EXPORT DomainTypeE<false>;
}

#endif
