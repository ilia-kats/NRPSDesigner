#include "domaintypee.h"
#include "global_enums.h"

using namespace nrps;

template<bool full>
DomainTypeE<full>::DomainTypeE()
: DomainBaseType<full>::type(DomainType::E)
{}

template<bool full>
DomainTypeE<full>::~DomainTypeE()
{}

template class DomainTypeE<true>;
template class DomainTypeE<false>;
