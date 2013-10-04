#include "domaintypeaoxa.h"

using namespace nrps;

DomainTypeAOxA::DomainTypeAOxA(uint32_t id)
: DomainTypeA(DomainType::AOxA, id)
{}

DomainTypeAOxA::DomainTypeAOxA(uint32_t id, uint32_t substrate)
: DomainTypeA(DomainType::AOxA, id, substrate)
{}

DomainTypeAOxA::~DomainTypeAOxA()
{}
