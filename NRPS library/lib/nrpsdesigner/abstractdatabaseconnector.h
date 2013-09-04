#ifndef NRPSDESIGNER_ABSTRACTDATABASECONNECTOR_H
#define NRPSDESIGNER_ABSTRACTDATABASECONNECTOR_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"
#include "domain.h"
#include "domains/domaintypea.h"
#include "domains/domaintypec.h"
#include "domains/domaintypet.h"
#include "domains/domaintypee.h"
#include "domains/domaintypete.h"
#include "generalpathway.h"
#include "monomer.h"

#include <vector>

namespace nrps
{
class NRPSDESIGNER_EXPORT AbstractDatabaseConnector
{
public:
    virtual ~AbstractDatabaseConnector();
    virtual void initialize() = 0; // TODO: write class for parameters (host, user, pw)
    virtual Monomer getMonomer(uint32_t) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeA>> getADomains(const Monomer&) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeC>> getCDomains(const Monomer&, Configuration) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeT>> getTDomains(DomainTPosition) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeE>> getEDomains() = 0;
    virtual std::vector<std::shared_ptr<DomainTypeTe>> getTeDomains() = 0;
    virtual std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> getPotentialDomains(const GeneralPathway &pathway) = 0;
    virtual void fillDomain(std::shared_ptr<Domain>) = 0;
    virtual void fillOrigin(Origin*) = 0;

    static AbstractDatabaseConnector* getInstance();

private:
    static AbstractDatabaseConnector *s_instance;
};
}

#endif
