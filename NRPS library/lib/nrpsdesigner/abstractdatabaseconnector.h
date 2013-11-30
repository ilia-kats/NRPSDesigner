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
#include "monomer.h"
#include "exceptions.h"

#include <vector>

#include <boost/program_options/options_description.hpp>

namespace nrps
{
class NRPSDESIGNER_EXPORT AbstractDatabaseConnector
{
public:
    virtual ~AbstractDatabaseConnector();
    virtual boost::program_options::options_description options() = 0;
    virtual void initialize() throw (DatabaseError) = 0;
    virtual Monomer getMonomer(uint32_t) throw (DatabaseError) = 0;
    virtual Monomer getMonomer(const std::string&) throw (DatabaseError) = 0;
    virtual std::vector<Monomer> searchMonomers(const std::string&) throw (DatabaseError) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeA>> getADomains(const Monomer&, bool axoa = false) throw (DatabaseError) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeC>> getCDomains(const Monomer&, Configuration) throw (DatabaseError) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeT>> getTDomains(DomainTPosition) throw (DatabaseError) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeE>> getEDomains() throw (DatabaseError) = 0;
    virtual std::vector<std::shared_ptr<DomainTypeTe>> getTeDomains() throw (DatabaseError) = 0;
    virtual std::shared_ptr<Domain> createDomain(uint32_t) throw (DatabaseError) = 0;
    virtual void fillDomain(const std::shared_ptr<Domain>&) throw (DatabaseError) = 0;
    virtual void fillOrigin(Origin*) throw (DatabaseError) = 0;
    virtual bool isDummy(const std::shared_ptr<Domain>&) = 0;

    static AbstractDatabaseConnector* getInstance();

private:
    static AbstractDatabaseConnector *s_instance;
};
}

#endif
