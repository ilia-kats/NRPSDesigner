#ifndef NRPSDESIGNER_MYSQLDATABASECONNECTOR_H
#define NRPSDESIGNER_MYSQLDATABASECONNECTOR_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"
#include "abstractdatabaseconnector.h"

namespace sql
{
class Connection;
class Statement;
class PreparedStatement;
class ResultSet;
}

namespace nrps
{
class Origin;
class MySQLDatabaseConnector : public AbstractDatabaseConnector
{
public:
    MySQLDatabaseConnector();
    virtual ~MySQLDatabaseConnector();
    virtual void initialize();
    virtual Monomer getMonomer(uint32_t);
    virtual std::vector<std::shared_ptr<DomainTypeA>> getADomains(const Monomer&);
    virtual std::vector<std::shared_ptr<DomainTypeC>> getCDomains(const Monomer&, Configuration);
    virtual std::vector<std::shared_ptr<DomainTypeT>> getTDomains(DomainTPosition);
    virtual std::vector<std::shared_ptr<DomainTypeTe>> getTeDomains();
    virtual std::vector<std::shared_ptr<DomainTypeE>> getEDomains();
    virtual void fillDomain(std::shared_ptr<Domain>);
    virtual void fillOrigin(Origin*);

private:
    bool testInitialized(bool except = true);
    template<class D, class initFunc>
    std::vector<std::shared_ptr<D>> getCoreDomains(const Monomer&, bool, DomainType, Configuration, const initFunc&);

    sql::Connection *m_connection;
    sql::PreparedStatement *m_stmtMonomer;
    sql::PreparedStatement *m_stmtCoreDomainsId;
    sql::PreparedStatement *m_stmtCoreDomainsEnantId;
    sql::PreparedStatement *m_stmtCoreDomains;
    sql::PreparedStatement *m_stmtDomain;
    sql::Statement *m_stmtOrigin;
};
}

#endif
