#ifndef NRPSDESIGNER_MYSQLDATABASECONNECTOR_H
#define NRPSDESIGNER_MYSQLDATABASECONNECTOR_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"
#include "abstractdatabaseconnector.h"

#include <stdexcept>

#include <cppconn/exception.h>

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
class NRPSDESIGNER_EXPORT MySQLDatabaseConnector : public AbstractDatabaseConnector
{
public:
    MySQLDatabaseConnector();
    virtual ~MySQLDatabaseConnector();
    virtual boost::program_options::options_description options();
    virtual void initialize() throw (DatabaseError);
    virtual Monomer getMonomer(uint32_t) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeA>> getADomains(const Monomer&) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeC>> getCDomains(const Monomer&, Configuration) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeT>> getTDomains(DomainTPosition) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeTe>> getTeDomains() throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeE>> getEDomains() throw (DatabaseError);
    virtual void fillDomain(std::shared_ptr<Domain>) throw (DatabaseError);
    virtual void fillOrigin(Origin*) throw (DatabaseError);

private:
    bool testInitialized(bool except = true) throw (DatabaseError);
    template<class D, class initFunc>
    std::vector<std::shared_ptr<D>> getCoreDomains(DomainType, const initFunc&) throw (DatabaseError);
    template<class D, class initFunc>
    std::vector<std::shared_ptr<D>> getCoreDomains(const Monomer&, DomainType, Configuration, const initFunc&) throw (DatabaseError);
    template<class D, class initFunc>
    std::vector<std::shared_ptr<D>> getCoreDomains(sql::PreparedStatement*, const initFunc&) throw (DatabaseError);
    DatabaseError makeException(const sql::SQLException &e) const;

    sql::Connection *m_connection;
    sql::PreparedStatement *m_stmtMonomer;
    sql::PreparedStatement *m_stmtCoreDomainsSubstrate;
    sql::PreparedStatement *m_stmtCoreDomainsNoSubstrate;
    sql::PreparedStatement *m_stmtDomain;
    sql::PreparedStatement *m_stmtProduct;
    sql::Statement *m_stmtOrigin;

    std::string m_host;
    uint16_t m_port;
    std::string m_user;
    std::string m_password;
};
}

#endif
