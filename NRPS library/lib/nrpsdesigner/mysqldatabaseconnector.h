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
    virtual Monomer getMonomer(const std::string&) throw (DatabaseError);
    virtual std::vector<Monomer> searchMonomers(const std::string&) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeA>> getADomains(const Monomer&, bool aoxa = false) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeC>> getCDomains(const Monomer&, Configuration) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeT>> getTDomains(DomainTPosition) throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeTe>> getTeDomains() throw (DatabaseError);
    virtual std::vector<std::shared_ptr<DomainTypeE>> getEDomains() throw (DatabaseError);
    virtual std::shared_ptr<Domain> createDomain(uint32_t) throw (DatabaseError);
    virtual void fillDomain(const std::shared_ptr<Domain>&) throw (DatabaseError);
    virtual void fillOrigin(Origin*) throw (DatabaseError);
    virtual bool  isDummy(const std::shared_ptr<Domain>&);

private:
    bool testInitialized(bool except = true) throw (DatabaseError);
    template<class D, class initFunc, class RD=D>
    std::vector<std::shared_ptr<RD>> getCoreDomains(const std::string&, const initFunc&) throw (DatabaseError);
    template<class D, class initFunc, class RD=D>
    std::vector<std::shared_ptr<RD>> getCoreDomains(const Monomer&, const std::string&, const initFunc&) throw (DatabaseError);
    template<class D, class initFunc, class RD=D>
    std::vector<std::shared_ptr<RD>> getCoreDomains(sql::PreparedStatement*, const initFunc&) throw (DatabaseError);
    Monomer makeMonomer(sql::ResultSet*);
    void fillDomain(const std::shared_ptr<Domain>&, sql::ResultSet*) throw (DatabaseError);
    DatabaseError makeException(const sql::SQLException &e) const;

    sql::Connection *m_connection;
    sql::PreparedStatement *m_stmtMonomerId;
    sql::PreparedStatement *m_stmtMonomerSmash;
    sql::PreparedStatement *m_stmtMonomerSearch;
    sql::PreparedStatement *m_stmtCoreDomainsSubstrate;
    sql::PreparedStatement *m_stmtCoreDomainsNoSubstrate;
    sql::PreparedStatement *m_stmtDomain;
    sql::PreparedStatement *m_stmtDomainSubstrate;
    sql::PreparedStatement *m_stmtProduct;
    sql::Statement *m_stmtOrigin;

    bool m_curatedOnly;
    std::string m_host;
    uint16_t m_port;
    std::string m_user;
    std::string m_password;
    std::string m_curationGroup;
};
}

#endif
