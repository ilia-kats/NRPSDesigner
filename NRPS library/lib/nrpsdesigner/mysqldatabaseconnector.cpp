#include "mysqldatabaseconnector.h"
#include "domaintypea.h"
#include "domaintypec.h"
#include "domaintypee.h"
#include "domaintypet.h"
#include "domaintypete.h"
#include "origin.h"
#include "product.h"

#include <unordered_map>

#include <cppconn/driver.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>

using namespace nrps;
namespace po = boost::program_options;

MySQLDatabaseConnector::MySQLDatabaseConnector()
: AbstractDatabaseConnector(), m_connection(nullptr), m_stmtMonomer(nullptr), m_stmtCoreDomainsSubstrate(nullptr), m_stmtCoreDomainsNoSubstrate(nullptr), m_stmtDomain(nullptr), m_stmtProduct(nullptr), m_stmtOrigin(nullptr)
{}

MySQLDatabaseConnector::~MySQLDatabaseConnector()
{
    if (m_connection != nullptr) {
        m_connection->close();
        delete m_connection;
    }
    if (m_stmtMonomer != nullptr) {
        delete m_stmtMonomer;
    }
    if (m_stmtCoreDomainsSubstrate != nullptr) {
        delete m_stmtCoreDomainsSubstrate;
    }
    if (m_stmtCoreDomainsNoSubstrate != nullptr) {
        delete m_stmtCoreDomainsNoSubstrate;
    }
    if (m_stmtDomain != nullptr) {
        delete m_stmtDomain;
    }
    if (m_stmtProduct != nullptr) {
        delete m_stmtProduct;
    }
    if (m_stmtOrigin != nullptr) {
        delete m_stmtOrigin;
    }
}

po::options_description MySQLDatabaseConnector::options()
{
    po::options_description options("MySQL options");
    options.add_options()("mysql-host", po::value<std::string>(&m_host)->default_value("127.0.0.1"), "host address")
                         ("mysql-port", po::value<uint16_t>(&m_port)->default_value(3306), "Port number to use for connection.")
                         ("mysql-user", po::value<std::string>(&m_user)->default_value("root"), "User for login.")
                         ("mysql-password", po::value<std::string>(&m_password), "Password to use when connecting to server.");
    return options;
}

void MySQLDatabaseConnector::initialize() throw (DatabaseError)
{
    if (testInitialized(false))
        return;
    if (m_host.empty() || m_user.empty())
        throw std::invalid_argument("MySQL database options not valid.");
    sql::Driver *driver = get_driver_instance();
    std::string server("tcp://");
    server.append(m_host).append(":").append(std::to_string(m_port));
    m_connection = driver->connect(server, m_user, m_password);
    m_connection->setSchema("nrps_designer");
    m_stmtMonomer = m_connection->prepareStatement("SELECT name, chirality, enantiomer_id, parent_id FROM `databaseInput_substrate` WHERE id = ?;");
    m_stmtCoreDomainsSubstrate = m_connection->prepareStatement("SELECT domain.id AS did, substr_xref.substrate_id AS sid, origin.id AS orid, product.id AS pid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_domain_substrateSpecificity AS substr_xref, databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype, databaseInput_product AS product) ON (domain.id = substr_xref.domain_id AND domain.cds_id = cds.id AND cds.origin_id = origin.id AND cds.product_id = product.id AND domain.domainType_id = dtype.id) WHERE (substr_xref.substrate_id = ? OR substr_xref.substrate_id = ?) AND dtype.name = ? ORDER BY IF(substr_xref.substrate_id = ?, 0, 1) ASC;");
    m_stmtCoreDomainsNoSubstrate = m_connection->prepareStatement("SELECT domain.id AS did, origin.id AS orid, product.id AS pid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype, databaseInput_product AS product) ON (domain.cds_id = cds.id AND cds.origin_id = origin.id AND cds.product_id = product.id AND domain.domainType_id = dtype.id) WHERE dtype.name = ?;");

    m_stmtDomain = m_connection->prepareStatement("SELECT d.id AS did, d.module AS dmodule, d.description AS ddesc, d.pfamLinkerStart AS dpfamlinkerstart, d.pfamLinkerStop AS dpfamlinkerstop, d.definedLinkerStart AS ddefinedlinkerstart, d.definedLinkerStop AS ddefinedlinkerstop, d.pfamStart AS dpfamstart, d.pfamStop AS dpfamstop, d.definedStart AS ddefinedstart, d.definedStop AS ddefinedstop, c.id AS cid, c.geneName AS cgenename, c.dnaSequence AS cdnaseq, c.description AS cdesc FROM databaseInput_domain AS d INNER JOIN databaseInput_cds AS c ON d.cds_id = c.id WHERE d.id = ?;");
    m_stmtProduct = m_connection->prepareStatement("SELECT id, name, description FROM databaseInput_product WHERE id = ?;");
    m_stmtOrigin = m_connection->createStatement();
}

Monomer MySQLDatabaseConnector::getMonomer(uint32_t id) throw (DatabaseError)
{
    testInitialized();
    sql::ResultSet *res = nullptr;
    try {
        m_stmtMonomer->setUInt(1, id);
        sql::ResultSet *res = m_stmtMonomer->executeQuery();
        res->next();
        Monomer m(id);
        m.setName(res->getString(1));
        m.setConfiguration(fromString<Configuration>(res->getString(2)));
        m.setEnantiomerId(res->getUInt(3));
        m.setParentId(res->getUInt(4));
        delete res;
        return m;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}

std::vector<std::shared_ptr<DomainTypeA>> MySQLDatabaseConnector::getADomains(const Monomer &m) throw (DatabaseError)
{
    return getCoreDomains<DomainTypeA>(m, toString(DomainType::A),
                                       [](DomainTypeA *d, sql::ResultSet *res){
                                           d->setSubstrate(res->getUInt("sid"));
                                       });
}

std::vector<std::shared_ptr<DomainTypeC>> MySQLDatabaseConnector::getCDomains(const Monomer &m, Configuration c) throw (DatabaseError)
{
    return getCoreDomains<DomainTypeC>(m, c == Configuration::D ? "C_D" : "C_L",
                                       [&c](DomainTypeC *d, sql::ResultSet *res){
                                           d->setSubstrate(res->getUInt("sid"));
                                           d->setChirality(c);
                                       });
}

std::vector<std::shared_ptr<DomainTypeT>> MySQLDatabaseConnector::getTDomains(DomainTPosition p) throw (DatabaseError)
{
    return getCoreDomains<DomainTypeT>(toString(DomainType::T),
                                       [&p](DomainTypeT *d, sql::ResultSet *res){
                                           d->setPosition(p);
                                       });
}

std::vector<std::shared_ptr<DomainTypeE>> MySQLDatabaseConnector::getEDomains() throw (DatabaseError)
{
    return getCoreDomains<DomainTypeE>(toString(DomainType::E),
                                       [](DomainTypeE *d, sql::ResultSet *res){});
}

std::vector<std::shared_ptr<DomainTypeTe>> MySQLDatabaseConnector::getTeDomains() throw (DatabaseError)
{
    return getCoreDomains<DomainTypeTe>(toString(DomainType::Te),
                                       [](DomainTypeTe *d, sql::ResultSet *res){});
}

bool MySQLDatabaseConnector::testInitialized(bool except) throw (DatabaseError)
{
    if (m_connection == nullptr || m_stmtMonomer == nullptr || m_stmtCoreDomainsSubstrate == nullptr || m_stmtCoreDomainsNoSubstrate == nullptr || m_stmtDomain == nullptr || m_stmtProduct == nullptr || m_stmtOrigin == nullptr)
        if (except)
            throw DatabaseError("Not initialized");
        else
            return false;
    return true;
}

template<class D, class initFunc>
std::vector<std::shared_ptr<D>> MySQLDatabaseConnector::getCoreDomains(const std::string& t, const initFunc &f) throw (DatabaseError)
{
    testInitialized();
    m_stmtCoreDomainsNoSubstrate->setString(1, t);
    return getCoreDomains<D>(m_stmtCoreDomainsNoSubstrate, f);
}

template<class D, class initFunc>
std::vector<std::shared_ptr<D>> MySQLDatabaseConnector::getCoreDomains(const Monomer &m, const std::string& t, const initFunc &f) throw (DatabaseError)
{
    testInitialized();
    try {
        m_stmtCoreDomainsSubstrate->setUInt(1, m.id());
        m_stmtCoreDomainsSubstrate->setUInt(2, m.enantiomerId());
        m_stmtCoreDomainsSubstrate->setString(3, t);
        m_stmtCoreDomainsSubstrate->setUInt(4, m.id());
        return getCoreDomains<D>(m_stmtCoreDomainsSubstrate, f);
    } catch (const sql::SQLException &e) {
        throw makeException(e);
    }
}

template<class D, class initFunc>
std::vector<std::shared_ptr<D>> MySQLDatabaseConnector::getCoreDomains(sql::PreparedStatement *stmt, const initFunc &f) throw (DatabaseError)
{
    testInitialized();
    sql::ResultSet *res = nullptr;
    try {
        std::vector<std::shared_ptr<D>> vec;
        res = stmt->executeQuery();
        if (res->rowsCount() == 0) {
            vec.emplace_back(new D(0));
        }
        else {
            while (res->next()) {
                D *d = new D(res->getUInt("did"));
                Origin *ori = d->setOrigin(res->getUInt("orid"));
                d->setProduct(res->getUInt("pid"));
                f(d, res);
                vec.emplace_back(d);
                fillOrigin(ori);
            }
        }
        delete res;
        return vec;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}

void MySQLDatabaseConnector::fillDomain(const std::shared_ptr<Domain> &d) throw (DatabaseError)
{
    if (d->id() == 0)
        return;
    sql::ResultSet *res = nullptr;
    try {
        m_stmtDomain->setUInt(1, d->id());
        res = m_stmtDomain->executeQuery();
        res->next();
        d->setModule(res->getUInt("dmodule"));
        d->setDescription(res->getString("ddesc"));
        d->setGeneName(res->getString("cgenename"));
        d->setGeneDescription(res->getString("cdesc"));
        std::string seq = res->getString("cdnaseq");
        uint32_t lstart = res->getUInt("dpfamlinkerstart"), lstop = res->getUInt("dpfamlinkerstop"), start = res->getUInt("dpfamstart"), stop = res->getUInt("dpfamstop");
        if (seq.size() > 0) {
            if (stop <= seq.size())
                d->setDnaSequencePfam(seq.substr(start, stop - start + 1));
            if (start <= seq.size())
                d->setNativePfamLinkerBefore(seq.substr(lstart, start - lstart));
            if (lstop <= seq.size())
                d->setNativePfamLinkerAfter(seq.substr(stop + 1, lstop - stop));
            lstart = res->getUInt("ddefinedlinkerstart");
            lstop = res->getUInt("ddefinedlinkerstop");
            start = res->getUInt("ddefinedstart");
            stop = res->getUInt("ddefinedstop");
            if (stop <= seq.size())
                d->setDnaSequenceDefined(seq.substr(start, stop - start + 1));
            if (lstart <= seq.size())
                d->setNativeDefinedLinkerBefore(seq.substr(lstart, start - lstart));
            if (lstop <= seq.size())
                d->setNativeDefinedLinkerAfter(seq.substr(stop + 1, lstop - stop));
        }
        delete res;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
    res = nullptr;
    try {
        if (d->product()->name().empty()) {
            m_stmtProduct->setUInt(1, d->product()->id());
            res = m_stmtProduct->executeQuery();
            res->next();
            d->product()->setName(res->getString("name"));
            d->product()->setDescription(res->getString("description"));
            delete res;
        }
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}

void MySQLDatabaseConnector::fillOrigin(Origin *ori) throw (DatabaseError)
{
    if (ori->taxId() > 0)
        return;
    sql::ResultSet *res = nullptr;
    try {
        std::string stmt("CALL get_origin_hierarchy(");
        stmt.append(std::to_string(ori->id())).append(");");
        m_stmtOrigin->execute(stmt);
        res = m_stmtOrigin->getResultSet();
        Origin *lastori = nullptr;
        while (res->next()) {
            ori = Origin::makeOrigin(res->getUInt("id"));
            ori->setSourceType(fromString<OriginSourceType>(res->getString("sourceType")));
            ori->setSource(res->getString("source"));
            ori->setSpecies(res->getString("species"));
            ori->setDescription(res->getString("description"));
            if (lastori != nullptr)
                lastori->setParent(ori);
            lastori = ori;
        }
        delete res;
        m_stmtOrigin->getMoreResults();
        delete m_stmtOrigin->getResultSet(); // extra resultset returned by CALL
    } catch (const sql::SQLException &e) {
        if (res != nullptr) {
            delete res;
            m_stmtOrigin->getMoreResults();
            delete m_stmtOrigin->getResultSet();
        }
        throw makeException(e);
    }
}

bool MySQLDatabaseConnector::isDummy(const std::shared_ptr<Domain> &d)
{
    return d->id() == 0;
}

DatabaseError MySQLDatabaseConnector::makeException(const sql::SQLException &e) const
{
    std::string what("MySQL Exception:");
    what.append(e.what());
    what.append(";    error code: ").append(std::to_string(e.getErrorCode()));
    what.append(";    SQL state: ").append(e.getSQLState());
    return DatabaseError(what);
}
