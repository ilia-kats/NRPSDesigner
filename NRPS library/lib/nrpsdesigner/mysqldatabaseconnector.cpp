#include "mysqldatabaseconnector.h"
#include "domaintypea.h"
#include "domaintypeaoxa.h"
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
: AbstractDatabaseConnector(), m_connection(nullptr), m_stmtMonomerId(nullptr), m_stmtMonomerSmash(nullptr), m_stmtMonomerSearch(nullptr), m_stmtCoreDomainsSubstrate(nullptr), m_stmtCoreDomainsNoSubstrate(nullptr), m_stmtDomain(nullptr), m_stmtProduct(nullptr), m_stmtOrigin(nullptr), m_curatedOnly(true)
{}

MySQLDatabaseConnector::~MySQLDatabaseConnector()
{
    if (m_connection != nullptr) {
        m_connection->close();
        delete m_connection;
    }
    if (m_stmtMonomerId != nullptr) {
        delete m_stmtMonomerId;
    }
    if (m_stmtMonomerSmash != nullptr) {
        delete m_stmtMonomerSmash;
    }
    if (m_stmtMonomerSearch != nullptr) {
        delete m_stmtMonomerSearch;
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
    if (m_stmtDomainSubstrate != nullptr) {
        delete m_stmtDomainSubstrate;
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
                         ("mysql-password", po::value<std::string>(&m_password), "Password to use when connecting to server.")
                         ("curated-only", po::bool_switch(&m_curatedOnly), "Use only curated domains for prediction.")
                         ("curation-group", po::value<std::string>(&m_curationGroup)->default_value("curator"), "Name of the curator user group.");
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
    m_stmtMonomerId = m_connection->prepareStatement("SELECT id, name, chirality, enantiomer_id, parent_id FROM `databaseInput_substrate` WHERE id = ?;");
    m_stmtMonomerSmash = m_connection->prepareStatement("SELECT id, name, chirality, enantiomer_id, parent_id FROM `databaseInput_substrate` WHERE smashName = ?;");
    m_stmtMonomerSearch = m_connection->prepareStatement("SELECT id, name, chirality, enantiomer_id, parent_id FROM `databaseInput_substrate` WHERE name LIKE ?;");
    std::string stmtCoreDomainsSubstrate = "SELECT domain.id AS did, domain.module AS module, substr_xref.substrate_id AS sid, origin.id AS orid, product.id AS pid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_domain_substrateSpecificity AS substr_xref, databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype, databaseInput_product AS product";
    std::string stmtCoreDomainsNoSubstrate = "SELECT domain.id AS did, domain.module AS module, origin.id AS orid, product.id AS pid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype, databaseInput_product AS product";
    if (m_curatedOnly) {
        std::string s = ", auth_user AS user, auth_group AS `group`, auth_user_groups AS user_xref";
        stmtCoreDomainsSubstrate.append(s);
        stmtCoreDomainsNoSubstrate.append(s);
    }
    stmtCoreDomainsSubstrate.append(") ON (domain.id = substr_xref.domain_id AND domain.cds_id = cds.id AND cds.origin_id = origin.id AND cds.product_id = product.id AND domain.domainType_id = dtype.id");
    stmtCoreDomainsNoSubstrate.append(") ON (domain.cds_id = cds.id AND cds.origin_id = origin.id AND cds.product_id = product.id AND domain.domainType_id = dtype.id");
    if (m_curatedOnly) {
        std::string s = " AND domain.user_id = user.id AND user.id = user_xref.user_id AND user_xref.group_id = group.id";
        stmtCoreDomainsSubstrate.append(s);
        stmtCoreDomainsNoSubstrate.append(s);
    }
    stmtCoreDomainsSubstrate.append(") WHERE (substr_xref.substrate_id = ? OR substr_xref.substrate_id = ?) AND dtype.name = ?");
    stmtCoreDomainsNoSubstrate.append(") WHERE dtype.name = ?");
    if (m_curatedOnly) {
        std::string s = " AND group.name = ?";
        stmtCoreDomainsSubstrate.append(s);
        stmtCoreDomainsNoSubstrate.append(s);
    }
    stmtCoreDomainsSubstrate.append(";");
    stmtCoreDomainsNoSubstrate.append(";");
    m_stmtCoreDomainsSubstrate = m_connection->prepareStatement(stmtCoreDomainsSubstrate);
    m_stmtCoreDomainsNoSubstrate = m_connection->prepareStatement(stmtCoreDomainsNoSubstrate);

    m_stmtDomain = m_connection->prepareStatement("SELECT d.id AS did, d.module AS dmodule, d.chirality AS chirality, d.description AS ddesc, d.pfamLinkerStart AS dpfamlinkerstart, d.pfamLinkerStop AS dpfamlinkerstop, d.definedLinkerStart AS ddefinedlinkerstart, d.definedLinkerStop AS ddefinedlinkerstop, d.pfamStart AS dpfamstart, d.pfamStop AS dpfamstop, d.definedStart AS ddefinedstart, d.definedStop AS ddefinedstop, c.id AS cid, c.geneName AS cgenename, c.dnaSequence AS cdnaseq, c.description AS cdesc, o.id AS oid, p.id AS pid, dtype.name AS type FROM databaseInput_domain AS d INNER JOIN (databaseInput_cds AS c, databaseInput_origin AS o, databaseInput_product AS p, databaseInput_type AS dtype) ON (d.cds_id = c.id AND c.origin_id = o.id AND c.product_id = p.id AND d.domainType_id = dtype.id) WHERE d.id = ?;");
    m_stmtDomainSubstrate = m_connection->prepareStatement("SELECT substrate_id AS sid FROM databaseInput_domain_substrateSpecificity WHERE domain_id = ?");
    m_stmtProduct = m_connection->prepareStatement("SELECT id, name, description FROM databaseInput_product WHERE id = ?;");
    m_stmtOrigin = m_connection->createStatement();
}

Monomer MySQLDatabaseConnector::getMonomer(uint32_t id) throw (DatabaseError)
{
    testInitialized();
    sql::ResultSet *res = nullptr;
    try {
        m_stmtMonomerId->setUInt(1, id);
        sql::ResultSet *res = m_stmtMonomerId->executeQuery();
        res->next();
        Monomer m = makeMonomer(res);
        delete res;
        return m;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}

Monomer MySQLDatabaseConnector::getMonomer(const std::string &name) throw (DatabaseError)
{
    testInitialized();
    sql::ResultSet *res = nullptr;
    try {
        m_stmtMonomerSmash->setString(1, name);
        sql::ResultSet *res = m_stmtMonomerSmash->executeQuery();
        res->next();
        Monomer m = makeMonomer(res);
        delete res;
        return m;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}

std::vector<Monomer> MySQLDatabaseConnector::searchMonomers(const std::string &name) throw (DatabaseError)
{
    testInitialized();
    sql::ResultSet *res = nullptr;
    std::string pattern = "%";
    pattern.append(name).append("%");
    std::vector<Monomer> ret;
    try {
        m_stmtMonomerSearch->setString(1, name);
        sql::ResultSet *res = m_stmtMonomerSmash->executeQuery();
        while(res->next()) {
            ret.push_back(makeMonomer(res));
        }
        delete res;
        return ret;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}

Monomer MySQLDatabaseConnector::makeMonomer(sql::ResultSet *res)
{
    Monomer m(res->getUInt(1));
    m.setName(res->getString(2));
    m.setConfiguration(fromString<Configuration>(res->getString(3)));
    m.setEnantiomerId(res->getUInt(4));
    m.setParentId(res->getUInt(5));
    return m;
}

std::vector<std::shared_ptr<DomainTypeA>> MySQLDatabaseConnector::getADomains(const Monomer &m, bool aoxa) throw (DatabaseError)
{
    auto res = [](DomainTypeA *d, sql::ResultSet *res){d->setSubstrate(res->getUInt("sid"));};
    if (aoxa)
        return getCoreDomains<DomainTypeAOxA, decltype(res), DomainTypeA>(m, toString(DomainType::AOxA), res);
    else
        return getCoreDomains<DomainTypeA>(m, toString(DomainType::A), res);
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
    return getCoreDomains<DomainTypeT>(p == DomainTPosition::BeforeC ? "Tstd" : "T_Ep",
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
    if (m_connection == nullptr || m_stmtMonomerId == nullptr || m_stmtMonomerSmash == nullptr || m_stmtMonomerSearch == nullptr|| m_stmtCoreDomainsSubstrate == nullptr || m_stmtCoreDomainsNoSubstrate == nullptr || m_stmtDomain == nullptr || m_stmtDomainSubstrate == nullptr || m_stmtProduct == nullptr || m_stmtOrigin == nullptr)
        if (except)
            throw DatabaseError("Not initialized");
        else
            return false;
    return true;
}

template<class D, class initFunc, class RD>
std::vector<std::shared_ptr<RD>> MySQLDatabaseConnector::getCoreDomains(const std::string& t, const initFunc &f) throw (DatabaseError)
{
    testInitialized();
    m_stmtCoreDomainsNoSubstrate->setString(1, t);
    if (m_curatedOnly)
        m_stmtCoreDomainsNoSubstrate->setString(2, m_curationGroup);
    return getCoreDomains<D, initFunc, RD>(m_stmtCoreDomainsNoSubstrate, f);
}

template<class D, class initFunc, class RD>
std::vector<std::shared_ptr<RD>> MySQLDatabaseConnector::getCoreDomains(const Monomer &m, const std::string& t, const initFunc &f) throw (DatabaseError)
{
    testInitialized();
    try {
        m_stmtCoreDomainsSubstrate->setUInt(1, m.id());
        m_stmtCoreDomainsSubstrate->setUInt(2, m.enantiomerId());
        m_stmtCoreDomainsSubstrate->setString(3, t);
        if (m_curatedOnly)
            m_stmtCoreDomainsSubstrate->setString(4, m_curationGroup);
        return getCoreDomains<D, initFunc, RD>(m_stmtCoreDomainsSubstrate, f);
    } catch (const sql::SQLException &e) {
        throw makeException(e);
    }
}

template<class D, class initFunc, class RD>
std::vector<std::shared_ptr<RD>> MySQLDatabaseConnector::getCoreDomains(sql::PreparedStatement *stmt, const initFunc &f) throw (DatabaseError)
{
    testInitialized();
    sql::ResultSet *res = nullptr;
    try {
        std::vector<std::shared_ptr<RD>> vec;
        res = stmt->executeQuery();
        if (!res->rowsCount()) {
            vec.emplace_back(new D(0));
        }
        else {
            while (res->next()) {
                D *d = new D(res->getUInt("did"));
                d->setModule(res->getUInt("module"));
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
    if (!d->id())
        return;
    sql::ResultSet *res = nullptr;
    try {
        m_stmtDomain->setUInt(1, d->id());
        res = m_stmtDomain->executeQuery();
        res->next();
        fillDomain(d, res);
        delete res;
    } catch (const sql::SQLException &e) {
        if (res != nullptr)
            delete res;
        throw makeException(e);
    }
}


std::shared_ptr<Domain> MySQLDatabaseConnector::createDomain(uint32_t id) throw (DatabaseError)
{
    sql::ResultSet *res = nullptr;
    sql::ResultSet *sres = nullptr;
    try {
        m_stmtDomain->setUInt(1, id);
        res = m_stmtDomain->executeQuery();
        res->next();
        std::string type = res->getString("type");
        std::shared_ptr<Domain> d;
        if (type[0] == 'C') {
            m_stmtDomainSubstrate->setUInt(1, id);
            sres = m_stmtDomainSubstrate->executeQuery();
            sres->next();
            d = std::shared_ptr<Domain>(new DomainTypeC(id, sres->getUInt("sid"), type == "C_L" ? Configuration::L : Configuration::D));
        } else if (type[0] == 'T' && (type == "T_Ep" || type == "Tstd")) {
            d = std::shared_ptr<Domain>(new DomainTypeT(id, type == "T_Ep" ? DomainTPosition::BeforeE : DomainTPosition::BeforeC));
        } else {
            switch (fromString<DomainType>(type)) {
                case DomainType::A:
                    m_stmtDomainSubstrate->setUInt(1, id);
                    sres = m_stmtDomainSubstrate->executeQuery();
                    sres->next();
                    d = std::shared_ptr<Domain>(new DomainTypeA(id, sres->getUInt("sid")));
                    break;
                case DomainType::E:
                    d = std::shared_ptr<Domain>(new DomainTypeE(id));
                    break;
                case DomainType::Te:
                    d = std::shared_ptr<Domain>(new DomainTypeTe(id));
                    break;
                case DomainType::AOxA:
                    m_stmtDomainSubstrate->setUInt(1, id);
                    sres = m_stmtDomainSubstrate->executeQuery();
                    sres->next();
                    d = std::shared_ptr<Domain>(new DomainTypeAOxA(id, sres->getUInt("sid")));
                    break;
            }
        }
        fillDomain(d, res);
        if (sres != nullptr)
            delete sres;
        delete res;
        return d;
    } catch (const sql::SQLException &e) {
        if (sres != nullptr)
            delete sres;
        if (res != nullptr)
            delete res;
        throw makeException(e);
    } catch (const DatabaseError &e) {
        if (sres != nullptr)
            delete sres;
        if (res != nullptr)
            delete res;
        throw;
    }
}

void MySQLDatabaseConnector::fillDomain(const std::shared_ptr<Domain> &d, sql::ResultSet *res) throw (DatabaseError)
{
    try {
        d->setModule(res->getUInt("dmodule"));
        d->setDescription(res->getString("ddesc"));
        d->setGeneName(res->getString("cgenename"));
        d->setGeneDescription(res->getString("cdesc"));
        std::string seq = res->getString("cdnaseq");
        uint32_t lstart = res->getUInt("dpfamlinkerstart") - 1, lstop = res->getUInt("dpfamlinkerstop") - 1, start = res->getUInt("dpfamstart") - 1, stop = res->getUInt("dpfamstop") - 1; // database 1-based
        if (seq.size()) {
            if (stop <= seq.size())
                d->setDnaSequencePfam(seq.substr(start, stop - start + 1));
            if (start <= seq.size() && !res->isNull("dpfamlinkerstart"))
                d->setNativePfamLinkerBefore(seq.substr(lstart, start - lstart));
            if (lstop <= seq.size() && !res->isNull("dpfamlinkerstop"))
                d->setNativePfamLinkerAfter(seq.substr(stop + 1, lstop - stop));
            lstart = res->getUInt("ddefinedlinkerstart") - 1;
            lstop = res->getUInt("ddefinedlinkerstop") - 1;
            start = res->getUInt("ddefinedstart") - 1;
            stop = res->getUInt("ddefinedstop") - 1;
            if (stop <= seq.size() && !res->isNull("ddefinedstart") && !res->isNull("ddefinedstop"))
                d->setDnaSequenceDefined(seq.substr(start, stop - start + 1));
            if (lstart <= seq.size() && !res->isNull("ddefinedstart") && !res->isNull("ddefinedlinkerstart"))
                d->setNativeDefinedLinkerBefore(seq.substr(lstart, start - lstart));
            if (lstop <= seq.size() && !res->isNull("ddefinedstop") && !res->isNull("ddefinedlinkerstart"))
                d->setNativeDefinedLinkerAfter(seq.substr(stop + 1, lstop - stop));
        }
        uint32_t id = res->getUInt("oid");
        if (d->origin() == nullptr || d->origin()->id() != id) {
            Origin *ori = d->setOrigin(id);
            fillOrigin(ori);
        }
        id = res->getUInt("pid");
        if (d->product() == nullptr || d->product()->id() != id) {
            d->setProduct(id);
        }
    } catch (const sql::SQLException &e) {
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
    if (ori->taxId())
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
    return !d->id();
}

DatabaseError MySQLDatabaseConnector::makeException(const sql::SQLException &e) const
{
    std::string what("MySQL Exception:");
    what.append(e.what());
    what.append("; error code: ").append(std::to_string(e.getErrorCode()));
    what.append("; SQL state: ").append(e.getSQLState());
    return DatabaseError(what);
}
