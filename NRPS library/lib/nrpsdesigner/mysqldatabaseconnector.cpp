#include "mysqldatabaseconnector.h"
#include "domaintypea.h"
#include "domaintypec.h"
#include "domaintypee.h"
#include "domaintypet.h"
#include "domaintypete.h"
#include "origin.h"

#include <unordered_map>
#include <stdexcept>

#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>

using namespace nrps;

MySQLDatabaseConnector::MySQLDatabaseConnector()
: AbstractDatabaseConnector(), m_connection(nullptr), m_stmtMonomer(nullptr), m_stmtCoreDomains(nullptr), m_stmtInitiationADomains(nullptr), m_stmtElongationACDomains(nullptr), m_stmtTDomains(nullptr), m_stmtEDomains(nullptr), m_stmtTeDomains(nullptr)
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
    if (m_stmtCoreDomains != nullptr) {
        delete m_stmtCoreDomains;
    }
    if (m_stmtOrigin != nullptr) {
        delete m_stmtOrigin;
    }
    if (m_stmtInitiationADomains != nullptr) {
        delete m_stmtInitiationADomains;
    }
    if (m_stmtElongationACDomains != nullptr) {
        delete m_stmtElongationACDomains;
    }
    if (m_stmtTDomains != nullptr) {
        delete m_stmtTDomains;
    }
    if (m_stmtEDomains != nullptr) {
        delete m_stmtEDomains;
    }
    if (m_stmtTeDomains != nullptr) {
        delete m_stmtTeDomains;
    }
}

void MySQLDatabaseConnector::initialize()
{
    if (testInitialized(false))
        return;
    sql::Driver *driver = get_driver_instance();
    m_connection = driver->connect("tcp://127.0.0.1:3306", "root", "");
    m_connection->setSchema("nrps_designer");
    m_stmtMonomer = m_connection->prepareStatement("SELECT name, chirality, enantiomer_id, parent_id FROM `databaseInput_substrate` WHERE id = ?;");
    m_stmtCoreDomainsId = m_connection->prepareStatement("SET @id := ?;");
    m_stmtCoreDomainsEnantId = m_connection->prepareStatement("SET @enant_id := ?;");
    m_stmtCoreDomains = m_connection->prepareStatement("SELECT domain.id AS did, substr_xref.substrate_id AS sid, origin.id AS orid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_domain_substrateSpecificity AS substr_xref, databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype) ON (domain.id = substr_xref.id AND domain.cds_id = cds.id AND cds.origin_id = origin.id AND domain.domainType_id = dtype.id) WHERE (? OR substr_xref.substrate_id = @id OR substr_xref.substrate_id = @enant_id) AND dtype.name = ? AND domain.chirality = ? ORDER BY IF(substr_xref.substrate_id = @id, 0, 1) ASC;");

    m_stmtDomain = m_connection->prepareStatement("SELECT d.id AS did, d.module AS dmodule, d.description AS ddesc, d.pfamLinkerStart AS dpfamlinkerstart, d.pfamLinkerStop AS dpfamlinkerstop, d.definedLinkerStart AS ddefinedlinkerstart, d.definedLinkerStop AS ddefinedlinkerstop, d.pfamStart AS dpfamstart, d.pfamStop AS dpfamstop, d.definedStart AS ddefinedstart, d.definedStop AS ddefinedstop, c.id AS cid, c.geneName AS cgenename, c.dnaSequence AS cdnaseq, c.description AS cdesc FROM databaseInput_domain AS d INNER JOIN databaseInput_cds AS c ON d.cds_id = c.id WHERE d.id = ?;");
    m_stmtOrigin = m_connection->createStatement();
    /*m_stmtInitiationADomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_a` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE domain.substrate_specificity_aa_id = ? AND (? OR domain.bb_id != '');");
    m_stmtElongationACDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_a_c` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE domain.substrate_specificity_aa_id = ? AND domain.chirality = ? AND (? OR domain.bb_id != '');");
    m_stmtTDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_t` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE domain.positioning = ? AND (? OR domain.bb_id != '');");
    m_stmtEDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_e` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE ? OR domain.bb_id != '';");
    m_stmtTeDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_te` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE ? OR domain.bb_id != '';");*/
}

Monomer MySQLDatabaseConnector::getMonomer(uint32_t id)
{
    testInitialized();
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
}

std::vector<std::shared_ptr<DomainTypeA>> MySQLDatabaseConnector::getADomains(const Monomer &m)
{
    return getCoreDomains<DomainTypeA>(m, false, DomainType::A, Configuration::None,
                                       [](DomainTypeA *d, sql::ResultSet *res){
                                           d->setSubstrate(res->getUInt("sid"));
                                       });
}

std::vector<std::shared_ptr<DomainTypeC>> MySQLDatabaseConnector::getCDomains(const Monomer &m, Configuration c)
{
    return getCoreDomains<DomainTypeC>(m, false, DomainType::C, c,
                                       [&c](DomainTypeC *d, sql::ResultSet *res){
                                           d->setSubstrate(res->getUInt("sid"));
                                           d->setChirality(c);
                                       });
}

std::vector<std::shared_ptr<DomainTypeT>> MySQLDatabaseConnector::getTDomains(DomainTPosition p)
{
    return getCoreDomains<DomainTypeT>(Monomer(), true, DomainType::T, Configuration::None,
                                       [&p](DomainTypeT *d, sql::ResultSet *res){
                                           d->setPosition(p);
                                       });
}

std::vector<std::shared_ptr<DomainTypeE>> MySQLDatabaseConnector::getEDomains()
{
    return getCoreDomains<DomainTypeE>(Monomer(), true, DomainType::E, Configuration::None,
                                       [](DomainTypeE *d, sql::ResultSet *res){});
}

std::vector<std::shared_ptr<DomainTypeTe>> MySQLDatabaseConnector::getTeDomains()
{
    return getCoreDomains<DomainTypeTe>(Monomer(), true, DomainType::Te, Configuration::None,
                                       [](DomainTypeTe *d, sql::ResultSet *res){});
}

std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> MySQLDatabaseConnector::getPotentialDomains(const GeneralPathway &pathway)
{
    testInitialized();
    std::unordered_map<AbstractDomainType*, std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> cache;
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> potentialDomains;
    /*for (const std::shared_ptr<AbstractDomainType> &domain : pathway) {
        if (cache.count(domain.get()) > 0)
            potentialDomains.push_back(cache[domain.get()]);
        else {
            const AbstractDomainType *d = domain.get();
            std::shared_ptr<std::vector<std::shared_ptr<Domain>>> fullDomain;
            switch (domain->type()) {
                case DomainType::A:
                    if (domain->full()) {
                        fullDomain = getADomains<true>(d);
                    } else {
                        fullDomain = getADomains<false>(d);
                    }
                    break;
                case DomainType::AC:
                    if (domain->full()) {
                        fullDomain = getACDomains<true>(d);
                    } else {
                        fullDomain = getACDomains<false>(d);
                    }
                    break;
                case DomainType::T:
                    if (domain->full()) {
                        fullDomain = getTDomains<true>(d);
                    } else {
                        fullDomain = getTDomains<false>(d);
                    }
                    break;
                case DomainType::E:
                    if (domain->full()) {
                        fullDomain = getEDomains<true>(d);
                    } else {
                        fullDomain = getEDomains<false>(d);
                    }
                    break;
                case DomainType::Te:
                    if (domain->full()) {
                        fullDomain = getTeDomains<true>(d);
                    } else {
                        fullDomain = getTeDomains<false>(d);
                    }
                    break;
            }
            potentialDomains.push_back(fullDomain);
        }
    }*/
    return potentialDomains;
}

bool MySQLDatabaseConnector::testInitialized(bool except)
{
    if (m_connection == nullptr || m_stmtMonomer == nullptr || m_stmtCoreDomains == nullptr || m_stmtOrigin == nullptr/* || m_stmtInitiationADomains == nullptr || m_stmtElongationACDomains == nullptr || m_stmtTDomains == nullptr || m_stmtEDomains == nullptr || m_stmtTeDomains == nullptr*/)
        if (except)
            throw std::logic_error("Not initialized");
        else
            return false;
    return true;
}

template<class D, class initFunc>
std::vector<std::shared_ptr<D>> MySQLDatabaseConnector::getCoreDomains(const Monomer &m, bool substr, DomainType t, Configuration c, const initFunc &f)
{
    testInitialized();
    std::vector<std::shared_ptr<D>> vec;
    m_stmtCoreDomainsId->setUInt(1, m.id());
    m_stmtCoreDomainsId->execute();
    m_stmtCoreDomainsEnantId->setUInt(1, m.enantiomerId());
    m_stmtCoreDomainsEnantId->execute();
    m_stmtCoreDomains->setBoolean(1, substr);
    m_stmtCoreDomains->setString(2, toString(t));
    m_stmtCoreDomains->setString(3, toString(c));
    //if (p != DomainTPosition::None)
    sql::ResultSet *res = m_stmtCoreDomains->executeQuery();
    while (res->next()) {
        D *d = new D(res->getUInt("did"));
        d->setOrigin(res->getUInt("orid"));
        f(d, res);
        vec.emplace_back(d);
    }
    delete res;
    return vec;
}

void MySQLDatabaseConnector::fillDomain(std::shared_ptr<Domain> d)
{
    m_stmtDomain->setUInt(1, d->id());
    sql::ResultSet *res = m_stmtDomain->executeQuery();
    res->next();
    d->setModule(res->getUInt("dmodule"));
    d->setDescription(res->getString("ddesc"));
    d->setGeneName(res->getString("cgenename"));
    d->setGeneDescription(res->getString("cdesc"));
    std::string seq = res->getString("cdnaseq");
    uint32_t lstart = res->getUInt("dpfamlinkerstart"), lstop = res->getUInt("dpfamlinkerstop"), start = res->getUInt("dpfamstart"), stop = res->getUInt("dpfamstop");
    d->setDnaSequencePfam(seq.substr(start, stop - start + 1));
    d->setNativePfamLinkerBefore(seq.substr(lstart, start - lstart));
    d->setNativePfamLinkerAfter(seq.substr(stop + 1, lstop - stop));
    lstart = res->getUInt("ddefinedlinkerstart");
    lstop = res->getUInt("ddefinedlinkerstop");
    start = res->getUInt("ddefinedstart");
    stop = res->getUInt("ddefinedstop");
    d->setDnaSequenceDefined(seq.substr(start, stop - start + 1));
    d->setNativeDefinedLinkerBefore(seq.substr(lstart, start - lstart));
    d->setNativeDefinedLinkerAfter(seq.substr(stop + 1, lstop - stop));
}

void MySQLDatabaseConnector::fillOrigin(Origin *ori)
{
    if (ori->taxId() > 0)
        return;
    std::string stmt("CALL get_origin_hierarchy(");
    stmt.append(std::to_string(ori->id())).append(");");
    m_stmtOrigin->execute(stmt);
    sql::ResultSet *res = m_stmtOrigin->getResultSet();
    Origin *lastori = nullptr;
    while (res->next()) {
        ori = Origin::makeOrigin(res->getUInt("id"));
        ori->setSourceType(fromString<OriginSourceType>(res->getString("sourceType")));
        ori->setSource(res->getString("source"));
        ori->setSpecies(res->getString("species"));
        ori->setProduct(res->getString("product"));
        ori->setDescription(res->getString("description"));
        if (lastori != nullptr)
            lastori->setParent(ori);
        lastori = ori;
    }
    delete res;
    m_stmtOrigin->getMoreResults();
    delete m_stmtOrigin->getResultSet(); // extra resultset returned by CALL
}

/*template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getADomains(const AbstractDomainType *domain)
{
    const DomainTypeA<full> *d = dynamic_cast<const DomainTypeA<full>*>(domain);
    m_stmtInitiationADomains->setUInt(1, d->substrate());
    m_stmtInitiationADomains->setBoolean(2, true);
    sql::ResultSet *res = m_stmtInitiationADomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeA>(res, d->substrate());
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getACDomains(const AbstractDomainType *domain)
{
    const DomainTypeAC<full> *d = dynamic_cast<const DomainTypeAC<full>*>(domain);
    m_stmtElongationACDomains->setUInt(1, d->substrate());
    m_stmtElongationACDomains->setBoolean(2, d->chirality() == Configuration::D);
    m_stmtElongationACDomains->setBoolean(3, true);
    sql::ResultSet *res = m_stmtElongationACDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeAC>(res, d->substrate(), d->chirality());
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getTDomains(const AbstractDomainType *domain)
{
    const DomainTypeT<full> *d = dynamic_cast<const DomainTypeT<full>*>(domain);
    m_stmtTDomains->setBoolean(1, d->position() == DomainTPosition::BeforeE);
    m_stmtTDomains->setBoolean(2, true);
    sql::ResultSet *res = m_stmtTDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeT>(res, d->position());
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getEDomains(const AbstractDomainType *domain)
{
    const DomainTypeE<full> *d = dynamic_cast<const DomainTypeE<full>*>(domain);
    m_stmtEDomains->setBoolean(1, true);
    sql::ResultSet *res = m_stmtEDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeE>(res);
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getTeDomains(const AbstractDomainType *domain)
{
    const DomainTypeTe<full> *d = dynamic_cast<const DomainTypeTe<full>*>(domain);
    m_stmtTeDomains->setBoolean(1, true);
    sql::ResultSet *res = m_stmtTeDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeTe>(res, d->circularizing());
    delete res;
    return ret;
}

template<template<bool> class D, typename... Args>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::populateInitialDomains(sql::ResultSet *res, const Args&... args)
{
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> domains(new std::vector<std::shared_ptr<Domain>>());
    while (res->next()) {
        D<true> *d = new D<true>(args...);
        d->setDomainId(res->getUInt("domain.id"));
        d->setPathway(res->getUInt("domain.pathway_id"))->setTaxId(res->getUInt("main.organism"));
        domains->emplace_back(d);
    }
    return domains;
}*/
