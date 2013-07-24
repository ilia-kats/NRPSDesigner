#include "mysqldatabaseconnector.h"
#include "domaintypea.h"
#include "domaintypeac.h"
#include "domaintypee.h"
#include "domaintypet.h"
#include "domaintypete.h"

#include <unordered_map>
#include <stdexcept>

#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/prepared_statement.h>

using namespace nrps;

MySQLDatabaseConnector::MySQLDatabaseConnector()
: AbstractDatabaseConnector(), m_connection(nullptr), m_stmtInitiationADomains(nullptr), m_stmtElongationACDomains(nullptr), m_stmtTDomains(nullptr), m_stmtEDomains(nullptr), m_stmtTeDomains(nullptr)
{}

MySQLDatabaseConnector::~MySQLDatabaseConnector()
{
    if (m_connection != nullptr) {
        m_connection->close();
        delete m_connection;
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
    sql::Driver *driver = get_driver_instance();
    m_connection = driver->connect("tcp://127.0.0.1:3306", "root", "");
    m_connection->setSchema("nrps_designer");
    m_stmtInitiationADomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_a` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE domain.substrate_specificity_aa_id = ? AND (? OR domain.bb_id != '');");
    m_stmtElongationACDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_a_c` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE domain.substrate_specificity_aa_id = ? AND domain.chirality = ? AND (? OR domain.bb_id != '');");
    m_stmtTDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_t` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE domain.positioning = ? AND (? OR domain.bb_id != '');");
    m_stmtEDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_e` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE ? OR domain.bb_id != '';");
    m_stmtTeDomains = m_connection->prepareStatement("SELECT domain.domain_id, domain.pathway_id, main.organism FROM `domain_types_te` AS domain INNER JOIN `main` on domain.pathway_id = main.pathway_id WHERE ? OR domain.bb_id != '';");
}

std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> MySQLDatabaseConnector::getPotentialDomains(const GeneralPathway &pathway)
{
    if (m_connection == nullptr || m_stmtInitiationADomains == nullptr || m_stmtElongationACDomains == nullptr || m_stmtTDomains == nullptr || m_stmtEDomains == nullptr || m_stmtTeDomains == nullptr)
        throw std::logic_error("Not inizialized");
    std::unordered_map<AbstractDomainType*, std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> cache;
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> potentialDomains;
    for (const std::shared_ptr<AbstractDomainType> &domain : pathway) {
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
    }
    return potentialDomains;
}

template <bool full>
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
        d->setPathway(Pathway(res->getUInt("domain.pathway_id"), res->getUInt("main.organism")));
        domains->emplace_back(d);
    }
    return domains;
}
