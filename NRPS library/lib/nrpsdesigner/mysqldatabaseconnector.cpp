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
    s_instance = nullptr;
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
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> sentinel;
    for (const std::shared_ptr<AbstractDomainType> &domain : pathway) {
        if (cache.count(domain.get()) > 0)
            potentialDomains.push_back(cache[domain.get()]);
        else {
            switch (domain->type()) {
                case DomainType::A:
                    if (domain->full()) {
                        const DomainTypeA<true> *d = dynamic_cast<const DomainTypeA<true>*>(domain.get());
                        potentialDomains.push_back(getADomains(d));
                    } else {
                        const DomainTypeA<false> *d = dynamic_cast<const DomainTypeA<false>*>(domain.get());
                        potentialDomains.push_back(getADomains(d));
                    }
                    break;
                case DomainType::AC:
                    if (domain->full()) {
                        const DomainTypeAC<true> *d = dynamic_cast<const DomainTypeAC<true>*>(domain.get());
                        potentialDomains.push_back(getACDomains(d));
                    } else {
                        const DomainTypeAC<false> *d = dynamic_cast<const DomainTypeAC<false>*>(domain.get());
                        potentialDomains.push_back(getACDomains(d));
                    }
                    break;
                case DomainType::T:
                    if (domain->full()) {
                        const DomainTypeT<true> *d = dynamic_cast<const DomainTypeT<true>*>(domain.get());
                        potentialDomains.push_back(getTDomains(d));
                    } else {
                        const DomainTypeT<false> *d = dynamic_cast<const DomainTypeT<false>*>(domain.get());
                        potentialDomains.push_back(getTDomains(d));
                    }
                    break;
                case DomainType::E:
                    if (domain->full()) {
                        const DomainTypeE<true> *d = dynamic_cast<const DomainTypeE<true>*>(domain.get());
                        potentialDomains.push_back(getEDomains(d));
                    } else {
                        const DomainTypeE<false> *d = dynamic_cast<const DomainTypeE<false>*>(domain.get());
                        potentialDomains.push_back(getEDomains(d));
                    }
                    break;
                case DomainType::Te:
                    if (domain->full()) {
                        const DomainTypeTe<true> *d = dynamic_cast<const DomainTypeTe<true>*>(domain.get());
                        potentialDomains.push_back(getTeDomains(d));
                    } else {
                        const DomainTypeTe<false> *d = dynamic_cast<const DomainTypeTe<false>*>(domain.get());
                        potentialDomains.push_back(getTeDomains(d));
                    }
                    break;
                default:
                    potentialDomains.push_back(sentinel);
                    break;
            }
        }
    }
    return potentialDomains;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getADomains(const DomainTypeA<full>* d)
{
    m_stmtInitiationADomains->setUInt(1, d->substrate());
    m_stmtInitiationADomains->setBoolean(2, true);
    sql::ResultSet *res = m_stmtInitiationADomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeA>(res, d->substrate());
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getACDomains(const DomainTypeAC<full>* d)
{
    m_stmtElongationACDomains->setUInt(1, d->substrate());
    m_stmtElongationACDomains->setBoolean(2, d->chirality() == Configuration::D);
    m_stmtElongationACDomains->setBoolean(3, true);
    sql::ResultSet *res = m_stmtElongationACDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeAC>(res, d->substrate(), d->chirality());
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getTDomains(const DomainTypeT<full>* d)
{
    m_stmtTDomains->setBoolean(1, d->position() == DomainTPosition::BeforeE);
    m_stmtTDomains->setBoolean(2, true);
    sql::ResultSet *res = m_stmtTDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeT>(res, d->position());
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getEDomains(const DomainTypeE<full>* d)
{
    m_stmtEDomains->setBoolean(1, true);
    sql::ResultSet *res = m_stmtEDomains->executeQuery();
    auto ret = populateInitialDomains<DomainTypeE>(res);
    delete res;
    return ret;
}

template <bool full>
std::shared_ptr<std::vector<std::shared_ptr<Domain>>> MySQLDatabaseConnector::getTeDomains(const DomainTypeTe<full>* d)
{
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
