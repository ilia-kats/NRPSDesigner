#ifndef NRPSDESIGNER_MYSQLDATABASECONNECTOR_H
#define NRPSDESIGNER_MYSQLDATABASECONNECTOR_H

#include "nrpsdesigner_export.h"
#include "abstractdatabaseconnector.h"

namespace sql
{
class Connection;
class PreparedStatement;
class ResultSet;
}

namespace nrps
{
template<bool>
class DomainTypeA;
template<bool>
class DomainTypeAC;
template<bool>
class DomainTypeT;
template<bool>
class DomainTypeE;
template<bool>
class DomainTypeTe;
class Domain;

class MySQLDatabaseConnector : public AbstractDatabaseConnector
{
public:
    MySQLDatabaseConnector();
    virtual ~MySQLDatabaseConnector();
    virtual void initialize();
    virtual std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> getPotentialDomains(const GeneralPathway &pathway);

private:
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getADomains(const DomainTypeA<full>*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getACDomains(const DomainTypeAC<full>*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getTDomains(const DomainTypeT<full>*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getEDomains(const DomainTypeE<full>*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getTeDomains(const DomainTypeTe<full>*);
    template<template<bool> class D, typename... Args>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> populateInitialDomains(sql::ResultSet *res, const Args&... args);
    sql::Connection *m_connection;
    sql::PreparedStatement *m_stmtInitiationADomains;
    sql::PreparedStatement *m_stmtElongationACDomains;
    sql::PreparedStatement *m_stmtTDomains;
    sql::PreparedStatement *m_stmtEDomains;
    sql::PreparedStatement *m_stmtTeDomains;
};
}

#endif
