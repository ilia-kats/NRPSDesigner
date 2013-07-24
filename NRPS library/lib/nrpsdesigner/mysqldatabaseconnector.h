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
class AbstractDomainType;
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
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getADomains(const AbstractDomainType*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getACDomains(const AbstractDomainType*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getTDomains(const AbstractDomainType*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getEDomains(const AbstractDomainType*);
    template <bool full>
    std::shared_ptr<std::vector<std::shared_ptr<Domain>>> getTeDomains(const AbstractDomainType*);
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
