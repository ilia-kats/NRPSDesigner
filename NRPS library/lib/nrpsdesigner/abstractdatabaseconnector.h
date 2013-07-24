#ifndef NRPSDESIGNER_ABSTRACTDATABASECONNECTOR_H
#define NRPSDESIGNER_ABSTRACTDATABASECONNECTOR_H

#include "nrpsdesigner_export.h"
#include "domain.h"
#include "generalpathway.h"

#include <vector>

namespace nrps
{
class NRPSDESIGNER_EXPORT AbstractDatabaseConnector
{
public:
    virtual ~AbstractDatabaseConnector();
    virtual void initialize() = 0; // TODO: write class for parameters (host, user, pw)
    virtual std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> getPotentialDomains(const GeneralPathway &pathway) = 0;

    static AbstractDatabaseConnector* getInstance();

private:
    static AbstractDatabaseConnector *s_instance;
};
}

#endif
