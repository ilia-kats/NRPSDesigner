#ifndef NRPSDESIGNER_GENERALPATHWAY_H
#define NRPSDESIGNER_GENERALPATHWAY_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"
#include "monomer.h"
#include "nrp.h"
#include "global_enums.h"

#include <memory>
#include <vector>

namespace nrps
{
class NRPSDESIGNER_EXPORT GeneralPathway : public std::vector<std::shared_ptr<AbstractDomainType>>
{
public:
    GeneralPathway(const Nrp&);
    ~GeneralPathway();

private:
    void processMonomer(const Monomer&, bool initial=false);

    Configuration m_lastConfiguration;
};
}

#endif
