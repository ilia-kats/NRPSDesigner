#include "generalpathway.h"
#include "domaintypea.h"
#include "domaintypeac.h"
#include "domaintypet.h"
#include "domaintypee.h"
#include "domaintypete.h"

using namespace nrps;

GeneralPathway::GeneralPathway(const Nrp &nrp)
: std::vector<std::shared_ptr<AbstractDomainType>>(), m_lastConfiguration(Configuration::L)
{
    auto iter = nrp.cbegin(), end = nrp.cend();
    processMonomer(*iter, true);
    for (++iter; iter != end; ++iter) {
        processMonomer(*iter);
    }
    emplace_back(new DomainTypeTe<>(nrp.type() == Nrp::Type::Circular));
}

GeneralPathway::~GeneralPathway()
{}

void GeneralPathway::processMonomer(const Monomer &monomer, bool initial)
{
    emplace_back(initial ? new DomainTypeA<>(monomer.id()) : new DomainTypeAC<>(monomer.id(), m_lastConfiguration));
    m_lastConfiguration = monomer.configuration();

    // tailoring before T
    /*if (monomer.hasModification(Monomer::Modification::Nmethyl)) {
        // add NMT-Domain
    }*/

    emplace_back(new DomainTypeT<>(m_lastConfiguration == Configuration::L ? DomainTPosition::BeforeC : DomainTPosition::BeforeE));

    // add tailoring domains

    if (m_lastConfiguration == Configuration::D)
        emplace_back(new DomainTypeE<>());
}
