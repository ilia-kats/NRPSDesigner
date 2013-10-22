#ifndef NRPSDESIGNER_NRPSBUILDER_H
#define NRPSDESIGNER_NRPSBUILDER_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"
#include "monomer.h"
#include "taxon.h"
#include "domain.h"
#include "nrps.h"
#include "exceptions.h"
#include "taxonbuilder.h"
#include "globals_internal.h"

#include <unordered_map>
#include <memory>

namespace nrps
{
class Node;
class AbstractDatabaseConnector;
class NRPSDESIGNER_EXPORT NrpsBuilder
{
public:
    NrpsBuilder();
    NrpsBuilder(const std::vector<Monomer>*, const Nrps*);
    void setScaffold(const std::vector<Monomer>*, const Nrps*);
    Nrps build(const std::vector<Monomer>&, bool indTag = false) throw (NetworkError, NCBITaxonomyError, TaxonomyDumpError, DatabaseError);

private:
    std::shared_ptr<std::vector<Node*>> makeCDomains(const Monomer&, std::shared_ptr<std::vector<Node*>>, std::shared_ptr<std::vector<Node*>>, Configuration);
    float makeWeight(Node*, Node*);
    Node* makeNode(const std::shared_ptr<Domain>&);

    const std::vector<Monomer> *m_scaffoldNrp;
    const Nrps *m_scaffoldNrps;
    std::unordered_map<uint32_t, std::shared_ptr<Taxon>> m_taxonCache;
    std::unordered_map<std::pair<uint32_t, uint32_t>, float> m_weightCache;
    std::vector<Node*> m_graph;
    AbstractDatabaseConnector *m_db;
    Node *m_startn;
    Node *m_endn;
    std::unordered_map<Node*, Node*> m_parents;
};
}

#endif
