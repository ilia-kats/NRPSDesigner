#ifndef NRPSDESIGNER_NRPSBUILDER_H
#define NRPSDESIGNER_NRPSBUILDER_H

#include "nrpsdesigner_export.h"
#include "monomer.h"
#include "taxon.h"
#include "domain.h"
#include "nrps.h"

#include <unordered_map>

namespace std
{
template<>
struct NRPSDESIGNER_EXPORT hash<pair<uint32_t, uint32_t>>
{
public:
    size_t operator()(const pair<uint32_t, uint32_t>&) const;
};
}

namespace nrps
{
class Node;
class NRPSDESIGNER_EXPORT NrpsBuilder
{
public:
    Nrps build(const std::vector<Monomer>&);

private:
    float makeWeight(Node*, Node*);
    Node* makeNode(std::shared_ptr<Domain>);

    std::unordered_map<uint32_t, Taxon> m_taxonCache;
    std::unordered_map<std::pair<uint32_t, uint32_t>, float> m_weightCache;
    Node *m_startn;
    Node *m_endn;
};
}

#endif
