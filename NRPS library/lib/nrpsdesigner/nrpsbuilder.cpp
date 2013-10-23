#include "nrpsbuilder.h"
#include "abstractdatabaseconnector.h"
#include "origin.h"
#include "product.h"
#include "globals_internal.h"

#include <queue>
#include <algorithm>
#include <tuple>

using namespace nrps;

namespace nrps
{
class Node
{
public:
    Node() {}
    Node(const std::shared_ptr<Domain> &d)
    : data(d) {}

    std::shared_ptr<Domain> data;
    std::shared_ptr<std::vector<Node*>> neighbors;
    std::vector<float> weights;
};

typedef std::tuple<float, Node*, Node*> dijkstra_weight;
}

namespace std
{
template<>
struct greater<dijkstra_weight>
{
public:
    bool operator()(const dijkstra_weight &lhs, const dijkstra_weight &rhs) const
    {
        return std::greater<float>()(std::get<0>(lhs), std::get<0>(rhs));
    }
};
}

NrpsBuilder::NrpsBuilder()
: m_scaffoldNrp(nullptr), m_scaffoldNrps(nullptr)
{}

NrpsBuilder::NrpsBuilder(const std::vector<Monomer> *scaffoldNrp, const Nrps *scaffoldNrps)
: m_scaffoldNrp(scaffoldNrp), m_scaffoldNrps(scaffoldNrps)
{}

void NrpsBuilder::setScaffold(const std::vector<Monomer> *scaffoldNrp, const Nrps *scaffoldNrps)
{
    m_scaffoldNrp = scaffoldNrp;
    m_scaffoldNrps = scaffoldNrps;
}

Nrps NrpsBuilder::build(const std::vector<Monomer> &nrp, bool indTag) throw (NetworkError, NCBITaxonomyError, TaxonomyDumpError, DatabaseError)
{
    m_db = AbstractDatabaseConnector::getInstance();
    Nrps ret;
    m_startn = new Node();
    m_endn = new Node();
    m_graph.push_back(m_startn);

    m_taxonCache.clear();
    m_weightCache.clear();
    m_parents.clear();

    std::vector<std::shared_ptr<DomainTypeT>> tdomainsc = m_db->getTDomains(DomainTPosition::BeforeC);
    std::vector<std::shared_ptr<DomainTypeT>> tdomainse = m_db->getTDomains(DomainTPosition::BeforeE);
    std::vector<std::shared_ptr<DomainTypeE>> edomains = m_db->getEDomains();
    std::vector<std::shared_ptr<DomainTypeTe>> tedomains = m_db->getTeDomains();

    std::vector<Node*> lastmodulenodes;
    Configuration lastmonomerconfiguration;

    auto iter = nrp.begin();
    auto end = nrp.end();
    Nrps::const_iterator nrpsiter;
    if (m_scaffoldNrp != nullptr && m_scaffoldNrps != nullptr) {
        auto nrpiter = m_scaffoldNrp->begin();
        nrpsiter = m_scaffoldNrps->begin();
        lastmodulenodes.push_back(m_startn);
        for (; iter != nrp.end() && nrpiter != m_scaffoldNrp->end() && iter->id() == nrpiter->id(); ++iter, ++nrpiter) {
            for (int i = 0; i == 0 || (*nrpsiter)->type() != DomainType::C; ++nrpsiter, ++i) {
                std::shared_ptr<std::vector<Node*>> neighs(new std::vector<Node*>());
                Node *n = makeNode(*nrpsiter);
                neighs->push_back(n);
                lastmodulenodes[0]->neighbors = neighs;
                lastmodulenodes[0] = n;
                m_graph.push_back(n);
            }
        }
        end = iter + 1;
        ++nrpiter;
        if (end < nrp.end() && end->configuration() != nrpiter->configuration())
            ++end;
        for (int i = 0; i < end - iter; ++i)
            ++nrpsiter; // skip C domain of module
            while (nrpsiter != m_scaffoldNrps->end() && (*nrpsiter)->type() != DomainType::C)
                ++nrpsiter;
    }
    for (; iter != end; ++iter) {
        auto adomains = m_db->getADomains(*iter, false);
        bool needtdomainsc = false, needtdomainse = false;
        std::shared_ptr<std::vector<Node*>> amatchingdomainnodes(new std::vector<Node*>());
        std::shared_ptr<std::vector<Node*>> anotmatchingdomainnodes(new std::vector<Node*>());
        std::shared_ptr<std::vector<Node*>> tdomainnodesc(new std::vector<Node*>());
        std::shared_ptr<std::vector<Node*>> tdomainnodese(new std::vector<Node*>());

        for (const auto &domain : adomains) {
            Node *d = makeNode(domain);
            if (domain->substrate() == iter->id()) {
                d->neighbors = tdomainnodesc;
                needtdomainsc = true;
                amatchingdomainnodes->push_back(d);
            }
            else {
                d->neighbors = tdomainnodese;
                needtdomainse = true;
                anotmatchingdomainnodes->push_back(d);
            }
            m_graph.push_back(d);
        }

        if (iter != nrp.begin()) {
            std::shared_ptr<std::vector<Node*>> cdomainnodes = makeCDomains(*iter, amatchingdomainnodes, anotmatchingdomainnodes, lastmonomerconfiguration);
            for (Node *n : lastmodulenodes) {
                n->neighbors = cdomainnodes;
            }
        } else {
            m_startn->neighbors = std::shared_ptr<std::vector<Node*>>(new std::vector<Node*>());
            for (Node *n : *amatchingdomainnodes)
                m_startn->neighbors->push_back(n);
            for (Node *n : *anotmatchingdomainnodes)
                m_startn->neighbors->push_back(n);
        }
        lastmodulenodes.clear();
        if (needtdomainsc)
            for (const auto &domain : tdomainsc) {
                Node *t = makeNode(domain);
                lastmodulenodes.push_back(t);
                tdomainnodesc->push_back(t);
                m_graph.push_back(t);
            }
        if (needtdomainse) {
            std::shared_ptr<std::vector<Node*>> edomainnodes(new std::vector<Node*>());
            for (const auto &edomain : edomains) {
                Node *e = makeNode(edomain);
                lastmodulenodes.push_back(e);
                edomainnodes->push_back(e);
                m_graph.push_back(e);
            }
            for (const auto &domain : tdomainse) {
                Node *t = makeNode(domain);
                t->neighbors = edomainnodes;
                tdomainnodese->push_back(t);
                m_graph.push_back(t);
            }
        }
        lastmonomerconfiguration = iter->configuration();
    }

    if (m_scaffoldNrp != nullptr && m_scaffoldNrps != nullptr) {
        for (; nrpsiter != m_scaffoldNrps->end(); ++nrpsiter) {
            std::shared_ptr<std::vector<Node*>> neighs(new std::vector<Node*>());
            Node *n = makeNode(*nrpsiter);
            neighs->push_back(n);
            for (Node *ln : lastmodulenodes) {
                ln->neighbors = neighs;
            }
            lastmodulenodes.clear();
            lastmodulenodes.push_back(n);
            m_graph.push_back(n);
        }
    } else if (indTag) {
        Monomer glu = m_db->getMonomer("gln");
        auto aoxa = m_db->getADomains(glu, true);
        std::shared_ptr<std::vector<Node*>> tdomainnodes(new std::vector<Node*>());
        std::shared_ptr<std::vector<Node*>> aoxamatchingdomainnodes(new std::vector<Node*>());
        std::shared_ptr<std::vector<Node*>> aoxanotmatchingdomainnodes(new std::vector<Node*>());
        for (const auto &domain : aoxa) {
            Node *d = makeNode(domain);
            d->neighbors = tdomainnodes;
            if (domain->substrate() == glu.id())
                aoxamatchingdomainnodes->push_back(d);
            else
                aoxanotmatchingdomainnodes->push_back(d);
            m_graph.push_back(d);
        }
        std::shared_ptr<std::vector<Node*>> cdomainnodes = makeCDomains(glu, aoxamatchingdomainnodes, aoxanotmatchingdomainnodes, lastmonomerconfiguration);
        for (Node *n : lastmodulenodes) {
            n->neighbors = cdomainnodes;
        }
        lastmodulenodes.clear();
        for (const auto &domain : tdomainsc) {
            Node *t = makeNode(domain);
            lastmodulenodes.push_back(t);
            tdomainnodes->push_back(t);
            m_graph.push_back(t);
        }
        ret.setIndigoidineTagged(true);
    }
    std::shared_ptr<std::vector<Node*>> goal(new std::vector<Node*>(1, m_endn));
    if (m_scaffoldNrp == nullptr || m_scaffoldNrps == nullptr || end == nrp.end()) {
        std::shared_ptr<std::vector<Node*>> tedomainnodes(new std::vector<Node*>());
        for (const auto &domain : tedomains) {
            Node *n = makeNode(domain);
            tedomainnodes->push_back(n);
            m_graph.push_back(n);
            n->neighbors = goal;
        }
        for (Node *n : lastmodulenodes)
            n->neighbors = tedomainnodes;
    } else {
        for (Node *n : lastmodulenodes)
            n->neighbors = goal;
    }
    m_graph.push_back(m_endn);
    TaxonBuilder::getInstance()->process();

    std::priority_queue<dijkstra_weight, std::vector<dijkstra_weight>, std::greater<dijkstra_weight>> heap;
    heap.push(dijkstra_weight(0, m_startn, m_startn));
    while (!heap.empty()) {
        dijkstra_weight n = heap.top();
        heap.pop();
        float weight = std::get<0>(n);
        Node *node = std::get<1>(n), *parent = std::get<2>(n);
        if (!m_parents.count(node)) {
            m_parents.emplace(node, parent);
            if (node == m_endn)
                break;
            if (node->neighbors.use_count())
                for (size_t i = 0; i < node->neighbors->size(); ++i) {
                    heap.push(std::make_tuple(weight + makeWeight(node, node->neighbors->at(i)), node->neighbors->at(i), node));
                }
        }
    }
    Node *cur = m_parents[m_endn];
    while (cur != m_startn) {
        m_db->fillDomain(cur->data);
        ret.push_back(cur->data);
        cur->data->setDeterminedLinkerBefore(cur->data->nativeDefinedLinkerBefore().empty() ? cur->data->nativePfamLinkerBefore() : cur->data->nativeDefinedLinkerBefore());
        cur->data->setDeterminedLinkerAfter(std::string());
        cur = m_parents[cur];
    }
    std::reverse(ret.begin(), ret.end());
    for (Node *n : m_graph)
        delete n;
    m_startn = nullptr;
    m_endn = nullptr;
    m_parents.clear();
    m_graph.clear();
    return ret;
}

std::shared_ptr<std::vector<Node*>> NrpsBuilder::makeCDomains(const Monomer &m, std::shared_ptr<std::vector<Node*>> amatchingdomainnodes, std::shared_ptr<std::vector<Node*>> anotmatchingdomainnodes, Configuration lastmonomerconfiguration)
{
    std::shared_ptr<std::vector<Node*>> cdomainnodes(new std::vector<Node*>());
    auto cdomains = m_db->getCDomains(m, lastmonomerconfiguration);
    for (const auto &domain : cdomains) {
        Node *d = makeNode(domain);
        if (domain->substrate() == m.id() && !amatchingdomainnodes->empty())
            d->neighbors = amatchingdomainnodes;
        else if (!anotmatchingdomainnodes->empty())
            d->neighbors = anotmatchingdomainnodes;
        else if (m_db->isDummy(domain)) { // only domain in list
            if (!amatchingdomainnodes->empty() && anotmatchingdomainnodes->empty())
                d->neighbors = amatchingdomainnodes;
            else if (amatchingdomainnodes->empty() && !anotmatchingdomainnodes->empty())
                d->neighbors = anotmatchingdomainnodes;
            else {
                amatchingdomainnodes->insert(amatchingdomainnodes->end(), anotmatchingdomainnodes->begin(), anotmatchingdomainnodes->end());
                d->neighbors = amatchingdomainnodes;
            }
        } else {
            delete d;
            d = nullptr;
        }
        if (d != nullptr) {
            m_graph.push_back(d);
            cdomainnodes->push_back(d);
        }
    }
    return cdomainnodes;
}

Node* NrpsBuilder::makeNode(const std::shared_ptr<Domain> &d)
{
    Node *n = new Node(d);
    if (d->origin() != nullptr && !m_taxonCache.count(d->origin()->taxId())) {
        m_taxonCache.emplace(d->origin()->taxId(), TaxonBuilder::getInstance()->buildMany(d->origin()->taxId()));
    }
    return n;
}

float NrpsBuilder::makeWeight(Node *lhs, Node *rhs)
{
    if (lhs == m_startn || rhs == m_endn || rhs->data->origin() == nullptr)
        return 0;
    if (lhs->data->origin() == nullptr && rhs->data->origin() != nullptr) {
        Node *n = lhs;
        while (n != m_startn && n->data->origin() == nullptr)
            n = m_parents[n];
        if (n == m_startn)
            return 0;
        lhs = n;
    }
    return nrps::makeWeight(lhs->data, rhs->data, &m_taxonCache, &m_weightCache);
}
