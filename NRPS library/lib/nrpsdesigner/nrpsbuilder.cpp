#include "nrpsbuilder.h"
#include "abstractdatabaseconnector.h"
#include "origin.h"

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
    size_t NRPSDESIGNER_EXPORT hash<pair<uint32_t, uint32_t>>::operator()(const pair<uint32_t, uint32_t> &s) const
{
    return hash<uint32_t>()(s.first) ^ (std::hash<uint32_t>()(s.second) << 1);
}

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

Nrps NrpsBuilder::build(const std::vector<Monomer> &nrp)
{
    AbstractDatabaseConnector *db = AbstractDatabaseConnector::getInstance();
    Nrps nrps(nrp);
    std::vector<Node*> graph;
    m_startn = new Node();
    m_endn = new Node();
    graph.push_back(m_startn);

    m_taxonCache.clear();
    m_weightCache.clear();

    std::vector<std::shared_ptr<DomainTypeT>> tdomainsc = db->getTDomains(DomainTPosition::BeforeC);
    std::vector<std::shared_ptr<DomainTypeT>> tdomainse = db->getTDomains(DomainTPosition::BeforeE);
    std::vector<std::shared_ptr<DomainTypeE>> edomains = db->getEDomains();
    std::vector<std::shared_ptr<DomainTypeTe>> tedomains = db->getTeDomains();

    std::vector<Node*> lastmodulenodes;
    Configuration lastmonomerconfiguration;

    for (auto iter = nrp.begin(); iter != nrp.end(); ++iter) {
        auto adomains = db->getADomains(*iter);
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
            graph.push_back(d);
        }

        if (iter != nrp.begin()) {
            std::shared_ptr<std::vector<Node*>> cdomainnodes(new std::vector<Node*>());
            auto cdomains = db->getCDomains(*iter, lastmonomerconfiguration);
            for (const auto &domain : cdomains) {
                Node *d = makeNode(domain);
                if (domain->substrate() == iter->id() && !amatchingdomainnodes->empty())
                    d->neighbors = amatchingdomainnodes;
                else if (!anotmatchingdomainnodes->empty())
                    d->neighbors = anotmatchingdomainnodes;
                else {
                    delete d;
                    d = nullptr;
                }
                if (d != nullptr) {
                    graph.push_back(d);
                    cdomainnodes->push_back(d);
                }
            }
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
                graph.push_back(t);
            }
        if (needtdomainse) {
            std::shared_ptr<std::vector<Node*>> edomainnodes(new std::vector<Node*>());
            for (const auto &edomain : edomains) {
                Node *e = makeNode(edomain);
                lastmodulenodes.push_back(e);
                edomainnodes->push_back(e);
                graph.push_back(e);
            }
            for (const auto &domain : tdomainse) {
                Node *t = makeNode(domain);
                t->neighbors = edomainnodes;
                tdomainnodese->push_back(t);
                graph.push_back(t);
            }
        }
        lastmonomerconfiguration = iter->configuration();
    }

    std::shared_ptr<std::vector<Node*>> goal(new std::vector<Node*>(1, m_endn));
    std::shared_ptr<std::vector<Node*>> tedomainnodes(new std::vector<Node*>());
    for (const auto &domain : tedomains) {
        Node *n = makeNode(domain);
        tedomainnodes->push_back(n);
        graph.push_back(n);
        n->neighbors = goal;
    }
    for (Node *n : lastmodulenodes)
        n->neighbors = tedomainnodes;

    graph.push_back(m_endn);

    std::priority_queue<dijkstra_weight, std::vector<dijkstra_weight>, std::greater<dijkstra_weight>> heap;
    std::unordered_map<Node*, Node*> parents;
    heap.push(dijkstra_weight(0, m_startn, m_startn));
    while (!heap.empty()) {
        dijkstra_weight n = heap.top();
        heap.pop();
        float weight = std::get<0>(n);
        Node *node = std::get<1>(n), *parent = std::get<2>(n);
        if (parents.count(node) == 0) {
            parents.emplace(node, parent);
            if (node == m_endn)
                break;
            if (node->neighbors.use_count() > 0)
                for (size_t i = 0; i < node->neighbors->size(); ++i) {
                    heap.push(std::make_tuple(weight + makeWeight(node, node->neighbors->at(i)), node->neighbors->at(i), node));
                }
        }
    }
    Node *cur = parents[m_endn];
    while (cur != m_startn) {
        nrps.push_back(cur->data);
        cur = parents[cur];
    }
    std::reverse(nrps.begin(), nrps.end());
    for (Node *n : graph)
        delete n;
    m_startn = nullptr;
    m_endn = nullptr;
    return nrps;
}

Node* NrpsBuilder::makeNode(std::shared_ptr<Domain> d)
{
    Node *n = new Node(d);
    if (m_taxonCache.count(d->origin()->taxId()) == 0) {
        m_taxonCache.emplace(d->origin()->taxId(), Taxon(d->origin()->taxId()));
    }
    return n;
}

float NrpsBuilder::makeWeight(Node *lhs, Node *rhs)
{
    if (lhs == m_startn || rhs == m_endn)
        return 0;
    uint32_t taxid1 = lhs->data->origin()->taxId(), taxid2 = rhs->data->origin()->taxId();
    std::pair<uint32_t, uint32_t> key(std::min(taxid1, taxid2), std::max(taxid1, taxid2));
    float weight;
    if (m_weightCache.count(key) == 0) {
        auto dist = m_taxonCache.at(taxid1) - m_taxonCache.at(taxid2);
        weight = dist[0] + dist[1];
        m_weightCache.emplace(key, weight);
    } else
        weight = m_weightCache[key];
    return weight;
}
