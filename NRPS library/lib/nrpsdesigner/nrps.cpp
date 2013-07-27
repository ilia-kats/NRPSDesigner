#include "nrps.h"
#include "taxon.h"

#include <unordered_map>
#include <queue>
#include <algorithm>

using namespace nrps;

class Node
{
public:
    Node() {}
    Node(const std::shared_ptr<Domain> &d)
    : data(d) {}

    std::shared_ptr<Domain> data;
    std::shared_ptr<std::vector<Node*>> neighbors;
    std::vector<uint8_t> weights;
};

typedef std::tuple<uint8_t, Node*, Node*> dijkstra_weight;

namespace std
{
    template<>
    struct hash<std::pair<uint32_t, uint32_t>>
    {
    public:
        std::size_t operator()(const std::pair<uint32_t, uint32_t> &s) const
        {
            return std::hash<uint32_t>()(s.first) ^ (std::hash<uint32_t>()(s.second) << 1);
        }
    };

    template<>
    struct greater<dijkstra_weight>
    {
    public:
        bool operator()(const dijkstra_weight &lhs, const dijkstra_weight &rhs) const
        {
            return std::greater<uint8_t>()(std::get<0>(lhs), std::get<0>(rhs));
        }
    };
}

Nrps::Nrps(const std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>> &potentialPathway)
: std::vector<std::shared_ptr<Domain>>()
{
    std::vector<std::shared_ptr<std::vector<Node*>>> graph;
    Node *startn = new Node(), *endn = new Node();
    graph.emplace_back(new std::vector<Node*>(1, startn));

    std::unordered_map<uint32_t, Taxon> taxon_cache;
    std::unordered_map<std::pair<uint32_t, uint32_t>, uint8_t> weight_cache;

    for (const auto &position : potentialPathway) {
        std::shared_ptr<std::vector<Node*>> nodes(new std::vector<Node*>());
        for (const auto &domain : *position) {
            nodes->push_back(new Node(domain));
            uint32_t taxid = domain->pathway().taxId();
            if (taxon_cache.count(taxid) == 0)
                taxon_cache.emplace(taxid, Taxon(taxid));
        }
        auto lastNodes = graph.back();
        for (auto iter = lastNodes->begin(); iter != lastNodes->end(); ++iter) {
            (*iter)->neighbors = nodes;
            uint32_t taxid1 = (*iter)->data->pathway().taxId();
            for (Node *node : *nodes) {
                uint32_t taxid2 = node->data->pathway().taxId();
                std::pair<uint32_t, uint32_t> key(std::min(taxid1, taxid2), std::max(taxid1, taxid2));
                uint8_t weight;
                if (weight_cache.count(key) == 0) {
                    auto dist = taxon_cache.at(taxid1) - taxon_cache.at(taxid2);
                    weight = dist[0] + dist[1];
                    weight_cache.emplace(key, weight);
                } else
                    weight = weight_cache[key];
                (*iter)->weights.push_back(weight);
            }
        }
        graph.emplace_back(std::move(nodes));
    }
    std::shared_ptr<std::vector<Node*>> goal(new std::vector<Node*>(1, endn));
    for (auto iter = graph.back()->begin(); iter != graph.back()->end(); ++iter) {
        (*iter)->neighbors = goal;
    }
    graph.emplace_back(std::move(goal));

    std::priority_queue<dijkstra_weight, std::vector<dijkstra_weight>, std::greater<dijkstra_weight>> heap;
    std::unordered_map<Node*, Node*> parents;
    heap.push(dijkstra_weight(0, startn, startn));
    while (!heap.empty()) {
        dijkstra_weight n = heap.top();
        heap.pop();
        uint8_t weight = std::get<0>(n);
        Node *node = std::get<1>(n), *parent = std::get<2>(n);
        if (parents.count(node) == 0) {
            parents.emplace(node, parent);
            if (node == endn)
                break;
            for (size_t i = 0; i < node->neighbors->size(); ++i) {
                heap.push(std::make_tuple(weight + node->weights[i], node->neighbors->at(i), node));
            }
        }
    }
    Node *cur = parents[endn];
    while (cur != startn) {
        push_back(cur->data);
        cur = parents[cur];
    }
    std::reverse(begin(), end());
    for (auto &pos : graph) {
        for (Node *n : *pos)
            delete n;
    }
}
