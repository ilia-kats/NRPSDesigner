#include "nrps.h"
#include "taxon.h"
#include "origin.h"
#include "abstractdatabaseconnector.h"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <tuple>

#include <unistd.h>

#define NRPS_NODE "nrps"
#define DOMAINS_NODE "domains"
#define ORIGINS_NODE "origins"

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
            return std::greater<float>()(std::get<0>(lhs), std::get<0>(rhs));
        }
    };
}

Nrps::Nrps(const std::vector<Monomer> &nrp)
: std::vector<std::shared_ptr<Domain>>(), m_nrp(nrp)
{
    AbstractDatabaseConnector *db = AbstractDatabaseConnector::getInstance();
    db->initialize();

    std::vector<Node*> graph;
    Node *startn = new Node(), *endn = new Node();
    graph.push_back(startn);

    std::unordered_map<uint32_t, Taxon> taxon_cache;
    std::unordered_map<std::pair<uint32_t, uint32_t>, float> weight_cache;

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
            Node *d = new Node(domain);
            if (taxon_cache.count(domain->origin()->taxId()) == 0) {
                taxon_cache.emplace(domain->origin()->taxId(), Taxon(domain->origin()->taxId()));
            }
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
                Node *d = new Node(domain);
                if (taxon_cache.count(domain->origin()->taxId()) == 0) {
                    taxon_cache.emplace(domain->origin()->taxId(), Taxon(domain->origin()->taxId()));
                }
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
            startn->neighbors = std::shared_ptr<std::vector<Node*>>(new std::vector<Node*>());
            for (Node *n : *amatchingdomainnodes)
                startn->neighbors->push_back(n);
            for (Node *n : *anotmatchingdomainnodes)
                startn->neighbors->push_back(n);
        }
        lastmodulenodes.clear();
        if (needtdomainsc)
            for (const auto &domain : tdomainsc) {
                Node *t = new Node(domain);
                if (taxon_cache.count(domain->origin()->taxId()) == 0) {
                    taxon_cache.emplace(domain->origin()->taxId(), Taxon(domain->origin()->taxId()));
                }
                lastmodulenodes.push_back(t);
                tdomainnodesc->push_back(t);
                graph.push_back(t);
            }
        if (needtdomainse) {
            std::shared_ptr<std::vector<Node*>> edomainnodes(new std::vector<Node*>());
            for (const auto &edomain : edomains) {
                Node *e = new Node(edomain);
                if (taxon_cache.count(edomain->origin()->taxId()) == 0) {
                    taxon_cache.emplace(edomain->origin()->taxId(), Taxon(edomain->origin()->taxId()));
                }
                lastmodulenodes.push_back(e);
                edomainnodes->push_back(e);
                graph.push_back(e);
            }
            for (const auto &domain : tdomainse) {
                Node *t = new Node(domain);
                if (taxon_cache.count(domain->origin()->taxId()) == 0) {
                    taxon_cache.emplace(domain->origin()->taxId(), Taxon(domain->origin()->taxId()));
                }
                t->neighbors = edomainnodes;
                tdomainnodese->push_back(t);
                graph.push_back(t);
            }
        }
        lastmonomerconfiguration = iter->configuration();
    }

    std::shared_ptr<std::vector<Node*>> goal(new std::vector<Node*>(1, endn));
    std::shared_ptr<std::vector<Node*>> tedomainnodes(new std::vector<Node*>());
    for (const auto &domain : tedomains) {
        Node *n = new Node(domain);
        if (taxon_cache.count(domain->origin()->taxId()) == 0) {
            taxon_cache.emplace(domain->origin()->taxId(), Taxon(domain->origin()->taxId()));
        }
        tedomainnodes->push_back(n);
        graph.push_back(n);
        n->neighbors = goal;
    }
    for (Node *n : lastmodulenodes)
        n->neighbors = tedomainnodes;

    /*for (const auto &position : potentialPathway) {
        std::shared_ptr<std::vector<Node*>> nodes(new std::vector<Node*>());
        for (const auto &domain : *position) {
            nodes->push_back(new Node(domain));
            uint32_t taxid = domain->pathway()->taxId();
            if (taxon_cache.count(taxid) == 0)
                taxon_cache.emplace(taxid, Taxon(taxid));
        }
        auto lastNodes = graph.back();
        for (auto iter = lastNodes->begin(); iter != lastNodes->end(); ++iter) {
            (*iter)->neighbors = nodes;
            uint32_t taxid1 = (*iter)->data->pathway()->taxId();
            for (Node *node : *nodes) {
                uint32_t taxid2 = node->data->pathway()->taxId();
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
    }*/
    graph.push_back(endn);

    std::priority_queue<dijkstra_weight, std::vector<dijkstra_weight>, std::greater<dijkstra_weight>> heap;
    std::unordered_map<Node*, Node*> parents;
    heap.push(dijkstra_weight(0, startn, startn));
    while (!heap.empty()) {
        dijkstra_weight n = heap.top();
        heap.pop();
        float weight = std::get<0>(n);
        Node *node = std::get<1>(n), *parent = std::get<2>(n);
        if (parents.count(node) == 0) {
            parents.emplace(node, parent);
            if (node == endn)
                break;
            if (node->neighbors.use_count() > 0)
                for (size_t i = 0; i < node->neighbors->size(); ++i) {
                    float newweight;
                    if (node == startn || node->neighbors->at(i) == endn)
                        newweight = 0;
                    else {
                        uint32_t taxid1 = node->data->origin()->taxId();
                        uint32_t taxid2 = node->neighbors->at(i)->data->origin()->taxId();
                        std::pair<uint32_t, uint32_t> key(std::min(taxid1, taxid2), std::max(taxid1, taxid2));
                        if (weight_cache.count(key) == 0) {
                            newweight = makeWeight(taxon_cache.at(taxid1), taxon_cache.at(taxid2));
                            weight_cache.emplace(key, newweight);
                        } else
                            newweight = weight_cache[key];
                    }
                    heap.push(std::make_tuple(weight + newweight, node->neighbors->at(i), node));
                }
        }
    }
    Node *cur = parents[endn];
    while (cur != startn) {
        push_back(cur->data);
        cur = parents[cur];
    }
    std::reverse(begin(), end());
    for (Node *n : graph)
        delete n;
}

std::string Nrps::toXml() const
{
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    toXml(writer);
    std::string xml((const char*)buf->content);
    xmlFreeTextWriter(writer);
    xmlBufferFree(buf);
    return xml;
}

void Nrps::toXml(std::ostream &of) const
{
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    toXml(writer);
    of << (const char*)buf->content;
    xmlFreeTextWriter(writer);
    xmlBufferFree(buf);
}

void Nrps::toXml(const std::string &file) const
{
    toXml(file.c_str());
}

void Nrps::toXml(const char *file) const
{
    xmlTextWriterPtr writer = xmlNewTextWriterFilename(file, 0);
    toXml(writer);
    xmlFreeTextWriter(writer);
}

void Nrps::toXml(int fd) const
{
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    toXml(writer);
    write(fd, buf->content, buf->use);
    xmlFreeTextWriter(writer);
    xmlBufferFree(buf);
}

void Nrps::toXml(xmlTextWriterPtr writer) const
{
    std::unordered_set<Origin*> seenOrigins;
    std::vector<Origin*> originsToWrite;
    xmlTextWriterSetIndent(writer, 1);
    xmlTextWriterSetIndentString(writer, BAD_CAST "    ");
    xmlTextWriterStartDocument(writer, nullptr, "UTF-8", nullptr);
    xmlTextWriterStartElement(writer, BAD_CAST NRPS_NODE);
    xmlTextWriterStartElement(writer, BAD_CAST DOMAINS_NODE);
    AbstractDatabaseConnector *dbconn = AbstractDatabaseConnector::getInstance();
    for (const auto &domain : *this) {
        dbconn->fillDomain(domain);
        originsToWrite.push_back(domain->origin());
        domain->toXml(writer);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterStartElement(writer, BAD_CAST ORIGINS_NODE);
    for (Origin *ori : originsToWrite) {
        if (seenOrigins.count(ori) == 0) {
            while (ori != nullptr && seenOrigins.count(ori) == 0) {
                ori->toXml(writer);
                seenOrigins.insert(ori);
                ori = ori->parent();
            }
        }
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterEndDocument(writer);
}

float Nrps::makeWeight(const Taxon &lhs, const Taxon &rhs)
{
    auto dist = lhs - rhs;
    float weight = dist[0] + dist[1];
    return weight;
}
