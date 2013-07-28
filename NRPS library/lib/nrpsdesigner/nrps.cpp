#include "nrps.h"
#include "taxon.h"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>

#include <unistd.h>

#define NRPS_NODE "nrps"
#define DOMAIN_NODE "domain"
#define TYPE_NODE "type"
#define ID_NODE "id"
#define MODULEID_NODE "moduleid"
#define BIOBRICK_NODE "biobrick"
#define DESCRIPTION_NODE "description"
#define SEQUENCE_NODE "sequence"
#define NATIVELINKERBEFORE_NODE "nativelinkerbefore"
#define NATIVELINKERAFTER_NODE "nativelinkerafter"
#define REFSEQID_NODE "refseqid"
#define UNIPROTID_NODE "uniprotid"
#define PATHWAYID_NODE "pathwayid"
#define TAXID_NODE "taxid"
#define PATHWAY_NODE "pathway"
#define LINKOUT_NODE "linkout"
#define PATHWAYNAME_NODE "product"
#define NORINEID_NODE "norineid"
#define DOMAINS_NODE "domains"
#define PATHWAYS_NODE "pathways"

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
    std::unordered_set<std::shared_ptr<Pathway>> seenPathways;
    xmlTextWriterSetIndent(writer, 1);
    xmlTextWriterSetIndentString(writer, BAD_CAST "    ");
    xmlTextWriterStartDocument(writer, nullptr, "UTF-8", nullptr);
    xmlTextWriterStartElement(writer, BAD_CAST NRPS_NODE);
    xmlTextWriterStartElement(writer, BAD_CAST DOMAINS_NODE);
    for (const auto &domain : *this) {
        seenPathways.insert(domain->pathway());
        xmlTextWriterStartElement(writer, BAD_CAST DOMAIN_NODE);
        xmlTextWriterWriteElement(writer, BAD_CAST ID_NODE, BAD_CAST std::to_string(domain->domainId()).c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST MODULEID_NODE, BAD_CAST std::to_string(domain->moduleId()).c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST PATHWAYID_NODE, BAD_CAST std::to_string(domain->pathway()->pathwayId()).c_str());
        if (!domain->bioBrickId().empty())
            xmlTextWriterWriteElement(writer, BAD_CAST BIOBRICK_NODE, BAD_CAST domain->bioBrickId().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST DESCRIPTION_NODE, BAD_CAST domain->description().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST SEQUENCE_NODE, BAD_CAST domain->dnaSequence().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST NATIVELINKERBEFORE_NODE, BAD_CAST domain->nativeLinkerBefore().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST NATIVELINKERAFTER_NODE, BAD_CAST domain->nativeLinkerAfter().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST REFSEQID_NODE, BAD_CAST domain->refSeqId().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST UNIPROTID_NODE, BAD_CAST domain->uniProtId().c_str());
        xmlTextWriterEndElement(writer);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterStartElement(writer, BAD_CAST PATHWAYS_NODE);
    for (const auto &pathway : seenPathways) {
        xmlTextWriterStartElement(writer, BAD_CAST PATHWAY_NODE);
        xmlTextWriterWriteElement(writer, BAD_CAST PATHWAYID_NODE, BAD_CAST std::to_string(pathway->pathwayId()).c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST TAXID_NODE, BAD_CAST std::to_string(pathway->taxId()).c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST PATHWAYNAME_NODE, BAD_CAST pathway->pathway().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST LINKOUT_NODE, BAD_CAST pathway->linkout().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST UNIPROTID_NODE, BAD_CAST pathway->uniProtId().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST NORINEID_NODE, BAD_CAST pathway->norineId().c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST DESCRIPTION_NODE, BAD_CAST pathway->description().c_str());
        xmlTextWriterEndElement(writer);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterEndDocument(writer);
}
