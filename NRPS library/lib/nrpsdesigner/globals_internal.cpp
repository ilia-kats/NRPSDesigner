#include "globals_internal.h"
#include "origin.h"
#include "product.h"
#include "taxonbuilder.h"

namespace std
{
size_t hash<pair<uint32_t, uint32_t>>::operator()(const pair<uint32_t, uint32_t> &s) const
{
    return hash<uint32_t>()(s.first) ^ (std::hash<uint32_t>()(s.second) << 1);
}
}

float nrps::makeWeight(const std::shared_ptr<Domain> &lhs, const std::shared_ptr<Domain> &rhs, std::unordered_map<uint32_t, std::shared_ptr<Taxon>> *taxonCache,  std::unordered_map<std::pair<uint32_t, uint32_t>, float> *weightCache)
{
    uint32_t taxid1 = lhs->origin()->taxId(), taxid2 = rhs->origin()->taxId();
    std::pair<uint32_t, uint32_t> key(std::min(taxid1, taxid2), std::max(taxid1, taxid2));
    float weight;
    if (weightCache == nullptr || !weightCache->count(key)) {
        if (!taxonCache->count(taxid1))
            taxonCache->emplace(taxid1, TaxonBuilder::getInstance()->buildOne(taxid1));
        if (!taxonCache->count(taxid2))
            taxonCache->emplace(taxid2, TaxonBuilder::getInstance()->buildOne(taxid2));
        auto dist = *taxonCache->at(taxid1) - *taxonCache->at(taxid2);
        weight = dist[0] + dist[1];
        if (weightCache != nullptr)
            weightCache->emplace(key, weight);
    } else
        weight = weightCache->at(key);
    if (!weight) // TODO: there could be multiple pathways producing the same product in the same organism
        weight += (lhs->product()->id() != rhs->product()->id()) * 0.7;
    if (!weight) {
        int diff = lhs->module() - rhs->module();
        if (diff > 0)
            weight += 0.6;
        else if (std::abs(diff) > 1)
            weight += 0.5;
        else if (std::abs(diff) == 1)
            weight += (rhs->type() != DomainType::C) ? 0.4 : 0.1;
        else if (std::abs(diff) == 0 && rhs->type() == DomainType::C)
            weight += 0.2;
    }
    return weight;
}

#ifdef WITH_INTERNAL_XML
void nrps::startXml(xmlTextWriterPtr writer)
{
    xmlTextWriterSetIndent(writer, 1);
    xmlTextWriterSetIndentString(writer, BAD_CAST "    ");
    xmlTextWriterStartDocument(writer, nullptr, "UTF-8", nullptr);
}

void nrps::endXml(xmlTextWriterPtr writer)
{
    xmlTextWriterEndDocument(writer);
}
#endif
