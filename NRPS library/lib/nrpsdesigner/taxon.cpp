#include "taxon.h"
#include "global_enums.h"
#include "networkoptions.h"

#include <system_error>
#include <stdexcept>
#include <sstream>
#include <ctime>

#include <libxml/parser.h>

#include <curl/curl.h>

#define TAXASET_NODE "TaxaSet"
#define TAXON_NODE "Taxon"
#define TAXID_NODE "TaxId"
#define SCIENTIFICNAME_NODE "ScientificName"
#define OTHERNAMES_NODE "OtherNames"
#define SYNONYM_NODE "Synonym"
#define INCLUDES_NODE "Includes"
#define NAME_NODE "Name"
#define CLASSCDE_NODE "ClassCDE"
#define DISPNAME_NODE "DispName"
#define PARENTTAXID_NODE "ParentTaxId"
#define RANK_NODE "Rank"
#define DIVISION_NODE "Division"
#define GENETICCODE_NODE "GeneticCode"
#define GCID_NODE "GCId"
#define GCNAME_NODE "GCName"
#define MITOGENETICCODE_NODE "MitoGeneticCode"
#define MGCID_NODE "MGCId"
#define MGCNAME_NODE "MGCName"
#define LINEAGE_NODE "Lineage"
#define LINEAGEEX_NODE "LineageEx"
#define CREATEDATE_NODE "CreateDate"
#define UPDATEDATE_NODE "UpdateDate"
#define PUBDATE_NODE "PubDate"

#define RANK_NORANK "no rank"
#define RANK_SUPERKINGDOM "superkingdom"
#define RANK_PHYLUM "phylum"
#define RANK_CLASS "class"
#define RANK_ORDER "order"
#define RANK_FAMILY "family"
#define RANK_GENUS "genus"
#define RANK_SPECIES "species"

#define TIME_FORMAT "%Y/%m/%d %H:%M:%S"

#define XMLCMP(n, y) xmlStrcmp(n->name, BAD_CAST y)
#define XMLTXT(n) (const char*)n->children->content

using namespace nrps;

Taxon::Taxon(uint32_t id)
: m_id(id), m_parentId(0), m_full(false)
{}

Taxon::Taxon()
: m_id(0), m_parentId(0), m_full(false)
{}

bool Taxon::full() const
{
    return m_full;
}

uint32_t Taxon::id() const
{
    return m_id;
}

uint32_t Taxon::parentId() const
{
    return m_parentId;
}

const std::string& Taxon::scientificName() const
{
    return m_scientificName;
}

const std::vector<std::string>& Taxon::synonyms() const
{
    return m_synonyms;
}

const std::vector<std::string>& Taxon::includes() const
{
    return m_includes;
}

const std::vector<Taxon::Name>& Taxon::names() const
{
    return m_names;
}

Taxon::Rank Taxon::rank() const
{
    return m_rank;
}

const std::string& Taxon::division() const
{
    return m_division;
}

const Taxon::GeneticCode& Taxon::geneticCode() const
{
    return m_gc;
}

const Taxon::GeneticCode& Taxon::mtGeneticCode() const
{
    return m_mitoGc;
}

const std::vector<std::shared_ptr<Taxon>>& Taxon::lineage() const
{
    return m_lineage;
}

std::chrono::system_clock::time_point Taxon::createDate() const
{
    return m_createDate;
}

std::chrono::system_clock::time_point Taxon::updateDate() const
{
    return m_updateDate;
}

std::chrono::system_clock::time_point Taxon::pubDate() const
{
    return m_pubDate;
}

std::array<uint8_t, 2> Taxon::operator-(const Taxon &o) const
{
    return diff(o, Rank::Species);
}

std::array<uint8_t, 2> Taxon::diff(const Taxon &o, Taxon::Rank r) const
{
    if (id() == o.id())
        return std::array<uint8_t, 2>({0, 0});
    std::array<uint8_t, 2> ret;
    int lca;
    bool tspecies = false, ospecies = false;
    for (lca = 0; lca < m_lineage.size() && lca < o.m_lineage.size() && m_lineage[lca]->id() == o.m_lineage[lca]->id(); ++lca) {
        if (m_lineage[lca]->rank() == r)
            tspecies = ospecies = true;
    }
    uint8_t tsp = 0, tssp = 0, osp = 0, ossp = 0;
    for (int i = lca; i < m_lineage.size(); ++i) {
        !tspecies ? ++tsp : ++tssp;
        if (m_lineage[i]->rank() == r)
            tspecies = true;
    }
    !tspecies ? ++tsp : ++tssp;
    for (int i = lca; i < o.m_lineage.size(); ++i) {
        !ospecies ? ++osp : ++ossp;
        if (o.m_lineage[i]->rank() == r)
            ospecies = true;
    }
    !ospecies ? ++osp : ++ossp;
    ret[0] = std::max(tsp, osp);
    ret[1] = std::max(tssp, ossp);
    return ret;
}
