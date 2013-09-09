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

void *Taxon::s_handle = nullptr;
std::string Taxon::s_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&rettype=xml&id=";

Taxon::Taxon(uint32_t id) throw (NCBITaxonomyError, std::logic_error, std::system_error)
: m_id(id), m_parentId(0), m_full(false)
{
    fetch();
}

Taxon::Taxon(xmlNodePtr node) throw (std::logic_error)
: m_id(0), m_parentId(0), m_full(false)
{
    parseTaxon(node);
}

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

void Taxon::fetch() throw (NCBITaxonomyError, std::logic_error, std::system_error)
{
    if (s_handle == nullptr) {
        CURLcode init = curl_global_init(CURL_GLOBAL_DEFAULT);
        if (init != CURLE_OK)
            throw std::system_error(init, std::system_category(), "curl initialization failed!");
        else
            s_handle = curl_easy_init( );
    }
    std::string url = s_url + std::to_string(m_id);
    xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(nullptr, nullptr, nullptr, 0, url.c_str());

    CURLcode ret = curl_easy_setopt(s_handle, CURLOPT_URL, url.c_str());
    ret = curl_easy_setopt(s_handle, CURLOPT_NOSIGNAL, 1);
    ret = curl_easy_setopt(s_handle, CURLOPT_TIMEOUT, NetworkOptions::getInstance()->timeout());
    ret = curl_easy_setopt(s_handle, CURLOPT_WRITEFUNCTION, static_cast<size_t (*)( char*, size_t, size_t, void*)>([](char *ptr, size_t size, size_t nmemb, void *userdata){xmlParseChunk((xmlParserCtxtPtr)userdata, ptr, size * nmemb, 0);return size * nmemb;}));
    ret = curl_easy_setopt(s_handle, CURLOPT_WRITEDATA, ctxt);
    ret = curl_easy_perform(s_handle);
    long httpcode;
    curl_easy_getinfo(s_handle, CURLINFO_RESPONSE_CODE, &httpcode);
    if (ret > 0)
        throw NCBITaxonomyError(std::string("Something went wrong when fetching the NCBI taxonomy information. Error code: ") + std::to_string(ret) + std::string(";    ") + curl_easy_strerror(ret));
    if (httpcode != 200 && ret != CURLE_ABORTED_BY_CALLBACK)
        throw NCBITaxonomyError("The NCBI taxonomy database seems to be experiencing problems.");
    xmlParseChunk(ctxt, nullptr, 0, 1);
    xmlDocPtr doc = ctxt->myDoc;
    xmlFreeParserCtxt(ctxt);
    xmlNodePtr node = xmlDocGetRootElement(doc);
    if (node == nullptr || XMLCMP(node, TAXASET_NODE)) {
        xmlFreeDoc(doc);
        throw NCBITaxonomyError("Invalid taxonomy XML");
    }
    if (xmlChildElementCount(node) > 1) {
        xmlFreeDoc(doc);
        throw NCBITaxonomyError("More than 1 taxon in taxaset");
    }
    node = xmlFirstElementChild(node);
    try {
        parseTaxon(node);
    } catch (std::exception &e) {
        xmlFreeDoc(doc);
        throw;
    }
    xmlFreeDoc(doc);
    m_full = true;
}

void Taxon::parseTaxon(xmlNodePtr node) throw (std::logic_error)
{
    if (XMLCMP(node, TAXON_NODE))
        throw NCBITaxonomyError("Not a taxon node.");
    node = xmlFirstElementChild(node);
    while (node != nullptr) {
        if (!XMLCMP(node, TAXID_NODE)) {
            uint32_t id = std::atoi(XMLTXT(node));
            if (m_id > 0 && id != m_id)
                throw std::logic_error("TaxID returned by server differs from requested.");
            else
                m_id = id;
        } else if (m_scientificName.empty() && !XMLCMP(node, SCIENTIFICNAME_NODE))
            m_scientificName = XMLTXT(node);
        else if (m_names.empty() && !XMLCMP(node, OTHERNAMES_NODE))
            parseOtherNames(node);
        else if (!m_parentId && !XMLCMP(node, PARENTTAXID_NODE))
            m_parentId = std::atoi(XMLTXT(node));
        else if (!XMLCMP(node, RANK_NODE))
            m_rank = parseRank(node->children->content);
        else if (!XMLCMP(node, DIVISION_NODE))
            m_division = XMLTXT(node);
        else if (m_gc.name.empty() && !XMLCMP(node, GENETICCODE_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            GeneticCode c;
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, GCID_NODE))
                    c.id = std::atoi(XMLTXT(nnode));
                else if (!XMLCMP(nnode, GCNAME_NODE))
                    c.name = XMLTXT(nnode);
                nnode = xmlNextElementSibling(nnode);
            }
            m_gc = c;
        } else if (m_mitoGc.name.empty() && !XMLCMP(node, MITOGENETICCODE_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            GeneticCode c;
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, MGCID_NODE))
                    c.id = std::atoi(XMLTXT(nnode));
                else if (!XMLCMP(nnode, MGCNAME_NODE))
                    c.name = XMLTXT(nnode);
                nnode = xmlNextElementSibling(nnode);
            }
            m_mitoGc = std::move(c);
        } else if (m_lineage.empty() && !XMLCMP(node, LINEAGEEX_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            while (nnode != nullptr) {
                m_lineage.emplace_back(new Taxon(nnode));
                nnode = xmlNextElementSibling(nnode);
            }
        } else if (!XMLCMP(node, CREATEDATE_NODE))
            m_createDate = parseDate(node);
        else if (!XMLCMP(node, UPDATEDATE_NODE))
            m_updateDate = parseDate(node);
        else if (!XMLCMP(node, PUBDATE_NODE))
            m_pubDate = parseDate(node);
        node = xmlNextElementSibling(node);
    }
}

void Taxon::parseOtherNames(xmlNodePtr node)
{
    node = xmlFirstElementChild(node);
    while (node != nullptr) {
        if (!XMLCMP(node, SYNONYM_NODE))
            m_synonyms.push_back(XMLTXT(node));
        else if (!XMLCMP(node, INCLUDES_NODE))
            m_includes.push_back(XMLTXT(node));
        else if (XMLCMP(node, NAME_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            Name n;
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, CLASSCDE_NODE))
                    n.classCDE = XMLTXT(nnode);
                else if (!XMLCMP(nnode, DISPNAME_NODE))
                    n.dispName = XMLTXT(nnode);
                nnode = xmlNextElementSibling(nnode);
            }
            m_names.push_back(std::move(n));
        }
        node = xmlNextElementSibling(node);
    }
}

Taxon::Rank Taxon::parseRank(xmlChar *r)
{
    if (!xmlStrcmp(r, BAD_CAST RANK_NORANK))
        return Rank::NoRank;
    else if (!xmlStrcmp(r, BAD_CAST RANK_SUPERKINGDOM))
        return Rank::Superkingdom;
    else if (!xmlStrcmp(r, BAD_CAST RANK_PHYLUM))
        return Rank::Phylum;
    else if (!xmlStrcmp(r, BAD_CAST RANK_CLASS))
        return Rank::Class;
    else if (!xmlStrcmp(r, BAD_CAST RANK_ORDER))
        return Rank::Order;
    else if (!xmlStrcmp(r, BAD_CAST RANK_FAMILY))
        return Rank::Family;
    else if (!xmlStrcmp(r, BAD_CAST RANK_GENUS))
        return Rank::Genus;
    else if (!xmlStrcmp(r, BAD_CAST RANK_SPECIES))
        return Rank::Species;
}

std::chrono::system_clock::time_point Taxon::parseDate(xmlNodePtr node)
{
    std::tm t;
#ifdef USE_STRPTIME
    auto test = strptime(XMLTXT(node), TIME_FORMAT, &t);
#else
    std::istringstream s(XMLTXT(node));
    s.imbue(std::locale("en_US"));
    s >> std::get_time(&t, TIME_FORMAT);
#endif
    auto test2 = std::mktime(&t);
    return std::chrono::system_clock::from_time_t(std::mktime(&t));
}

std::array<uint8_t, 2> Taxon::operator-(const Taxon &o) const
{
    return diff(o, Rank::Species);
}

std::array<uint8_t, 2> Taxon::diff(const Taxon &o, Taxon::Rank r) const
{
    std::array<uint8_t, 2> ret;
    int lca;
    bool species = false;
    for (lca = 0; lca < m_lineage.size() && lca < o.m_lineage.size() && m_lineage[lca]->id() == o.m_lineage[lca]->id(); ++lca) {
        if (m_lineage[lca]->rank() == r)
            species = true;
    }
    uint8_t tsp = 0, tssp = 0, osp = 0, ossp = 0;
    for (int i = lca; i < m_lineage.size(); ++i) {
        !species ? ++tsp : ++tssp;
        if (m_lineage[i]->rank() == r)
            species = true;
    }
    for (int i = lca; i < o.m_lineage.size(); ++i) {
        !species ? ++osp : ++ossp;
        if (o.m_lineage[i]->rank() == r)
            species = true;
    }
    ret[0] = std::max(tsp, osp);
    ret[1] = std::max(tssp, ossp);
    return ret;
}
