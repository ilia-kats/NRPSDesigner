#include "taxonbuilder.h"
#include "global_enums.h"
#include "networkoptions.h"
#include "globals_internal.h"

#include <iostream>
#include <locale>
#include <cctype>
#include <stack>
#include <unordered_map>
#include <algorithm>

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
#define RANK_SPECIESGROUP "species group"
#define RANK_SPECIES "species"

#define TIME_FORMAT "%Y/%m/%d %H:%M:%S"

using namespace nrps;
namespace po = boost::program_options;

std::string TaxonBuilder::s_url("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&rettype=xml&id=");
TaxonBuilder *TaxonBuilder::s_instance = nullptr;

static std::ctype_base::mask whitespacePipeTable[std::ctype<char>::table_size];

TaxonBuilder::TaxonBuilder()
: m_handle(nullptr)
{}

TaxonBuilder::~TaxonBuilder()
{
    s_instance = nullptr;
    curl_global_cleanup();
    m_fallbackFile.close();
}

void TaxonBuilder::init()
{
    if (m_handle == nullptr) {
        CURLcode ret = curl_global_init(CURL_GLOBAL_DEFAULT);
        if (ret != CURLE_OK)
            throw std::system_error(ret, std::system_category(), "curl initialization failed!");
        else
            m_handle = curl_easy_init( );
        ret = curl_easy_setopt(m_handle, CURLOPT_NOSIGNAL, 1);
        ret = curl_easy_setopt(m_handle, CURLOPT_TIMEOUT, NetworkOptions::getInstance()->timeout());
        ret = curl_easy_setopt(m_handle, CURLOPT_WRITEFUNCTION, static_cast<size_t (*)( char*, size_t, size_t, void*)>([](char *ptr, size_t size, size_t nmemb, void *userdata){xmlParseChunk((xmlParserCtxtPtr)userdata, ptr, size * nmemb, 0);return size * nmemb;}));
    }
    if (!m_fallbackFile.is_open() && !m_fallbackFilePath.empty()) {
        m_fallbackFile.open(m_fallbackFilePath);
        if (m_fallbackFile.fail())
            std::clog << "WARNING: taxonomy dump file could not be opened" << std::endl;
        else {
            std::copy_n(std::ctype<char>::classic_table(),  std::ctype<char>::table_size, whitespacePipeTable);
            for (int i = 0; i <  std::ctype<char>::table_size; ++i) {
                if (whitespacePipeTable[i] & std::ctype_base::space) {
                    whitespacePipeTable[i] &= ~std::ctype_base::space;
                }
            }
            whitespacePipeTable['|'] = std::ctype_base::space;
            whitespacePipeTable['\t'] = std::ctype_base::space;
            std::locale l(m_fallbackFile.getloc(), new std::ctype<char>(whitespacePipeTable));
            m_fallbackFile.imbue(l);
        }
    }
}

TaxonBuilder* TaxonBuilder::getInstance()
{
    if (s_instance == nullptr) {
        s_instance = new TaxonBuilder();
    }
    return s_instance;
}

po::options_description TaxonBuilder::options()
{
    po::options_description options("NCBI taxonomy options");
    options.add_options()("dumpfile", po::value<std::string>(&m_fallbackFilePath), "Path to the NCBI taxonomy nodes.dmp file to be used as fallback when no internet connection is available or the NCBI server is experiencing technical difficulties. Set the NRPSDESIGNER_DUMPFILE environment variable to avoid specifying this on the command-line.");
    return options;
}

std::string TaxonBuilder::mapEnvironment(std::string env) {
    if (env == "NRPSDESIGNER_DUMPFILE")
        return "dumpfile";
    else
        return std::string();
}

std::shared_ptr<Taxon> TaxonBuilder::buildOne(uint32_t id) throw (NCBITaxonomyError, std::logic_error, std::system_error)
{
    init();
    std::shared_ptr<Taxon> ret(new Taxon(id));
    process(ret);
    return ret;
}

std::shared_ptr<Taxon> TaxonBuilder::buildMany(uint32_t id)
{
    init();
    std::shared_ptr<Taxon> ret(new Taxon(id));
    try {
        fetch(ret);
    } catch (const std::exception &e) {
        std::clog << "WARNING: Could not fetch taxon " + std::to_string(id) + " from NCBI taxonomy, falling back to local copy. Reason: " << e.what() << std::endl;
        m_toProcess.emplace(id, ret);
    }
    return ret;
}

void TaxonBuilder::process() throw (TaxonomyDumpError, std::logic_error)
{
    init();
    if (!m_fallbackFile.is_open() && !m_toProcess.empty())
        throw TaxonomyDumpError("Fallback from taxonomy dump file requested, but it is not open");
    std::unordered_map<uint32_t, std::shared_ptr<Taxon>> taxa;
    while (!m_toProcess.empty()) {
        for (auto iter = m_toProcess.begin(); iter != m_toProcess.end(); ++iter) {
            seek(iter->first);
            parseDump(iter->second);
            uint32_t parent = iter->second->parentId();
            if (parent != 1)
                m_toProcess.emplace(parent, std::shared_ptr<Taxon>(new Taxon(parent)));
            taxa.emplace(iter->first, iter->second);
            m_toProcess.erase(iter);
        }
    }
    for (const auto &taxon : taxa) {
        uint32_t id = taxon.second->id();
        uint32_t pid = taxon.second->parentId();
        if (taxon.second->lineage().empty() && pid != 1) {
            std::stack<std::shared_ptr<Taxon>> stck;
            std::shared_ptr<Taxon> parent = taxon.second;
            stck.push(parent);
            while (parent->lineage().empty() && pid != 1) {
                parent = taxa[pid];
                stck.push(parent);
                pid = parent->parentId();
                id = parent->id();
            }
            std::vector<std::shared_ptr<Taxon>> lineage = parent->lineage();
            std::shared_ptr<Taxon> t = stck.top();
            stck.pop();
            lineage.push_back(t);
            while (!stck.empty()) {
                t = stck.top();
                stck.pop();
                t->m_lineage = lineage;
                lineage.push_back(t);
            }
        }
    }
}

void TaxonBuilder::process(const std::shared_ptr<Taxon> &t) throw (NCBITaxonomyError, std::logic_error, std::system_error)
{
    if (!t->full()) {
        init();
        try {
            fetch(t);
        } catch (const std::exception&) {
            if (m_fallbackFile.is_open()) {
                std::vector<std::shared_ptr<Taxon>> lineage;
                std::shared_ptr<Taxon> taxon = t;
                while (taxon->id() != 1) {
                    seek(taxon->id());
                    parseDump(taxon);
                    lineage.push_back(taxon);
                    taxon = std::shared_ptr<Taxon>(new Taxon(taxon->parentId()));
                }
                for (auto i = lineage.end() - 2; i >= lineage.begin(); --i) {
                    (*i)->m_lineage.resize(std::distance(i, lineage.end()));
                    std::reverse_copy(i, lineage.end(), (*i)->m_lineage.begin());
                }
            }
            else
                throw;
        }
    }
}

void TaxonBuilder::seek(uint32_t id) throw (TaxonomyDumpError)
{
    auto iter = m_filePos.find(id);
    if (iter != m_filePos.end()) {
        m_fallbackFile.seekg(iter->second);
    } else {
        decltype(m_fallbackFile)::int_type ch;
        if (!m_filePos.empty() && m_fallbackFile.tellg() < m_filePos.rbegin()->second) {
            m_fallbackFile.seekg(m_filePos.rbegin()->second);
            while ((ch = m_fallbackFile.get()) != '\n' && ch != decltype(m_fallbackFile)::traits_type::eof());
        }
        else
            m_fallbackFile.seekg(0);
        uint32_t fid = 0;
        size_t pos;
        pos = m_fallbackFile.tellg();
        m_fallbackFile >> fid;
        m_filePos.emplace_hint(m_filePos.end(), fid, pos);
        while (fid < id && ch != decltype(m_fallbackFile)::traits_type::eof()) {
            while ((ch = m_fallbackFile.get()) != '\n' && ch != decltype(m_fallbackFile)::traits_type::eof());
            pos = m_fallbackFile.tellg();
            m_fallbackFile >> fid;
            m_filePos.emplace_hint(m_filePos.end(), fid, pos);
        }
        if (fid < id) // EOF
            throw TaxonomyDumpError("Requested ID could not be found in file. Requested: " + std::to_string(id) + "; last in file: " + std::to_string(fid));
        if (m_fallbackFile.bad())
            throw TaxonomyDumpError("Error while seeking in taxonomy dump file");
        m_fallbackFile.seekg(pos);
    }
}

void TaxonBuilder::parseDump(const std::shared_ptr<Taxon> &t) throw (TaxonomyDumpError, std::logic_error)
{
    uint32_t id;
    m_fallbackFile >> id;
    if (id != t->m_id)
        throw std::logic_error("Found TaxID differs from requested. Requested: " + std::to_string(t->m_id) + "; found: " + std::to_string(id));
    if (m_fallbackFile.fail() || m_fallbackFile.bad())
        throw TaxonomyDumpError("TaxID could not be parsed from dump");
    std::locale l = m_fallbackFile.getloc();
    m_fallbackFile >> t->m_parentId;
    if (m_fallbackFile.fail() || m_fallbackFile.bad())
        throw TaxonomyDumpError("Taxon parent ID could not be parsed from dump");
    std::string buf;
    m_fallbackFile >> buf;
    if (m_fallbackFile.fail() || m_fallbackFile.bad())
        throw TaxonomyDumpError("Taxon rank could not be parsed from dump");
    t->m_rank = parseRank(buf.c_str());
    if (m_fallbackFile.peek() != '\t')
        m_fallbackFile >> buf;
    m_fallbackFile >> buf >> buf;
    m_fallbackFile >> t->m_gc.id;
    if (m_fallbackFile.fail() || m_fallbackFile.bad())
        throw TaxonomyDumpError("Taxon genetic code could not be parsed from dump");
    m_fallbackFile >> buf;
    m_fallbackFile >> t->m_mitoGc.id;
    if (m_fallbackFile.fail() || m_fallbackFile.bad())
        throw TaxonomyDumpError("Taxon mt genetic code could not be parsed from dump");
    t->m_full = false;
}

void TaxonBuilder::fetch(const std::shared_ptr<Taxon> &t) throw (NCBITaxonomyError, std::logic_error, std::system_error)
{
    std::string url = s_url + std::to_string(t->id());
    xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(nullptr, nullptr, nullptr, 0, url.c_str());
    CURLcode ret = curl_easy_setopt(m_handle, CURLOPT_URL, url.c_str());
    ret = curl_easy_setopt(m_handle, CURLOPT_WRITEDATA, ctxt);
    ret = curl_easy_perform(m_handle);
    long httpcode;
    curl_easy_getinfo(m_handle, CURLINFO_RESPONSE_CODE, &httpcode);
    if (ret > 0)
        throw NCBITaxonomyError(std::string("Something went wrong when fetching the NCBI taxonomy information. Error code: ") + std::to_string(ret) + std::string("; ") + curl_easy_strerror(ret));
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
        parseTaxon(node, t);
    } catch (std::exception &e) {
        xmlFreeDoc(doc);
        throw;
    }
    xmlFreeDoc(doc);
    t->m_full = true;
}

void TaxonBuilder::parseTaxon(xmlNodePtr node, const std::shared_ptr<Taxon> &t) throw (std::logic_error)
{
    if (XMLCMP(node, TAXON_NODE))
        throw NCBITaxonomyError("Not a taxon node.");
    node = xmlFirstElementChild(node);
    while (node != nullptr) {
        if (!XMLCMP(node, TAXID_NODE)) {
            uint32_t id = std::atoi(XMLTXT(node));
            if (t->m_id > 0 && id != t->m_id)
                throw std::logic_error("TaxID returned by server differs from requested. Requested: " + std::to_string(t->m_id) + "; found: " + std::to_string(id));
            else
                t->m_id = id;
        } else if (t->m_scientificName.empty() && !XMLCMP(node, SCIENTIFICNAME_NODE))
            t->m_scientificName = XMLTXT(node);
        else if (t->m_names.empty() && !XMLCMP(node, OTHERNAMES_NODE))
            parseOtherNames(node, t);
        else if (!t->m_parentId && !XMLCMP(node, PARENTTAXID_NODE))
            t->m_parentId = std::atoi(XMLTXT(node));
        else if (!XMLCMP(node, RANK_NODE))
            t->m_rank = parseRank((const char*)node->children->content);
        else if (!XMLCMP(node, DIVISION_NODE))
            t->m_division = XMLTXT(node);
        else if (t->m_gc.name.empty() && !XMLCMP(node, GENETICCODE_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            Taxon::GeneticCode c;
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, GCID_NODE))
                    c.id = std::atoi(XMLTXT(nnode));
                else if (!XMLCMP(nnode, GCNAME_NODE))
                    c.name = XMLTXT(nnode);
                nnode = xmlNextElementSibling(nnode);
            }
            t->m_gc = c;
        } else if (t->m_mitoGc.name.empty() && !XMLCMP(node, MITOGENETICCODE_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            Taxon::GeneticCode c;
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, MGCID_NODE))
                    c.id = std::atoi(XMLTXT(nnode));
                else if (!XMLCMP(nnode, MGCNAME_NODE))
                    c.name = XMLTXT(nnode);
                nnode = xmlNextElementSibling(nnode);
            }
            t->m_mitoGc = std::move(c);
        } else if (t->m_lineage.empty() && !XMLCMP(node, LINEAGEEX_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            while (nnode != nullptr) {
                std::shared_ptr<Taxon> lt(new Taxon());
                parseTaxon(nnode, lt);
                t->m_lineage.push_back(lt);
                nnode = xmlNextElementSibling(nnode);
            }
        } else if (!XMLCMP(node, CREATEDATE_NODE))
            t->m_createDate = parseDate(node);
        else if (!XMLCMP(node, UPDATEDATE_NODE))
            t->m_updateDate = parseDate(node);
        else if (!XMLCMP(node, PUBDATE_NODE))
            t->m_pubDate = parseDate(node);
        node = xmlNextElementSibling(node);
    }
}

void TaxonBuilder::parseOtherNames(xmlNodePtr node, const std::shared_ptr<Taxon> &t)
{
    node = xmlFirstElementChild(node);
    while (node != nullptr) {
        if (!XMLCMP(node, SYNONYM_NODE))
            t->m_synonyms.push_back(XMLTXT(node));
        else if (!XMLCMP(node, INCLUDES_NODE))
            t->m_includes.push_back(XMLTXT(node));
        else if (XMLCMP(node, NAME_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            Taxon::Name n;
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, CLASSCDE_NODE))
                    n.classCDE = XMLTXT(nnode);
                else if (!XMLCMP(nnode, DISPNAME_NODE))
                    n.dispName = XMLTXT(nnode);
                nnode = xmlNextElementSibling(nnode);
            }
            t->m_names.push_back(std::move(n));
        }
        node = xmlNextElementSibling(node);
    }
}

Taxon::Rank TaxonBuilder::parseRank(const char *r)
{
    if (!strcmp(r, RANK_NORANK))
        return Taxon::Rank::NoRank;
    else if (!strcmp(r, RANK_SUPERKINGDOM))
        return Taxon::Rank::Superkingdom;
    else if (!strcmp(r, RANK_PHYLUM))
        return Taxon::Rank::Phylum;
    else if (!strcmp(r, RANK_CLASS))
        return Taxon::Rank::Class;
    else if (!strcmp(r, RANK_ORDER))
        return Taxon::Rank::Order;
    else if (!strcmp(r, RANK_FAMILY))
        return Taxon::Rank::Family;
    else if (!strcmp(r, RANK_GENUS))
        return Taxon::Rank::Genus;
    else if (!strcmp(r, RANK_SPECIESGROUP))
        return Taxon::Rank::SpeciesGroup;
    else if (!strcmp(r, RANK_SPECIES))
        return Taxon::Rank::Species;
}

std::chrono::system_clock::time_point TaxonBuilder::parseDate(xmlNodePtr node)
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
