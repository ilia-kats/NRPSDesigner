#ifndef NRPSDESIGNER_TAXONBUILDER_H
#define NRPSDESIGNER_TAXONBUILDER_H

#include "nrpsdesigner_export.h"
#include "taxon.h"

#include "exceptions.h"

#include <vector>
#include <string>
#include <chrono>
#include <memory>
#include <map>
#include <array>
#include <fstream>
#ifdef _CXX_CLANG // Clang has std::array defined in tuple header
#include <tuple>
#endif

#include <libxml/tree.h>

#include <boost/program_options/options_description.hpp>

namespace nrps
{
class NRPSDESIGNER_EXPORT TaxonBuilder
{
public:
    ~TaxonBuilder();
    std::shared_ptr<Taxon> buildOne(uint32_t) throw (NCBITaxonomyError, std::logic_error, std::system_error);
    std::shared_ptr<Taxon> buildMany(uint32_t);
    bool toProcess() const;
    void process() throw (TaxonomyDumpError, std::logic_error);
    void process(const std::shared_ptr<Taxon>&) throw (NCBITaxonomyError, std::logic_error, std::system_error);

    boost::program_options::options_description options();
    std::string mapEnvironment(std::string);

    static TaxonBuilder* getInstance();

private:
    TaxonBuilder();
    void init();
    void fetch(const std::shared_ptr<Taxon>&) throw (NCBITaxonomyError, std::logic_error, std::system_error);
    void parseTaxon(xmlNodePtr, const std::shared_ptr<Taxon>&) throw (std::logic_error);
    void parseOtherNames(xmlNodePtr, const std::shared_ptr<Taxon>&);
    Taxon::Rank parseRank(const char*);
    std::chrono::system_clock::time_point parseDate(xmlNodePtr);
    void seek(uint32_t) throw (TaxonomyDumpError);
    void parseDump(const std::shared_ptr<Taxon>&) throw (TaxonomyDumpError, std::logic_error);

    void *m_handle;
    std::map<uint32_t, std::shared_ptr<Taxon>> m_toProcess;
    std::string m_fallbackFilePath;
    std::ifstream m_fallbackFile;
    std::map<uint32_t, size_t> m_filePos;

    static std::string s_url;
    static TaxonBuilder* s_instance;
};
}

#endif
