#ifndef NRPSDESIGNER_TAXON_H
#define NRPSDESIGNER_TAXON_H

#include "nrpsdesigner_export.h"
//#include "taxonbuilder.h"
#include "exceptions.h"

#include <vector>
#include <string>
#include <chrono>
#include <memory>
#include <array>
#ifdef _CXX_CLANG // Clang has std::array defined in tuple header
#include <tuple>
#endif
#include <system_error>

#include <libxml/tree.h>

namespace nrps
{
class NRPSDESIGNER_EXPORT Taxon
{
friend class TaxonBuilder;
public:
    struct Name
    {
    public:
        std::string classCDE;
        std::string dispName;
    };
    struct GeneticCode
    {
    public:
        uint32_t id;
        std::string name;
    };
    enum class Rank {NoRank, Superkingdom, Phylum, Class, Order, Family, Genus, Species};

    bool full() const;
    uint32_t id() const;
    uint32_t parentId() const;
    const std::string& scientificName() const;
    const std::vector<std::string>& synonyms() const;
    const std::vector<std::string>& includes() const;
    const std::vector<Name>& names() const;
    Rank rank() const;
    const std::string& division() const;
    const GeneticCode& geneticCode() const;
    const GeneticCode& mtGeneticCode() const;
    const std::vector<std::shared_ptr<Taxon>>& lineage() const;
    std::chrono::system_clock::time_point createDate() const;
    std::chrono::system_clock::time_point updateDate() const;
    std::chrono::system_clock::time_point pubDate() const;

    void fetch() throw (NCBITaxonomyError, std::logic_error, std::system_error);

    std::array<uint8_t, 2> operator-(const Taxon&) const;
    std::array<uint8_t, 2> diff(const Taxon&, Rank r = Rank::Species) const;

private:
    Taxon();
    Taxon(uint32_t);

    uint32_t m_id;
    uint32_t m_parentId;
    bool m_full;
    std::string m_scientificName;
    std::vector<std::string> m_synonyms;
    std::vector<std::string> m_includes;
    std::vector<Name> m_names;
    Rank m_rank;
    std::string m_division;
    GeneticCode m_gc;
    GeneticCode m_mitoGc;
    std::vector<std::shared_ptr<Taxon>> m_lineage;
    std::chrono::system_clock::time_point m_createDate;
    std::chrono::system_clock::time_point m_updateDate;
    std::chrono::system_clock::time_point m_pubDate;
};
}

#endif
