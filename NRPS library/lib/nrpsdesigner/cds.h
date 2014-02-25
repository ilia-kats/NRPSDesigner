#ifndef NRPSDESIGNER_CDS_H
#define NRPSDESIGNER_CDS_H

#include "config.h"
#include "nrpsdesigner_export.h"

#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>

#ifdef WITH_INTERNAL_XML
#include <libxml/xmlwriter.h>
#endif

namespace nrps
{
class Origin;
class Product;
class NRPSDESIGNER_EXPORT Cds
{
public:
    ~Cds();
    uint32_t id() const;
    const std::string& geneName() const;
    const std::string& dnaSequence() const;
    const std::string& description() const;
    Origin* origin() const;
    Product* product() const;

    void setId(uint32_t);
    void setGeneName(const std::string&);
    void setGeneName(std::string&&);
    void setDnaSequence(const std::string&);
    void setDnaSequence(std::string&&);
    void setDescription(const std::string&);
    void setDescription(std::string&&);
    Origin* setOrigin(uint32_t);
    Product* setProduct(uint32_t);

#ifdef WITH_INTERNAL_XML
    void toXml(xmlTextWriterPtr) const;
#endif
    static Cds* makeCds(uint32_t);

private:
    Cds(uint32_t);

    static std::unordered_map<uint32_t, Cds*> s_cdss;

    uint32_t m_id;
    Origin *m_origin;
    Product *m_product;
    std::string m_geneName;
    std::string m_dnaSequence;
    std::string m_description;
};
}

#endif

