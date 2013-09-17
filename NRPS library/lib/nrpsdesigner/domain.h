#ifndef NRPSDESIGNER_DOMAIN_H
#define NRPSDESIGNER_DOMAIN_H

#include "nrpsdesigner_export.h"
#include "global_enums.h"

#include <cstdint>
#include <string>
#include <memory>
#include <unordered_map>

#include <libxml/xmlwriter.h>

namespace nrps
{
class Origin;
class Product;
class NRPSDESIGNER_EXPORT Domain
{
public:
    virtual std::size_t hash() const;

    DomainType type() const;
    uint32_t id() const;
    uint32_t module() const;
    Origin* origin() const;
    Product* product() const;
    const std::string& description() const;
    const std::string& dnaSequencePfam() const;
    const std::string& dnaSequenceDefined() const;
    const std::string& geneName() const;
    const std::string& geneDescription() const;
    const std::string& nativePfamLinkerBefore() const;
    const std::string& nativePfamLinkerAfter() const;
    const std::string& nativeDefinedLinkerBefore() const;
    const std::string& nativeDefinedLinkerAfter() const;

    void toXml(xmlTextWriterPtr) const;

    void setId(uint32_t);
    void setModule(uint32_t);
    Origin* setOrigin(uint32_t);
    Product* setProduct(uint32_t);
    void setDescription(const std::string&);
    void setDescription(std::string&&);
    void setDnaSequencePfam(const std::string&);
    void setDnaSequencePfam(std::string&&);
    void setDnaSequenceDefined(const std::string&);
    void setDnaSequenceDefined(std::string&&);
    void setGeneDescription(const std::string&);
    void setGeneDescription(std::string&&);
    void setGeneName(const std::string&);
    void setGeneName(std::string&&);
    void setNativePfamLinkerBefore(const std::string&);
    void setNativePfamLinkerBefore(std::string&&);
    void setNativePfamLinkerAfter(const std::string&);
    void setNativePfamLinkerAfter(std::string&&);
    void setNativeDefinedLinkerBefore(const std::string&);
    void setNativeDefinedLinkerBefore(std::string&&);
    void setNativeDefinedLinkerAfter(const std::string&);
    void setNativeDefinedLinkerAfter(std::string&&);

protected:
    Domain(DomainType, uint32_t);
    virtual void writeXml(xmlTextWriterPtr) const;


private:
    void startXml(xmlTextWriterPtr) const;
    void endXml(xmlTextWriterPtr) const;

    DomainType m_type;
    uint32_t m_id;
    uint32_t m_module;
    Origin *m_origin;
    Product *m_product;
    std::string m_description;
    std::string m_dnaSeqPfam;
    std::string m_dnaSeqDefined;
    std::string m_geneDescription;
    std::string m_geneName;
    std::string m_nativePfamLinkerBefore;
    std::string m_nativePfamLinkerAfter;
    std::string m_nativeDefinedLinkerBefore;
    std::string m_nativeDefinedLinkerAfter;
};
}

namespace std
{
    template<>
    struct NRPSDESIGNER_EXPORT hash<nrps::Domain>
    {
    public:
        std::size_t operator()(const nrps::Domain &d) const;
    };
}

#endif
