#ifndef NRPSDESIGNER_DOMAIN_H
#define NRPSDESIGNER_DOMAIN_H

#include "config.h"
#include "nrpsdesigner_export.h"
#include "global_enums.h"

#include <cstdint>
#include <string>
#include <memory>
#include <unordered_map>

#ifdef WITH_INTERNAL_XML
#include <libxml/xmlwriter.h>
#endif

namespace nrps
{
class Origin;
class Product;
class Cds;
class NRPSDESIGNER_EXPORT Domain
{
public:
    typedef struct {
        uint32_t stop;
        uint32_t next_start;
    } workingBorders;
    virtual std::size_t hash() const;

    DomainType type() const;
    uint32_t id() const;
    uint32_t module() const;
    Origin* origin() const;
    Product* product() const;
    Cds* cds() const;
    const std::string& description() const;
    std::string dnaSequencePfam() const;
    std::string dnaSequenceDefined() const;
    std::string dnaSequence(Domain *prev = nullptr, Domain *next = nullptr) const;
    std::string pfamLinkerBefore() const;
    std::string pfamLinkerAfter() const;
    std::string definedLinkerBefore() const;
    std::string definedLinkerAfter() const;
    uint32_t pfamStart() const;
    uint32_t pfamStop() const;
    uint32_t definedStart() const;
    uint32_t definedStop() const;
    uint32_t pfamLinkerStart() const;
    uint32_t pfamLinkerStop() const;
    uint32_t definedLinkerStart() const;
    uint32_t definedLinkerStop() const;
    const std::unordered_map<uint32_t, workingBorders>& workingNextDomains() const;
    bool worksWithNextDomain(uint32_t) const;
    const workingBorders& workingDomainBorders(uint32_t) const;

#ifdef WITH_INTERNAL_XML
    void toXml(xmlTextWriterPtr) const;
#endif
    std::string toString() const;

    void setId(uint32_t);
    void setModule(uint32_t);
    Origin* setOrigin(uint32_t);
    Product* setProduct(uint32_t);
    Cds* setCds(uint32_t);
    void setDescription(const std::string&);
    void setDescription(std::string&&);
    void setPfamStart(uint32_t);
    void setPfamStop(uint32_t);
    void setDefinedStart(uint32_t);
    void setDefinedStop(uint32_t);
    void setPfamLinkerStart(uint32_t);
    void setPfamLinkerStop(uint32_t);
    void setDefinedLinkerStart(uint32_t);
    void setDefinedLinkerStop(uint32_t);
    void addWorkingNextDomain(uint32_t, workingBorders);
    void removeWorkingNextDomain(uint32_t);

protected:
    Domain(DomainType, uint32_t);
#ifdef WITH_INTERNAL_XML
    virtual void writeXml(xmlTextWriterPtr) const;
#endif

private:
#ifdef WITH_INTERNAL_XML
    void startXml(xmlTextWriterPtr) const;
    void endXml(xmlTextWriterPtr) const;
#endif

    DomainType m_type;
    uint32_t m_id;
    uint32_t m_module;
    Cds *m_cds;
    std::string m_description;
    uint32_t m_pfamStart;
    uint32_t m_pfamStop;
    uint32_t m_definedStart;
    uint32_t m_definedStop;
    uint32_t m_pfamLinkerStart;
    uint32_t m_pfamLinkerStop;
    uint32_t m_definedLinkerStart;
    uint32_t m_definedLinkerStop;
    std::unordered_map<uint32_t, workingBorders> m_workingNextDomains;
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
