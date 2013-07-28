#ifndef NRPSDESIGNER_DOMAIN_H
#define NRPSDESIGNER_DOMAIN_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"
#include "pathway.h"
#include "global_enums.h"

#include <cstdint>
#include <string>
#include <memory>
#include <unordered_map>

namespace nrps
{
class NRPSDESIGNER_EXPORT Domain : public AbstractDomainType
{
public:
    virtual std::size_t hash() const;
    virtual bool full() const;

    uint32_t domainId() const;
    uint32_t moduleId() const;
    const std::shared_ptr<Pathway>& pathway() const;
    const std::string& bioBrickId() const;
    const std::string& description() const;
    const std::string& dnaSequence() const;
    const std::string& nativeLinkerBefore() const;
    const std::string& nativeLinkerAfter() const;
    const std::string& refSeqId() const;
    const std::string& uniProtId() const;

    void setDomainId(uint32_t);
    void setModuleId(uint32_t);
    const std::shared_ptr<Pathway>& setPathway(uint32_t);
    void setBioBrickId(const std::string&);
    void setBioBrickId(std::string&&);
    void setDescription(const std::string&);
    void setDescription(std::string&&);
    void setDnaSequence(const std::string&);
    void setDnaSequence(std::string&&);
    void setNativeLinkerBefore(const std::string&);
    void setNativeLinkerBefore(std::string&&);
    void setNativeLinkerAfter(const std::string&);
    void setNativeLinkerAfter(std::string&&);
    void setRefSeqId(const std::string&);
    void setRefSeqId(std::string&&);
    void setUniProtId(const std::string&);
    void setUniProtId(std::string&&);


protected:
    Domain(DomainType type);

private:
    uint32_t m_domainId;
    uint32_t m_moduleId;
    uint32_t m_pathway;
    std::string m_bioBrickId;
    std::string m_description;
    std::string m_dnaSeq;
    std::string m_nativeLinkerBefore;
    std::string m_nativeLinkerAfter;
    std::string m_refSeqId;
    std::string m_uniProtId;

    static std::unordered_map<uint32_t, std::shared_ptr<Pathway>> s_pathways;
};

template<bool full>
class DomainBaseType
{
public:
    typedef typename std::conditional<full, Domain, AbstractDomainType>::type type;
};
}

#endif
