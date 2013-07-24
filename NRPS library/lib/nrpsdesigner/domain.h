#ifndef NRPSDESIGNER_DOMAIN_H
#define NRPSDESIGNER_DOMAIN_H

#include "nrpsdesigner_export.h"
#include "abstractdomaintype.h"
#include "pathway.h"
#include "global_enums.h"

#include <cstdint>
#include <string>
#include <memory>

namespace nrps
{
class NRPSDESIGNER_EXPORT Domain : public AbstractDomainType
{
public:
    virtual std::size_t hash() const;
    virtual bool full() const;

    uint32_t domainId() const;
    uint32_t moduleId() const;
    const Pathway& pathway() const;
    const std::string& bioBrickId() const;
    const std::string& description() const;
    const std::string& dnaSequence() const;
    const std::string& nativeLinkerBefore() const;
    const std::string& nativeLinkerAfter() const;
    const std::string& refSeqId() const;
    const std::string& uniProtId() const;

    Pathway &pathway();
    void setDomainId(uint32_t);
    void setModuleId(uint32_t);
    void setPathway(const Pathway&);
    void setPathway(Pathway&&);
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
    Pathway m_pathway;
    std::string m_bioBrickId;
    std::string m_description;
    std::string m_dnaSeq;
    std::string m_nativeLinkerBefore;
    std::string m_nativeLinkerAfter;
    std::string m_refSeqId;
    std::string m_uniProtId;
};

template<bool full>
class DomainBaseType
{
public:
    typedef typename std::conditional<full, Domain, AbstractDomainType>::type type;
};
}

#endif
