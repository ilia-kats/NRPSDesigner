#include "domain.h"

using namespace nrps;

std::unordered_map<uint32_t, std::shared_ptr<Pathway>> Domain::s_pathways = std::unordered_map<uint32_t, std::shared_ptr<Pathway>>();

Domain::Domain(DomainType type)
: AbstractDomainType(type)
{}

std::size_t Domain::hash() const
{
    return static_cast<std::size_t>(domainId());
}

bool Domain::full() const
{
    return true;
}

uint32_t Domain::domainId() const
{
    return m_domainId;
}

void Domain::setDomainId(uint32_t id)
{
    m_domainId = id;
}

uint32_t Domain::moduleId() const
{
    return m_moduleId;
}

void Domain::setModuleId(uint32_t id)
{
    m_moduleId = id;
}

const std::shared_ptr<Pathway>& Domain::pathway() const
{
    return s_pathways[m_pathway];
}

const std::shared_ptr<Pathway>& Domain::setPathway(uint32_t pathwayid)
{
    if (s_pathways.count(pathwayid) == 0) {
        s_pathways.emplace(pathwayid, std::shared_ptr<Pathway>(new Pathway(pathwayid)));
    }
    m_pathway = pathwayid;
    return s_pathways[m_pathway];
}

const std::string& Domain::bioBrickId() const
{
    return m_bioBrickId;
}

void Domain::setBioBrickId(const std::string &id)
{
    m_bioBrickId = id;
}

void Domain::setBioBrickId(std::string &&id)
{
    m_bioBrickId = std::move(id);
}

const std::string& Domain::description() const
{
    return m_description;
}

void Domain::setDescription(const std::string &description)
{
    m_description = description;
}

void Domain::setDescription(std::string &&description)
{
    m_description = std::move(description);
}

const std::string& Domain::dnaSequence() const
{
    return m_dnaSeq;
}

void Domain::setDnaSequence(const std::string &seq)
{
    m_dnaSeq = seq;
}

void Domain::setDnaSequence(std::string &&seq)
{
    m_dnaSeq = std::move(seq);
}

const std::string& Domain::nativeLinkerBefore() const
{
    return m_nativeLinkerBefore;
}

void Domain::setNativeLinkerBefore(const std::string &linker)
{
    m_nativeLinkerBefore = linker;
}

void Domain::setNativeLinkerBefore(std::string &&linker)
{
    m_nativeLinkerBefore = std::move(linker);
}

const std::string& Domain::nativeLinkerAfter() const
{
    return m_nativeLinkerAfter;
}

void Domain::setNativeLinkerAfter(const std::string &linker)
{
    m_nativeLinkerAfter = linker;
}

void Domain::setNativeLinkerAfter(std::string &&linker)
{
    m_nativeLinkerAfter = std::move(linker);
}

const std::string& Domain::refSeqId() const
{
    return m_refSeqId;
}

void Domain::setRefSeqId(const std::string &id)
{
    m_refSeqId = id;
}

void Domain::setRefSeqId(std::string &&id)
{
    m_refSeqId = std::move(id);
}

const std::string& Domain::uniProtId() const
{
    return m_uniProtId;
}

void Domain::setUniProtId(const std::string &id)
{
    m_uniProtId = id;
}

void Domain::setUniProtId(std::string &&id)
{
    m_uniProtId = std::move(id);
}

