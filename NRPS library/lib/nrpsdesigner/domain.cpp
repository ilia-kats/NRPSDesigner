#include "domain.h"
#include "origin.h"
#include "product.h"
#include "cds.h"

#define NRPS_NODE "nrps"
#define DOMAIN_NODE "domain"
#define TYPE_NODE "type"
#define ID_NODE "id"
#define MODULE_NODE "module"
#define DESCRIPTION_NODE "description"
#define PFAMSTART_NODE "pfamstart"
#define PFAMSTOP_NODE "pfamstop"
#define DEFINEDSTART_NODE "definedstart"
#define DEFINEDSTOP_NODE "definedstop"
#define PFAMLINKERSTART_NODE "pfamlinkerstart"
#define PFAMLINKERSTOP_NODE "pfamlinkerstop"
#define DEFINEDLINKERSTART_NODE "definedlinkerstart"
#define DEFINEDLINKERSTOP_NODE "definedlinkerstop"
#define EXPERIMENTALWORKED_NODE "workedwithnextdomains"
#define STOP_NODE "stop"
#define NEXTSTART_NODE "nextdomainstart"
#define CDSID_NODE "cdsid"

using namespace nrps;

Domain::Domain(DomainType t, uint32_t id)
: m_type(t), m_id(id), m_module(0), m_cds(nullptr)
{}

std::size_t Domain::hash() const
{
    return static_cast<std::size_t>(id());
}

DomainType Domain::type() const
{
    return m_type;
}

uint32_t Domain::id() const
{
    return m_id;
}

void Domain::setId(uint32_t id)
{
    m_id = id;
}

uint32_t Domain::module() const
{
    return m_module;
}

void Domain::setModule(uint32_t id)
{
    m_module = id;
}

Origin* Domain::origin() const
{
    if (m_cds != nullptr)
        return m_cds->origin();
    else
        return nullptr;
}

Origin* Domain::setOrigin(uint32_t id)
{
    if (m_cds != nullptr)
        return m_cds->setOrigin(id);
    else
        return nullptr;
}

Product* Domain::product() const
{
    if (m_cds != nullptr)
        return m_cds->product();
    else
        return nullptr;
}

Product* Domain::setProduct(uint32_t id)
{
    if (m_cds != nullptr)
        return m_cds->setProduct(id);
    else
        return nullptr;
}

Cds* Domain::cds() const
{
    return m_cds;
}

Cds* Domain::setCds(uint32_t id)
{
    m_cds = Cds::makeCds(id);
    return m_cds;
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

std::string Domain::dnaSequencePfam() const
{
    if (m_cds != nullptr && m_pfamStop <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(m_pfamStart, m_pfamStop - m_pfamStart + 1);
    else
        return std::string();
}

std::string Domain::dnaSequenceDefined() const
{
    if (m_cds != nullptr && m_definedStop <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(m_definedStart, m_definedStart - m_definedStop + 1);
    else
        return std::string();
}

std::string Domain::pfamLinkerBefore() const
{
    if (m_cds != nullptr && m_pfamLinkerStart <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(m_pfamLinkerStart, m_pfamStart - m_pfamLinkerStart);
    else
        return std::string();
}

std::string Domain::pfamLinkerAfter() const
{
    if (m_cds != nullptr && m_pfamLinkerStop <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(m_pfamStop + 1, m_pfamLinkerStop - m_pfamStop);
    else
        return std::string();
}

std::string Domain::definedLinkerBefore() const
{
    if (m_cds != nullptr && m_definedLinkerStart <= m_cds->dnaSequence().size() && m_definedStart <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(m_definedLinkerStart, m_definedStart - m_definedLinkerStart);
    else
        return std::string();
}

std::string Domain::definedLinkerAfter() const
{
    if (m_cds != nullptr && m_definedLinkerStop <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(m_definedStop + 1, m_definedLinkerStop - m_definedStop);
    else
        return std::string();
}


std::string Domain::dnaSequence(Domain *prev, Domain *next) const
{
    uint32_t start, stop;
    if (prev != nullptr && prev->worksWithNextDomain(id()))
        start = prev->workingDomainBorders(id()).next_start;
    else
        start = definedLinkerStart();
    if (next != nullptr && worksWithNextDomain(next->id()))
        stop = workingDomainBorders(next->id()).stop;
    else
        stop = definedStop();
    if (m_cds != nullptr && start <= m_cds->dnaSequence().size() && stop <= m_cds->dnaSequence().size())
        return m_cds->dnaSequence().substr(start, stop - start + 1);
    else
        return std::string();
}

uint32_t Domain::pfamStart() const
{
    return m_pfamStart;
}

void Domain::setPfamStart(uint32_t start)
{
    m_pfamStart = start;
}

uint32_t Domain::pfamStop() const
{
    return m_pfamStop;
}

void Domain::setPfamStop(uint32_t stop)
{
    m_pfamStop = stop;
}

uint32_t Domain::definedStart() const
{
    return m_definedStart;
}

void Domain::setDefinedStart(uint32_t start)
{
    m_definedStart = start;
}

uint32_t Domain::definedStop() const
{
    return m_definedStop;
}

void Domain::setDefinedStop(uint32_t stop)
{
    m_definedStop = stop;
}

uint32_t Domain::pfamLinkerStart() const
{
    return m_pfamLinkerStart;
}

void Domain::setPfamLinkerStart(uint32_t start)
{
    m_pfamLinkerStart = start;
}

uint32_t Domain::pfamLinkerStop() const
{
    return m_pfamLinkerStop;
}

void Domain::setPfamLinkerStop(uint32_t stop)
{
    m_pfamLinkerStop = stop;
}

uint32_t Domain::definedLinkerStart() const
{
    return m_definedLinkerStart;
}

void Domain::setDefinedLinkerStart(uint32_t start)
{
    m_definedLinkerStart = start;
}

uint32_t Domain::definedLinkerStop() const
{
    return m_definedLinkerStop;
}

void Domain::setDefinedLinkerStop(uint32_t stop)
{
    m_definedLinkerStop = stop;
}

void Domain::addWorkingNextDomain(uint32_t did, Domain::workingBorders borders)
{
    m_workingNextDomains.emplace(did, borders);
}

void Domain::removeWorkingNextDomain(uint32_t did)
{
    m_workingNextDomains.erase(did);
}

const std::unordered_map<uint32_t, Domain::workingBorders>& Domain::workingNextDomains() const
{
    return m_workingNextDomains;
}

bool Domain::worksWithNextDomain(uint32_t id) const
{
    return m_workingNextDomains.count(id) > 0;
}

const Domain::workingBorders& Domain::workingDomainBorders(uint32_t id) const
{
    return m_workingNextDomains.at(id);
}

#ifdef WITH_INTERNAL_XML
void Domain::toXml(xmlTextWriterPtr writer) const
{
    startXml(writer);
    writeXml(writer);
    endXml(writer);
}

void Domain::startXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterStartElement(writer, BAD_CAST DOMAIN_NODE);
}

void Domain::writeXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterWriteElement(writer, BAD_CAST ID_NODE, BAD_CAST std::to_string(id()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST TYPE_NODE, BAD_CAST nrps::toString(type()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST MODULE_NODE, BAD_CAST std::to_string(module()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DESCRIPTION_NODE, BAD_CAST description().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST PFAMSTART_NODE, BAD_CAST std::to_string(pfamStart()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST PFAMSTOP_NODE, BAD_CAST std::to_string(pfamStop()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DEFINEDSTART_NODE, BAD_CAST std::to_string(definedStart()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DEFINEDSTOP_NODE, BAD_CAST std::to_string(definedStop()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST PFAMLINKERSTART_NODE, BAD_CAST std::to_string(pfamLinkerStart()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST PFAMLINKERSTOP_NODE, BAD_CAST std::to_string(pfamLinkerStop()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DEFINEDLINKERSTART_NODE, BAD_CAST std::to_string(definedLinkerStart()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DEFINEDLINKERSTOP_NODE, BAD_CAST std::to_string(definedLinkerStop()).c_str());

    if (!m_workingNextDomains.empty()) {
        xmlTextWriterStartElement(writer, BAD_CAST EXPERIMENTALWORKED_NODE);
        for (const auto &it : m_workingNextDomains) {
            xmlTextWriterStartElement(writer, BAD_CAST DOMAIN_NODE);
            xmlTextWriterWriteElement(writer, BAD_CAST ID_NODE, BAD_CAST std::to_string(it.first).c_str());
            xmlTextWriterWriteElement(writer, BAD_CAST STOP_NODE, BAD_CAST std::to_string(it.second.stop).c_str());
            xmlTextWriterWriteElement(writer, BAD_CAST NEXTSTART_NODE, BAD_CAST std::to_string(it.second.next_start).c_str());
            xmlTextWriterEndElement(writer);
        }
        xmlTextWriterEndElement(writer);
    }

    if (cds() != nullptr)
        xmlTextWriterWriteElement(writer, BAD_CAST CDSID_NODE, BAD_CAST std::to_string(cds()->id()).c_str());
}

void Domain::endXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterEndElement(writer);
}
#endif

std::string Domain::toString() const
{
    std::string retVal("Domain of type ");
    retVal.append(nrps::toString(type())).append(";");
    if (product() != nullptr)
        retVal.append("\nmodule ").append(std::to_string(module())).append(" of ").append(product()->name()).append(" pathway");
    if (origin() != nullptr)
        retVal.append("\nfrom ").append(origin()->toString());
    return retVal;
}
