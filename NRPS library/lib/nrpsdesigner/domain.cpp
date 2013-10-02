#include "domain.h"
#include "origin.h"
#include "product.h"

#define NRPS_NODE "nrps"
#define DOMAIN_NODE "domain"
#define TYPE_NODE "type"
#define ID_NODE "id"
#define MODULE_NODE "module"
#define DESCRIPTION_NODE "description"
#define GENENAME_NODE "genename"
#define GENEDESCRIPTION_NODE "genedescription"
#define SEQUENCEPFAM_NODE "sequencepfam"
#define SEQUENCEDEFINED_NODE "sequencedefined"
#define NATIVEPFAMLINKERBEFORE_NODE "nativepfamlinkerbefore"
#define NATIVEPFAMLINKERAFTER_NODE "nativepfamlinkerafter"
#define NATIVEDEFINEDLINKERBEFORE_NODE "nativedefinedlinkerbefore"
#define NATIVEDEFINEDLINKERAFTER_NODE "nativedefinedlinkerafter"
#define DETERMINEDLINKERBEFORE_NODE "determinedlinkerbefore"
#define DETERMINEDLINKERAFTER_NODE "determinedlinkerafter"
#define ORIGINID_NODE "originid"
#define PRODUCTID_NODE "productid"

using namespace nrps;

Domain::Domain(DomainType t, uint32_t id)
: m_type(t), m_id(id), m_module(0), m_origin(nullptr), m_product(nullptr)
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
    return m_origin;
}

Origin* Domain::setOrigin(uint32_t id)
{
    m_origin = Origin::makeOrigin(id);
    return m_origin;
}

Product* Domain::product() const
{
    return m_product;
}

Product* Domain::setProduct(uint32_t id)
{
    m_product = Product::makeProduct(id);
    return m_product;
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

const std::string& Domain::dnaSequencePfam() const
{
    return m_dnaSeqPfam;
}

void Domain::setDnaSequencePfam(const std::string &seq)
{
    m_dnaSeqPfam = seq;
}

void Domain::setDnaSequencePfam(std::string &&seq)
{
    m_dnaSeqPfam = std::move(seq);
}

const std::string& Domain::dnaSequenceDefined() const
{
    return m_dnaSeqDefined;
}

void Domain::setDnaSequenceDefined(const std::string &seq)
{
    m_dnaSeqDefined = seq;
}

void Domain::setDnaSequenceDefined(std::string &&seq)
{
    m_dnaSeqDefined = std::move(seq);
}

const std::string& Domain::geneDescription() const
{
    return m_geneDescription;
}

void Domain::setGeneDescription(const std::string &description)
{
    m_geneDescription = description;
}

void Domain::setGeneDescription(std::string &&description)
{
    m_geneDescription = std::move(description);
}

const std::string& Domain::geneName() const
{
    return m_geneName;
}

void Domain::setGeneName(const std::string &name)
{
    m_geneName = name;
}

void Domain::setGeneName(std::string &&name)
{
    m_geneName = std::move(name);
}

const std::string& Domain::nativePfamLinkerBefore() const
{
    return m_nativePfamLinkerBefore;
}

void Domain::setNativePfamLinkerBefore(const std::string &linker)
{
    m_nativePfamLinkerBefore = linker;
}

void Domain::setNativePfamLinkerBefore(std::string &&linker)
{
    m_nativePfamLinkerBefore = std::move(linker);
}

const std::string& Domain::nativePfamLinkerAfter() const
{
    return m_nativePfamLinkerAfter;
}

void Domain::setNativePfamLinkerAfter(const std::string &linker)
{
    m_nativePfamLinkerAfter = linker;
}

void Domain::setNativePfamLinkerAfter(std::string &&linker)
{
    m_nativePfamLinkerAfter = std::move(linker);
}

const std::string& Domain::nativeDefinedLinkerBefore() const
{
    return m_nativeDefinedLinkerBefore;
}

void Domain::setNativeDefinedLinkerBefore(const std::string &linker)
{
    m_nativeDefinedLinkerBefore = linker;
}

void Domain::setNativeDefinedLinkerBefore(std::string &&linker)
{
    m_nativeDefinedLinkerBefore = std::move(linker);
}

const std::string& Domain::nativeDefinedLinkerAfter() const
{
    return m_nativeDefinedLinkerAfter;
}

void Domain::setNativeDefinedLinkerAfter(const std::string &linker)
{
    m_nativeDefinedLinkerAfter = linker;
}

void Domain::setNativeDefinedLinkerAfter(std::string &&linker)
{
    m_nativeDefinedLinkerAfter = std::move(linker);
}

const std::string& Domain::determinedLinkerBefore() const
{
    return m_determinedLinkerBefore;
}

void Domain::setDeterminedLinkerBefore(const std::string &linker)
{
    m_determinedLinkerBefore = linker;
}

void Domain::setDeterminedLinkerBefore(std::string &&linker)
{
    m_determinedLinkerBefore = std::move(linker);
}

const std::string& Domain::determinedLinkerAfter() const
{
    return m_determinedLinkerAfter;
}

void Domain::setDeterminedLinkerAfter(const std::string &linker)
{
    m_determinedLinkerAfter = linker;
}

void Domain::setDeterminedLinkerAfter(std::string &&linker)
{
    m_determinedLinkerAfter = std::move(linker);
}

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
    xmlTextWriterWriteElement(writer, BAD_CAST GENEDESCRIPTION_NODE, BAD_CAST geneDescription().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST GENENAME_NODE, BAD_CAST geneName().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST SEQUENCEPFAM_NODE, BAD_CAST dnaSequencePfam().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST SEQUENCEDEFINED_NODE, BAD_CAST dnaSequenceDefined().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST NATIVEPFAMLINKERBEFORE_NODE, BAD_CAST nativePfamLinkerBefore().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST NATIVEPFAMLINKERAFTER_NODE, BAD_CAST nativePfamLinkerAfter().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST NATIVEDEFINEDLINKERBEFORE_NODE, BAD_CAST nativeDefinedLinkerBefore().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST NATIVEDEFINEDLINKERAFTER_NODE, BAD_CAST nativeDefinedLinkerAfter().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DETERMINEDLINKERBEFORE_NODE, BAD_CAST determinedLinkerBefore().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DETERMINEDLINKERAFTER_NODE, BAD_CAST determinedLinkerAfter().c_str());

    if (origin() != nullptr)
        xmlTextWriterWriteElement(writer, BAD_CAST ORIGINID_NODE, BAD_CAST std::to_string(origin()->id()).c_str());
    if (product() != nullptr)
        xmlTextWriterWriteElement(writer, BAD_CAST PRODUCTID_NODE, BAD_CAST std::to_string(product()->id()).c_str());
}

void Domain::endXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterEndElement(writer);
}

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
