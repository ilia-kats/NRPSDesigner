#include "origin.h"

#include <cstdlib>

#define ORIGIN_NODE "origin"
#define ID_NODE "id"
#define SOURCETYPE_NODE "sourcetype"
#define SOURCE_NODE "source"
#define SPECIES_NODE "species"
#define DESCRIPTION_NODE "description"
#define TAXID_NODE "taxid"
#define PARENT_NODE "parent"

using namespace nrps;

std::unordered_map<uint32_t, Origin*> Origin::s_origins = std::unordered_map<uint32_t, Origin*>();

Origin::Origin(uint32_t id)
: m_id(id), m_taxid(0), m_parent(nullptr), m_sourceType(OriginSourceType::DNA)
{}

Origin::~Origin()
{
    if (m_parent != nullptr)
        m_parent->m_children.erase(this);
    for (Origin *child : m_children)
        delete child;
    s_origins.erase(m_id);
}

uint32_t Origin::id() const
{
    return m_id;
}

void Origin::setId(uint32_t id)
{
    m_id = id;
}

uint32_t Origin::taxId() const
{
    return m_taxid;
}

OriginSourceType Origin::sourceType() const
{
    return m_sourceType;
}

void Origin::setSourceType(OriginSourceType srctype)
{
    m_sourceType = srctype;
    if (srctype == OriginSourceType::Species && !m_source.empty()) {
        setTaxId(std::atoi(m_source.c_str()));
    }
}

const std::string& Origin::source() const
{
    return m_source;
}

void Origin::setSource(const std::string &src)
{
    m_source = src;
    if (m_sourceType == OriginSourceType::Species && !m_source.empty()) {
        setTaxId(std::atoi(m_source.c_str()));
    }
}

void Origin::setSource(std::string &&src)
{
    m_source = std::move(src);
    if (m_sourceType == OriginSourceType::Species && !m_source.empty()) {
        setTaxId(std::atoi(m_source.c_str()));
    }
}

const std::string& Origin::species() const
{
    return m_species;
}

void Origin::setSpecies(const std::string &species)
{
    m_species = species;
}

void Origin::setSpecies(std::string &&species)
{
    m_species = std::move(species);
}

const std::string& Origin::description() const
{
    return m_description;
}

void Origin::setDescription(const std::string &description)
{
    m_description = description;
}

void Origin::setDescription(std::string &&description)
{
    m_description = std::move(description);
}

Origin* Origin::parent() const
{
    return m_parent;
}

void Origin::setParent(Origin* parent)
{
    m_parent = parent;
    m_parent->m_children.insert(this);
    setTaxId(m_parent->m_taxid);
}

Origin* Origin::makeOrigin(uint32_t id)
{
    Origin *ori;
    if (s_origins.count(id) == 0) {
        ori = new Origin(id);
        s_origins.emplace(id, ori);
    } else {
        ori = s_origins[id];
    }
    return ori;
}

void Origin::setTaxId(uint32_t taxid)
{
    m_taxid = taxid;
    for (Origin *child : m_children) {
        child->setTaxId(taxid);
    }
}

#ifdef WITH_INTERNAL_XML
void Origin::toXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterStartElement(writer, BAD_CAST ORIGIN_NODE);
    xmlTextWriterWriteElement(writer, BAD_CAST ID_NODE, BAD_CAST std::to_string(id()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST SOURCETYPE_NODE, BAD_CAST nrps::toString(sourceType()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST SOURCE_NODE, BAD_CAST source().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST SPECIES_NODE, BAD_CAST species().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DESCRIPTION_NODE, BAD_CAST description().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST TAXID_NODE, BAD_CAST std::to_string(taxId()).c_str());
    if (parent() != nullptr)
        xmlTextWriterWriteElement(writer, BAD_CAST PARENT_NODE, BAD_CAST std::to_string(parent()->id()).c_str());
    xmlTextWriterEndElement(writer);
}
#endif

std::string Origin::toString() const
{
    std::string retVal = nrps::toString(sourceType());
    retVal.append("; taxID ").append(std::to_string(taxId())).append(": ").append(species());
    if (sourceType() == OriginSourceType::DNA)
        retVal.append("; source: ").append(source());
    return retVal;
}
