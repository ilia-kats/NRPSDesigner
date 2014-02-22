#include "cds.h"
#include "origin.h"
#include "product.h"

#include <cstdlib>

#define CDS_NODE "cds"
#define ID_NODE "id"
#define GENENAME_NODE "genename"
#define DNASEQUENCE_NODE "dnasequence"
#define ORIGIN_NODE "origin"
#define PRODUCT_NODE "product"
#define DESCRIPTION_NODE "description"

using namespace nrps;

std::unordered_map<uint32_t, Cds*> Cds::s_cdss = std::unordered_map<uint32_t, Cds*>();

Cds::Cds(uint32_t id)
: m_id(id), m_origin(nullptr), m_product(nullptr)
{}

Cds::~Cds()
{
    s_cdss.erase(m_id);
}

uint32_t Cds::id() const
{
    return m_id;
}

void Cds::setId(uint32_t id)
{
    m_id = id;
}

const std::string& Cds::geneName() const
{
    return m_geneName;
}

void Cds::setGeneName(const std::string &name)
{
    m_geneName = name;
}

void Cds::setGeneName(std::string &&name)
{
    m_geneName = std::move(name);
}

const std::string& Cds::dnaSequence() const
{
    return m_dnaSequence;
}

void Cds::setDnaSequence(const std::string &seq)
{
    m_dnaSequence = seq;
}

void Cds::setDnaSequence(std::string &&seq)
{
    m_dnaSequence = std::move(seq);
}

const std::string& Cds::description() const
{
    return m_description;
}

void Cds::setDescription(const std::string &description)
{
    m_description = description;
}

void Cds::setDescription(std::string &&description)
{
    m_description = std::move(description);
}

Origin* Cds::origin() const
{
    return m_origin;
}

Origin* Cds::setOrigin(uint32_t id)
{
    m_origin = Origin::makeOrigin(id);
    return m_origin;
}

Product* Cds::product() const
{
    return m_product;
}

Product* Cds::setProduct(uint32_t id)
{
    m_product = Product::makeProduct(id);
    return m_product;
}

Cds* Cds::makeCds(uint32_t id)
{
    Cds *cds;
    if (s_cdss.count(id) == 0) {
        cds = new Cds(id);
        s_cdss.emplace(id, cds);
    } else {
        cds = s_cdss[id];
    }
    return cds;
}

#ifdef WITH_INTERNAL_XML
void Cds::toXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterStartElement(writer, BAD_CAST CDS_NODE);
    xmlTextWriterWriteElement(writer, BAD_CAST ID_NODE, BAD_CAST std::to_string(id()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST GENENAME_NODE, BAD_CAST geneName().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DNASEQUENCE_NODE, BAD_CAST dnaSequence().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST ORIGIN_NODE, BAD_CAST std::to_string(origin()->id()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST PRODUCT_NODE, BAD_CAST std::to_string(product()->id()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DESCRIPTION_NODE, BAD_CAST description().c_str());
    xmlTextWriterEndElement(writer);
}
#endif
