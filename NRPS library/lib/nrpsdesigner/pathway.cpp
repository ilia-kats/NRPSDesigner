#include "pathway.h"

using namespace nrps;

Pathway::Pathway(uint32_t pathwayId, uint32_t taxId)
: m_pathwayId(pathwayId), m_taxId(taxId)
{}

Pathway::Pathway(uint32_t pathwayId, uint32_t taxId, const std::string &pathway, const std::string &linkout, const std::string &uniprotId, const std::string &norineId, const std::string &description)
: m_pathwayId(pathwayId), m_taxId(taxId), m_pathway(pathway), m_linkout(linkout), m_uniProtId(uniprotId), m_norineId(norineId), m_description(description)
{}

uint32_t Pathway::pathwayId() const
{
    return m_pathwayId;
}

void Pathway::setPathwayId(uint32_t id)
{
    m_pathwayId = id;
}

uint32_t Pathway::taxId() const
{
    return m_taxId;
}

void Pathway::setTaxId(uint32_t id)
{
    m_taxId = id;
}

const std::string& Pathway::pathway() const
{
    return m_pathway;
}

void Pathway::setPathway(const std::string &pathway)
{
    m_pathway = pathway;
}

void Pathway::setPathway(std::string &&pathway)
{
    m_pathway = std::move(pathway);
}

const std::string& Pathway::linkout() const
{
    return m_linkout;
}

void Pathway::setLinkout(const std::string &linkout)
{
    m_linkout = linkout;
}

void Pathway::setLinkout(std::string &&linkout)
{
    m_linkout = std::move(linkout);
}

const std::string& Pathway::uniProtId() const
{
    return m_uniProtId;
}

void Pathway::setUniProtId(const std::string &id)
{
    m_uniProtId = id;
}

void Pathway::setUniProtId(std::string &&id)
{
    m_uniProtId = std::move(id);
}

const std::string& Pathway::norineId() const
{
    return m_norineId;
}

void Pathway::setNorineId(const std::string &id)
{
    m_norineId = id;
}

void Pathway::setNorineId(std::string &&id)
{
    m_norineId = std::move(id);
}

const std::string& Pathway::description() const
{
    return m_description;
}

void Pathway::setDescription(const std::string &description)
{
    m_description = description;
}

void Pathway::setDescription(std::string &&description)
{
    m_description = std::move(description);
}
