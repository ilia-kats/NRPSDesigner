#include "nrps.h"
#include "taxon.h"
#include "origin.h"
#include "product.h"
#include "abstractdatabaseconnector.h"

#include <unordered_set>

#include <unistd.h>

#define NRPS_NODE "nrps"
#define DOMAINS_NODE "domains"
#define ORIGINS_NODE "origins"
#define PRODUCTS_NODE "products"

using namespace nrps;

Nrps::Nrps(const std::vector<Monomer> &nrp)
: std::vector<std::shared_ptr<Domain>>(), m_nrp(nrp)
{}

std::string Nrps::toXml() const
{
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    toXml(writer);
    std::string xml((const char*)buf->content);
    xmlFreeTextWriter(writer);
    xmlBufferFree(buf);
    return xml;
}

void Nrps::toXml(std::ostream &of) const
{
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    toXml(writer);
    of << (const char*)buf->content;
    xmlFreeTextWriter(writer);
    xmlBufferFree(buf);
}

void Nrps::toXml(const std::string &file) const
{
    toXml(file.c_str());
}

void Nrps::toXml(const char *file) const
{
    xmlTextWriterPtr writer = xmlNewTextWriterFilename(file, 0);
    toXml(writer);
    xmlFreeTextWriter(writer);
}

void Nrps::toXml(int fd) const
{
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    toXml(writer);
    write(fd, buf->content, buf->use);
    xmlFreeTextWriter(writer);
    xmlBufferFree(buf);
}

void Nrps::toXml(xmlTextWriterPtr writer) const
{
    std::unordered_set<Origin*> seenOrigins;
    std::unordered_set<Product*> seenProducts;
    std::vector<Origin*> originsToWrite;
    xmlTextWriterSetIndent(writer, 1);
    xmlTextWriterSetIndentString(writer, BAD_CAST "    ");
    xmlTextWriterStartDocument(writer, nullptr, "UTF-8", nullptr);
    xmlTextWriterStartElement(writer, BAD_CAST NRPS_NODE);
    xmlTextWriterStartElement(writer, BAD_CAST DOMAINS_NODE);
    AbstractDatabaseConnector *dbconn = AbstractDatabaseConnector::getInstance();
    for (const auto &domain : *this) {
        dbconn->fillDomain(domain);
        originsToWrite.push_back(domain->origin());
        seenProducts.insert(domain->product());
        domain->toXml(writer);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterStartElement(writer, BAD_CAST ORIGINS_NODE);
    for (Origin *ori : originsToWrite) {
        if (seenOrigins.count(ori) == 0) {
            while (ori != nullptr && seenOrigins.count(ori) == 0) {
                ori->toXml(writer);
                seenOrigins.insert(ori);
                ori = ori->parent();
            }
        }
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterStartElement(writer, BAD_CAST PRODUCTS_NODE);
    for (Product *pr : seenProducts) {
        pr->toXml(writer);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterEndElement(writer);
    xmlTextWriterEndDocument(writer);
}
