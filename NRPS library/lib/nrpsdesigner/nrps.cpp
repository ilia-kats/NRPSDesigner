#include "nrps.h"
#include "taxon.h"
#include "origin.h"
#include "product.h"
#include "abstractdatabaseconnector.h"

#include <unordered_set>
#include <fstream>
#include <unistd.h>

#ifdef WITH_SBOL
extern "C" {
#include <sbol.h>
}
#endif

#define NRPS_NODE "nrps"
#define DOMAINS_NODE "domains"
#define ORIGINS_NODE "origins"
#define PRODUCTS_NODE "products"
#define INDIGOIDINETAGGED_ATTR "indigoidinetagged"

using namespace nrps;

#ifdef WITH_SBOL
class SbolBuilder
{
public:
    SbolBuilder(const Nrps &n)
    : m_nrps(n), m_start(1) {}
    char* build()
    {
        m_doc = createDocument();
        DNAComponent *nrps = createDNAComponent(m_doc, "#nrps");
        setDNAComponentDisplayID(nrps, "NRPS");
        for (size_t i = 0; i < m_nrps.size(); ++i) {
            const auto &domain = m_nrps[i];
            m_did = std::to_string(domain->id());
            m_id = std::to_string(i);
            m_dc = createDNAComponent(m_doc, std::string("#").append(m_id).c_str());
            setDNAComponentDisplayID(m_dc, std::string("_").append(m_did).c_str());
            setDNAComponentDescription(m_dc, domain->toString().c_str());
            SequenceAnnotation *rootsa = createSequenceAnnotation(m_doc, std::string("#sal1_").append(m_id).c_str());
            setSequenceAnnotationStart(rootsa, m_start);
            setSequenceAnnotationSubComponent(rootsa, m_dc);
            setSequenceAnnotationStrand(rootsa, STRAND_FORWARD);
            addSequenceAnnotation(nrps, rootsa);
            if (!domain->determinedLinkerBefore().empty())
                makeSubAnnotation(domain->determinedLinkerBefore(), "lb", "linker before domain ");
            makeSubAnnotation(domain->dnaSequenceDefined().empty() ? domain->dnaSequencePfam() : domain->dnaSequenceDefined(), "d", "domain ");
            if (!domain->determinedLinkerAfter().empty())
                makeSubAnnotation(domain->determinedLinkerAfter(), "la", "linker after domain ");
            setSequenceAnnotationEnd(rootsa, m_start++);
        }
        DNASequence *ds = createDNASequence(m_doc, "#seq");
        setDNASequenceNucleotides(ds, m_seq.c_str());
        setDNAComponentSequence(nrps, ds);
        char *ret = writeDocumentToString(m_doc);
        deleteDocument(m_doc);
        return ret;
    }

private:
    void makeSubAnnotation(const std::string &seq, const std::string &id, const std::string &desc)
    {
        DNAComponent *dc = createDNAComponent(m_doc, std::string("#").append(id).append(m_id).c_str());
        setDNAComponentDisplayID(dc, (id + m_did).c_str());
        setDNAComponentDescription(dc, (desc + m_did).c_str());
        SequenceAnnotation *sa = createSequenceAnnotation(m_doc, std::string("#sa").append(id).append(m_id).c_str());
        setSequenceAnnotationStart(sa, m_start);
        m_start += seq.size();
        setSequenceAnnotationEnd(sa, m_start);
        setSequenceAnnotationSubComponent(sa, dc);
        addSequenceAnnotation(m_dc, sa);
        m_seq.append(seq);
    }
    const Nrps &m_nrps;
    size_t m_start;
    Document *m_doc;
    DNAComponent *m_dc;
    std::string m_did;
    std::string m_id;
    std::string m_seq;
};
#endif

Nrps::Nrps(const std::vector<Monomer> &nrp)
: std::vector<std::shared_ptr<Domain>>(), m_nrp(nrp), m_indigoidineTagged(false)
{}

bool Nrps::isIndigoidineTagged() const
{
    return m_indigoidineTagged;
}

void Nrps::setIndigoidineTagged(bool tagged)
{
    m_indigoidineTagged = tagged;
}

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
    if (m_indigoidineTagged)
        xmlTextWriterWriteAttribute(writer, BAD_CAST INDIGOIDINETAGGED_ATTR, BAD_CAST INDIGOIDINETAGGED_ATTR);
    xmlTextWriterStartElement(writer, BAD_CAST DOMAINS_NODE);
    AbstractDatabaseConnector *dbconn = AbstractDatabaseConnector::getInstance();
    for (const auto &domain : *this) {
        if (domain->origin() != nullptr)
            originsToWrite.push_back(domain->origin());
        if (domain->product() != nullptr)
            seenProducts.insert(domain->product());
        domain->toXml(writer);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterStartElement(writer, BAD_CAST ORIGINS_NODE);
    for (Origin *ori : originsToWrite) {
        if (!seenOrigins.count(ori)) {
            while (ori != nullptr && !seenOrigins.count(ori)) {
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

#ifdef WITH_SBOL

std::string Nrps::toSbol() const
{
    SbolBuilder builder(*this);
    char *sbol = builder.build();
    std::string seq(sbol);
    std::free(sbol);
    return seq;
}

void Nrps::toSbol(std::ostream &of) const
{
    SbolBuilder builder(*this);
    char *sbol = builder.build();
    of << sbol;
    std::free(sbol);
}

void Nrps::toSbol(const std::string &file) const
{
    toSbol(file.c_str());
}

void Nrps::toSbol(const char *file) const
{
    std::ofstream of(file);
    toSbol(of);
    of.close();
}

void Nrps::toSbol(int fd) const
{
    SbolBuilder builder(*this);
    char *sbol = builder.build();
    write(fd, sbol, strlen(sbol));
    std::free(sbol);
}

#endif
