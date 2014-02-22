#include "nrps.h"
#include "taxon.h"
#include "origin.h"
#include "product.h"
#include "cds.h"
#include "abstractdatabaseconnector.h"
#include "globals_internal.h"

#include <unordered_set>
#include <fstream>
#include <chrono>
#include <unistd.h>

#define NRPS_NODE "nrps"
#define CONSTRUCTDOMAIN_NODE "constructdomain"
#define START_NODE "start"
#define STOP_NODE "stop"
#define DOMAINS_NODE "domains"
#define ORIGINS_NODE "origins"
#define PRODUCTS_NODE "products"
#define CDS_NODE "cdss"
#define CONSTRUCTSEQUENCE_NODE "constructsequence"
#define INDIGOIDINETAGGED_ATTR "indigoidinetagged"

using namespace nrps;

#ifdef WITH_SBOL
class SbolBuilder
{
public:
    SbolBuilder(const Nrps &n)
    : m_nrps(n), m_start(1), m_uniqueId("_") {}

    DNAComponent* build(Document *doc)
    {
        if (doc == nullptr)
            return nullptr;
        m_doc = doc;
        m_uniqueId.append(std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        DNAComponent *nrps = createDNAComponent(m_doc, std::string("#nrps").append(m_uniqueId).c_str());
        setDNAComponentDisplayID(nrps, "NRPS");
        size_t start = 1;
        Domain *lastDomain = nullptr, *nextDomain = nullptr;
        for (size_t i = 0; i < m_nrps.size(); ++i) {
            m_start = 1;
            const auto &domain = m_nrps[i];
            if (i < m_nrps.size() - 1)
                nextDomain = m_nrps[i + 1].get();
            else
                nextDomain = nullptr;
            m_did = std::to_string(domain->id());
            m_id = std::to_string(i);
            m_dc = createDNAComponent(m_doc, std::string("#").append(m_id).append(m_uniqueId).c_str());
            setDNAComponentDisplayID(m_dc, std::string("_").append(m_did).c_str());
            setDNAComponentDescription(m_dc, domain->toString().c_str());
            SequenceAnnotation *rootsa = createSequenceAnnotation(m_doc, std::string("#sal1_").append(m_id).append(m_uniqueId).c_str());
            setSequenceAnnotationStart(rootsa, start);
            setSequenceAnnotationSubComponent(rootsa, m_dc);
            setSequenceAnnotationStrand(rootsa, STRAND_FORWARD);
            addSequenceAnnotation(nrps, rootsa);
            makeSubAnnotation(domain->dnaSequence(lastDomain, nextDomain), "d", "domain ");
            start += m_start;
            setSequenceAnnotationEnd(rootsa, start - 1);
            lastDomain = domain.get();
        }
        DNASequence *ds = createDNASequence(m_doc, std::string("#seq").append(m_uniqueId).c_str());
        setDNASequenceNucleotides(ds, m_seq.c_str());
        setDNAComponentSequence(nrps, ds);
        m_doc = nullptr;
        return nrps;
    }

    char* build()
    {
        Document *doc = createDocument();
        build(doc);
        char *ret = writeDocumentToString(doc);
        deleteDocument(doc);
        return ret;
    }

private:
    void makeSubAnnotation(const std::string &seq, const std::string &id, const std::string &desc)
    {
        DNAComponent *dc = createDNAComponent(m_doc, std::string("#").append(id).append(m_id).append(m_uniqueId).c_str());
        setDNAComponentDisplayID(dc, (id + m_did).c_str());
        setDNAComponentDescription(dc, (desc + m_did).c_str());
        SequenceAnnotation *sa = createSequenceAnnotation(m_doc, std::string("#sa").append(id).append(m_id).append(m_uniqueId).c_str());
        setSequenceAnnotationStart(sa, m_start);
        m_start += seq.size() - 1;
        setSequenceAnnotationEnd(sa, m_start);
        setSequenceAnnotationSubComponent(sa, dc);
        addSequenceAnnotation(m_dc, sa);
        m_seq.append(seq);
    }
    const Nrps &m_nrps;
    size_t m_start;
    std::string m_uniqueId;
    Document *m_doc;
    DNAComponent *m_dc;
    std::string m_did;
    std::string m_id;
    std::string m_seq;
};
#endif

Nrps::Nrps()
: std::vector<std::shared_ptr<Domain>>(), m_indigoidineTagged(false)
{}

bool Nrps::isIndigoidineTagged() const
{
    return m_indigoidineTagged;
}

void Nrps::setIndigoidineTagged(bool tagged)
{
    m_indigoidineTagged = tagged;
}

#ifdef WITH_INTERNAL_XML
std::string Nrps::toXml() const
{
    return nrps::toXml(static_cast<void (Nrps::*)(xmlTextWriterPtr)const>(&Nrps::toXml), this);
}

void Nrps::toXml(std::ostream &of) const
{
    nrps::toXml(of, static_cast<void (Nrps::*)(xmlTextWriterPtr)const>(&Nrps::toXml), this);
}

void Nrps::toXml(const std::string &file) const
{
    toXml(file.c_str());
}

void Nrps::toXml(const char *file) const
{
    nrps::toXml(file, static_cast<void (Nrps::*)(xmlTextWriterPtr)const>(&Nrps::toXml), this);
}

void Nrps::toXml(int fd) const
{
    nrps::toXml(fd, static_cast<void (Nrps::*)(xmlTextWriterPtr)const>(&Nrps::toXml), this);
}

void Nrps::toXml(xmlTextWriterPtr writer) const
{
    std::string sequence;
    std::unordered_set<Cds*> seenCds;
    std::unordered_set<Origin*> seenOrigins;
    std::unordered_set<Product*> seenProducts;
    std::vector<Origin*> originsToWrite;
    xmlTextWriterStartElement(writer, BAD_CAST NRPS_NODE);
    if (m_indigoidineTagged)
        xmlTextWriterWriteAttribute(writer, BAD_CAST INDIGOIDINETAGGED_ATTR, BAD_CAST INDIGOIDINETAGGED_ATTR);
    xmlTextWriterStartElement(writer, BAD_CAST DOMAINS_NODE);
    AbstractDatabaseConnector *dbconn = AbstractDatabaseConnector::getInstance();
    Domain *lastDomain = nullptr, *nextDomain = nullptr;
    for (auto it = cbegin(); it != end(); ++it) {
        if (it < end() - 1)
            nextDomain = (*(it + 1)).get();
        else
            nextDomain = nullptr;
        std::string currseq = (*it)->dnaSequence(lastDomain, nextDomain);
        lastDomain = it->get();
        if ((*it)->cds() != nullptr)
            seenCds.insert((*it)->cds());
        xmlTextWriterStartElement(writer, BAD_CAST CONSTRUCTDOMAIN_NODE);
        xmlTextWriterWriteElement(writer, BAD_CAST START_NODE, BAD_CAST std::to_string(sequence.size()).c_str());
        xmlTextWriterWriteElement(writer, BAD_CAST STOP_NODE, BAD_CAST std::to_string(sequence.size() + currseq.size() - 1).c_str());
        (*it)->toXml(writer);
        xmlTextWriterEndElement(writer);
        sequence.append(currseq);
    }
    xmlTextWriterEndElement(writer);
    xmlTextWriterWriteElement(writer, BAD_CAST CONSTRUCTSEQUENCE_NODE, BAD_CAST sequence.c_str());
    xmlTextWriterStartElement(writer, BAD_CAST CDS_NODE);
    for (Cds *cds : seenCds) {
        if (cds->origin() != nullptr)
            originsToWrite.push_back(cds->origin());
        if (cds->product() != nullptr)
            seenProducts.insert(cds->product());
        cds->toXml(writer);
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
}
#endif

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

DNAComponent* Nrps::toSbol(Document *doc) const
{
    SbolBuilder builder(*this);
    return builder.build(doc);
}
#endif
