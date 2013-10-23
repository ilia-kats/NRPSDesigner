#include "nrpslibrary.h"
#include "abstractdatabaseconnector.h"
#include "nrpsbuilder.h"
#include "globals_internal.h"

#include <ios>
#include <iterator>
#include <stdexcept>

#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>

#include <libxml/xmlreader.h>

#define NRP_NODE "nrp"
#define MONOMERS_NODE "monomers"
#define MONOMER_NODE "monomer"
#define NRPS_NODE "nrps"
#define DOMAIN_NODE "domain"
#define ID_NODE "id"
#define NRPSLIBRARY_NODE "nrpslibrary"


using namespace nrps;

int skip_whitespace(xmlTextReaderPtr reader, bool skip_significant = true)
{
    int ret;
    while ((ret = xmlTextReaderRead(reader)) == 1) {
        int type = xmlTextReaderNodeType(reader);
        if (type == XML_READER_TYPE_WHITESPACE || skip_significant && type == XML_READER_TYPE_SIGNIFICANT_WHITESPACE)
            break;
    }
    ret = xmlTextReaderRead(reader);
    return ret;
}

NrpsLibrary::NrpsLibrary(const char *file)
{
    fromFile(file);
}

NrpsLibrary::NrpsLibrary(const std::string &file)
{
    fromFile(file);
}

NrpsLibrary::NrpsLibrary(int fd)
{
    fromFile(fd);
}

NrpsLibrary::NrpsLibrary(std::istream &stream)
{
    fromFile(stream);
}

NrpsLibrary::NrpsLibrary()
{}

NrpsLibrary::~NrpsLibrary()
{}

void NrpsLibrary::fromFile(const char *file)
{
    int fd = open(file, O_RDONLY);
    fromFile(fd);
    close(fd);
}

void NrpsLibrary::fromFile(const std::string &file)
{
    fromFile(file.c_str());
}

void NrpsLibrary::fromFile(int fd)
{
    xmlDocPtr doc = xmlReadFd(fd, nullptr, nullptr, XML_PARSE_RECOVER | XML_PARSE_NONET | XML_PARSE_COMPACT);
    readXml(doc);
    xmlFreeDoc(doc);
}

void NrpsLibrary::fromFile(std::istream &stream)
{
    std::string xml((std::istreambuf_iterator<std::string::value_type>(stream)), std::istreambuf_iterator<std::string::value_type>());
    xmlDocPtr doc = xmlReadMemory(xml.c_str(), xml.size(), nullptr, nullptr, XML_PARSE_RECOVER | XML_PARSE_NONET | XML_PARSE_COMPACT);
    readXml(doc);
    xmlFreeDoc(doc);
}

void NrpsLibrary::makeLibrary()
{
    if (m_nrp.empty() || m_nrps.empty())
        throw std::invalid_argument("No scaffold to build library upon.");
    std::vector<Monomer> scaffoldnrp;
    int switchedpos = -1;
    int i = 0;
    for (const auto &monomers : m_nrp) {
        scaffoldnrp.push_back(monomers.front());
        if (monomers.size() > 1)
            switchedpos = i;
        ++i;
    }
    if (switchedpos == -1)
        throw std::invalid_argument("No variants for the library defined.");
    std::vector<Monomer> nrp(scaffoldnrp);
    NrpsBuilder builder;
    for (int i = 1; i < m_nrp[switchedpos].size(); ++i) {
        nrp[switchedpos] = m_nrp[switchedpos][i];
        builder.setScaffold(&scaffoldnrp, &m_nrps);
        m_library.push_back(builder.build(nrp));
    }
}

#ifdef WITH_INTERNAL_XML
std::string NrpsLibrary::toXml() const
{
    return nrps::toXml(static_cast<void (NrpsLibrary::*)(xmlTextWriterPtr)const>(&NrpsLibrary::toXml), this);
}

void NrpsLibrary::toXml(std::ostream &of) const
{
    nrps::toXml(of, static_cast<void (NrpsLibrary::*)(xmlTextWriterPtr)const>(&NrpsLibrary::toXml), this);
}

void NrpsLibrary::toXml(const std::string &file) const
{
    toXml(file.c_str());
}

void NrpsLibrary::toXml(const char *file) const
{
    nrps::toXml(file, static_cast<void (NrpsLibrary::*)(xmlTextWriterPtr)const>(&NrpsLibrary::toXml), this);
}

void NrpsLibrary::toXml(int fd) const
{
    nrps::toXml(fd, static_cast<void (NrpsLibrary::*)(xmlTextWriterPtr)const>(&NrpsLibrary::toXml), this);
}

void NrpsLibrary::toXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterStartElement(writer, BAD_CAST NRPSLIBRARY_NODE);
    for (const auto &nrps : m_library)
        nrps.toXml(writer);
    xmlTextWriterEndElement(writer);
}
#endif

#ifdef WITH_SBOL
std::string NrpsLibrary::toSbol() const
{
    Document *doc = createDocument();
    toSbol(doc);
    char *sbol = writeDocumentToString(doc);
    deleteDocument(doc);
    std::string seq(sbol);
    std::free(sbol);
    return seq;
}

void NrpsLibrary::toSbol(std::ostream &of) const
{
    of << toSbol();
}

void NrpsLibrary::toSbol(const std::string &file) const
{
    toSbol(file.c_str());
}

void NrpsLibrary::toSbol(const char *file) const
{
    std::ofstream of(file);
    toSbol(of);
    of.close();
}

void NrpsLibrary::toSbol(int fd) const
{
    std::string sbol = toSbol();
    write(fd, sbol.c_str(), sbol.size());
}

Collection* NrpsLibrary::toSbol(Document *doc) const
{
    std::string uniqueId("_");
    uniqueId.append(std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    Collection *col = createCollection(doc, std::string("#library").append(uniqueId).c_str());
    for (const auto &nrps : m_library) {
        addDNAComponentToCollection(col, nrps.toSbol(doc));
    }
    return col;
}
#endif

void NrpsLibrary::readXml(xmlDocPtr doc)
{
    AbstractDatabaseConnector *dbConn = AbstractDatabaseConnector::getInstance();
    std::vector<uint32_t> domains;
    xmlNodePtr node = xmlDocGetRootElement(doc);
    if (node == nullptr || XMLCMP(node, NRP_NODE))
        throw std::invalid_argument("Invalid NRP XML.");
    node = xmlFirstElementChild(node);
    int monomer = 0;
    while (node != nullptr) {
        if (!XMLCMP(node, MONOMERS_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, MONOMER_NODE)) {
                    xmlNodePtr mnode = xmlFirstElementChild(nnode);
                    m_nrp.emplace_back();
                    auto &monomers = m_nrp.back();
                    while (mnode != nullptr) {
                        if (!XMLCMP(mnode, ID_NODE)) {
                            monomers.push_back(dbConn->getMonomer(std::atoi(XMLTXT(mnode))));
                        }
                        mnode = xmlNextElementSibling(mnode);
                    }
                }
                nnode = xmlNextElementSibling(nnode);
            }
        } else if (!XMLCMP(node, NRPS_NODE)) {
            xmlNodePtr nnode = xmlFirstElementChild(node);
            while (nnode != nullptr) {
                if (!XMLCMP(nnode, DOMAIN_NODE)) {
                    xmlNodePtr dnode = xmlFirstElementChild(nnode);
                    while (dnode != nullptr) {
                        if (!XMLCMP(dnode, ID_NODE)) {
                            m_nrps.push_back(dbConn->createDomain(std::atoi(XMLTXT(dnode))));
                            break;
                        }
                        dnode = xmlNextElementSibling(dnode);
                    }
                    dbConn->fillDomain(m_nrps.back());
                }
                nnode = xmlNextElementSibling(nnode);
            }
        }
        node = xmlNextElementSibling(node);
    }
}
