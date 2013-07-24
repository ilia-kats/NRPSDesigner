#include "nrp.h"
#include "monomer.h"

#include <ios>
#include <stdexcept>

#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>

#include <libxml/xmlreader.h>

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

using namespace nrps;

Nrp::Nrp(Nrp::Type type)
: std::vector<Monomer>(), m_type(type)
{}

Nrp::Nrp(const char *file)
: std::vector<Monomer>(), m_type(Nrp::Type::Linear)
{
    fromFile(file);
}

Nrp::Nrp(const std::string &file)
: std::vector<Monomer>(), m_type(Nrp::Type::Linear)
{
    fromFile(file);
}

Nrp::Nrp(int fd)
: std::vector<Monomer>(), m_type(Nrp::Type::Linear)
{
    fromFile(fd);
}

Nrp::~Nrp()
{}

Nrp::Type Nrp::type() const
{
    return m_type;
}

void Nrp::fromFile(const char *file)
{
    int fd = open(file, O_RDONLY);
    fromFile(fd);
    close(fd);
}

void Nrp::fromFile(const std::string &file)
{
    fromFile(file.c_str());
}

void Nrp::fromFile(int fd)
{
    xmlTextReaderPtr reader = xmlReaderForFd(fd, nullptr, nullptr, XML_PARSE_RECOVER | XML_PARSE_NONET | XML_PARSE_COMPACT);
    if (reader == nullptr)
        throw std::ios_base::failure("Could not initialize XML reader");
    int ret;
    bool nrp_found = false;
    xmlChar *name;
    xmlChar *value;
    while (!nrp_found && (ret = xmlTextReaderRead(reader)) == 1) {
        name = xmlTextReaderName(reader);
        nrp_found = !xmlStrcmp(name, s_nrp_node);
        xmlFree(name);
    }
    if (!nrp_found)
        throw std::invalid_argument("Could not find nrp root node");
    clear();
    value = xmlTextReaderGetAttribute(reader, s_type_attr);
    if (!xmlStrcmp(value, s_type_circular))
        m_type = Type::Circular;
    xmlFree(value);
    ret = skip_whitespace(reader);
    name = xmlTextReaderName(reader);
    bool in_monomer = !xmlStrcmp(xmlTextReaderName(reader), s_monomer_node);
    while (ret == 1 && in_monomer) {
        Monomer monomer;
        in_monomer = false;
        ret = skip_whitespace(reader);
        name = xmlTextReaderName(reader);
        while (ret == 1 && !in_monomer) {
            if (!xmlStrcmp(name, s_id_node)) {
                value = xmlTextReaderReadInnerXml(reader);
                monomer.setId(std::atoi((const char*)value));
                xmlFree(value);
            } else if (monomer.name().empty() && !xmlStrcmp(name, s_name_node)) {
                value = xmlTextReaderReadInnerXml(reader);
                monomer.setName((const char*)value);
                xmlFree(value);
            } else if (!xmlStrcmp(name, s_configuration_node)) {
                value = xmlTextReaderReadInnerXml(reader);
                monomer.setConfiguration(*((const char*)value) == 'D' ? Configuration::D : Configuration::L);
                xmlFree(value);
            } else if (!xmlStrcmp(name, s_modification_node)) {
                value = xmlTextReaderReadInnerXml(reader);
                if (!xmlStrcmp(value, s_modification_nmethyl))
                    monomer.addModification(Monomer::Modification::Nmethyl);
                xmlFree(value);
            }
            xmlFree(name);
            ret = skip_whitespace(reader);
            name = xmlTextReaderName(reader);
            in_monomer = !xmlStrcmp(name, s_monomer_node) && !(xmlTextReaderNodeType(reader) == XML_READER_TYPE_END_ELEMENT);
        }
        xmlFree(name);
        push_back(std::move(monomer));
    }
    xmlFreeTextReader(reader);
}

#define NRP_NODE "nrp"
#define TYPE_ATTR "type"
#define MONOMER_NODE "monomer"
#define NAME_NODE "name"
#define ID_NODE "id"
#define CONFIGURATION_NODE "configuration"
#define MODIFICATION_NODE "modification"

xmlChar *Nrp::s_nrp_node           = xmlCharStrdup(NRP_NODE);
xmlChar *Nrp::s_type_attr          = xmlCharStrdup(TYPE_ATTR);
xmlChar *Nrp::s_monomer_node       = xmlCharStrdup(MONOMER_NODE);
xmlChar *Nrp::s_name_node          = xmlCharStrdup(NAME_NODE);
xmlChar *Nrp::s_id_node            = xmlCharStrdup(ID_NODE);
xmlChar *Nrp::s_configuration_node = xmlCharStrdup(CONFIGURATION_NODE);
xmlChar *Nrp::s_modification_node  = xmlCharStrdup(MODIFICATION_NODE);

#define TYPE_LINEAR "linear"
#define TYPE_CIRCULAR "circular"
#define MODIFICATION_NMETHYL "N-methylation"
xmlChar *Nrp::s_type_linear          = xmlCharStrdup(TYPE_LINEAR);
xmlChar *Nrp::s_type_circular        = xmlCharStrdup(TYPE_CIRCULAR);
xmlChar *Nrp::s_modification_nmethyl = xmlCharStrdup(MODIFICATION_NMETHYL);
