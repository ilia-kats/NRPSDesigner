#ifndef NRPSDESIGNER_NRPSLIBRARY_H
#define NRPSDESIGNER_NRPSLIBRARY_H

#include "config.h"
#include "nrpsdesigner_export.h"
#include "monomer.h"
#include "nrps.h"

#include <vector>
#include <iostream>

#include <libxml/tree.h>
#ifdef WITH_INTERNAL_XML
#include <libxml/xmlwriter.h>
#endif

namespace nrps
{
class NRPSDESIGNER_EXPORT NrpsLibrary
{
public:
    NrpsLibrary(const char*);
    NrpsLibrary(const std::string&);
    NrpsLibrary(int);
    NrpsLibrary(std::istream&);
    NrpsLibrary();
    ~NrpsLibrary();

    void fromFile(const char*);
    void fromFile(const std::string&);
    void fromFile(int);
    void fromFile(std::istream&);

#ifdef WITH_INTERNAL_XML
    std::string toXml() const;
    void toXml(std::ostream&) const;
    void toXml(const std::string&) const;
    void toXml(const char*) const;
    void toXml(int) const;
    void toXml(xmlTextWriterPtr) const;
#endif
#ifdef WITH_SBOL
    std::string toSbol() const;
    void toSbol(std::ostream&) const;
    void toSbol(const std::string&) const;
    void toSbol(const char*) const;
    void toSbol(int) const;
    Collection* toSbol(Document*) const;
#endif

    void makeLibrary();

private:
    void readXml(xmlDocPtr);

    std::vector<std::vector<Monomer>> m_nrp;
    Nrps m_nrps;
    std::vector<Nrps> m_library;
};
}

#endif
