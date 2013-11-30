#ifndef NRPSDESIGNER_GLOBALS_INTERNAL
#define NRPSDESIGNER_GLOBALS_INTERNAL

#include "config.h"
#include "domain.h"
#include "taxon.h"

#include <memory>
#include <unordered_map>
#include <fstream>

#ifdef WITH_INTERNAL_XML
#include <libxml/xmlwriter.h>
#include <unistd.h>
#endif

#define XMLCMP(n, y) xmlStrcmp(n->name, BAD_CAST y)
#define XMLTXT(n) (const char*)n->children->content


namespace std
{
template<>
struct NRPSDESIGNER_EXPORT hash<pair<uint32_t, uint32_t>>
{
public:
    size_t operator()(const pair<uint32_t, uint32_t>&) const;
};
}


namespace nrps
{
    NRPSDESIGNER_EXPORT float makeWeight(const std::shared_ptr<Domain>&, const std::shared_ptr<Domain>&, std::unordered_map<uint32_t, std::shared_ptr<Taxon>> *taxonCache = nullptr, std::unordered_map<std::pair<uint32_t, uint32_t>, float> *weightCache = nullptr);
#ifdef WITH_INTERNAL_XML
    void startXml(xmlTextWriterPtr);
    void endXml(xmlTextWriterPtr);

    template<class Func, typename... Args>
    std::string toXml(Func &&func, Args&&... args)
    {
        xmlBufferPtr buf = xmlBufferCreate();
        xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
        startXml(writer);
        std::bind(func, args..., writer)();
        endXml(writer);
        std::string xml((const char*)buf->content);
        xmlFreeTextWriter(writer);
        xmlBufferFree(buf);
        return xml;
    }
    template<class Func, typename... Args>
    void toXml(std::ostream &of, Func &&func, Args&&... args)
    {
        xmlBufferPtr buf = xmlBufferCreate();
        xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
        startXml(writer);
        std::bind(func, args..., writer)();
        endXml(writer);
        of << (const char*)buf->content;
        xmlFreeTextWriter(writer);
        xmlBufferFree(buf);
    }
    template<class Func, typename... Args>
    void toXml(const std::string &of, Func &&func, Args&&... args)
    {
        toXml(of.c_str, func, args...);
    }
    template<class Func, typename... Args>
    void toXml(const char *file, Func &&func, Args&&... args)
    {
        xmlTextWriterPtr writer = xmlNewTextWriterFilename(file, 0);
        startXml(writer);
        std::bind(func, args..., writer)();
        endXml(writer);
        xmlFreeTextWriter(writer);
    }
    template<class Func, typename... Args>
    void toXml(int fd, Func &&func, Args&&... args)
    {
        xmlBufferPtr buf = xmlBufferCreate();
        xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
        startXml(writer);
        std::bind(func, args..., writer)();
        endXml(writer);
        write(fd, buf->content, buf->use);
        xmlFreeTextWriter(writer);
        xmlBufferFree(buf);
    }
#endif
}

#endif
