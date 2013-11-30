#ifndef NRPSDESIGNER_PRODUCT_H
#define NRPSDESIGNER_PRODUCT_H

#include "config.h"
#include "nrpsdesigner_export.h"

#include <cstdint>
#include <string>
#include <unordered_map>

#ifdef WITH_INTERNAL_XML
#include <libxml/xmlwriter.h>
#endif

namespace nrps
{
class NRPSDESIGNER_EXPORT Product
{
public:
    ~Product();
   uint32_t id() const;
   const std::string& name() const;
   const std::string& description() const;

   void setId(uint32_t);
   void setName(const std::string&);
   void setName(std::string&&);
   void setDescription(const std::string&);
   void setDescription(std::string&&);

#ifdef WITH_INTERNAL_XML
   void toXml(xmlTextWriterPtr) const;
#endif

   static Product* makeProduct(uint32_t);

private:
    Product(uint32_t);

    static std::unordered_map<uint32_t, Product*> s_products;

    uint32_t m_id;
    std::string m_name;
    std::string m_description;
};
}

#endif
