#include "product.h"

#define PRODUCT_NODE "product"
#define ID_NODE "id"
#define NAME_NODE "name"
#define DESCRIPTION_NODE "description"

using namespace nrps;

std::unordered_map<uint32_t, Product*> Product::s_products = std::unordered_map<uint32_t, Product*>();

Product::Product(uint32_t id)
: m_id(id)
{}

Product::~Product()
{
    s_products.erase(m_id);
}

uint32_t Product::id() const
{
    return m_id;
}

void Product::setId(uint32_t id)
{
    m_id = id;
}

const std::string& Product::name() const
{
    return m_name;
}

void Product::setName(const std::string &name)
{
    m_name = name;
}

void Product::setName(std::string &&name)
{
    m_name = std::move(name);
}

const std::string& Product::description() const
{
    return m_description;
}

void Product::setDescription(const std::string &desc)
{
    m_description = desc;
}

void Product::setDescription(std::string &&desc)
{
    m_description = std::move(desc);
}

Product* Product::makeProduct(uint32_t id)
{
    Product *pw;
    if (s_products.count(id) == 0) {
        pw = new Product(id);
        s_products.emplace(id, pw);
    } else {
        pw = s_products[id];
    }
    return pw;
}

#ifdef WITH_INTERNAL_XML
void Product::toXml(xmlTextWriterPtr writer) const
{
    xmlTextWriterStartElement(writer, BAD_CAST PRODUCT_NODE);
    xmlTextWriterWriteElement(writer, BAD_CAST ID_NODE, BAD_CAST std::to_string(id()).c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST NAME_NODE, BAD_CAST name().c_str());
    xmlTextWriterWriteElement(writer, BAD_CAST DESCRIPTION_NODE, BAD_CAST description().c_str());
    xmlTextWriterEndElement(writer);
}
#endif
