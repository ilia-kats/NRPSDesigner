#ifndef NRPSDESIGNER_NRP_H
#define NRPSDESIGNER_NRP_H

#include "nrpsdesigner_export.h"
#include "monomer.h"

#include <vector>

#include <libxml/xmlstring.h>

namespace nrps
{
class NRPSDESIGNER_EXPORT Nrp : public std::vector<Monomer>
{
public:
    enum class Type {Linear, Circular};
    
    Nrp(Type type = Type::Linear);
    Nrp(const char*);
    Nrp(const std::string&);
    Nrp(int);
    ~Nrp();
    
    Type type() const;
    
    void fromFile(const char*);
    void fromFile(const std::string&);
    void fromFile(int);
    
private:
    Type m_type;
    
    static xmlChar *s_nrp_node;
    static xmlChar *s_type_attr;
    static xmlChar *s_monomer_node;
    static xmlChar *s_name_node;
    static xmlChar *s_id_node;
    static xmlChar *s_configuration_node;
    static xmlChar *s_modification_node;
    
    static xmlChar *s_type_linear;
    static xmlChar *s_type_circular;
    static xmlChar *s_modification_nmethyl;
};
}

#endif
