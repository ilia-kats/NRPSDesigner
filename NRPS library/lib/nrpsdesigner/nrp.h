#ifndef NRPSDESIGNER_NRP_H
#define NRPSDESIGNER_NRP_H

#include "nrpsdesigner_export.h"
#include "monomer.h"

#include <vector>

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
};
}

#endif
