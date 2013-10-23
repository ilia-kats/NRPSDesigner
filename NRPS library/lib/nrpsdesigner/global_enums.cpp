#include "global_enums.h"

using namespace nrps;

template<>
std::string nrps::toString(DomainType t)
{
    switch (t) {
        case DomainType::A:
            return "A";
            break;
        case DomainType::C:
            return "C";
            break;
        case DomainType::T:
            return "T";
            break;
        case DomainType::E:
            return "E";
            break;
        case DomainType::Te:
            return "TE";
            break;
        case DomainType::AOxA:
            return "AOxA";
            break;
        default:
            return std::string();
            break;
    }
}

template<>
DomainType NRPSDESIGNER_EXPORT nrps::fromString(const std::string &s)
{
    if (s == "A")
        return DomainType::A;
    else if (s == "C")
        return DomainType::C;
    else if (s == "T")
        return DomainType::T;
    else if (s == "E")
        return DomainType::E;
    else if (s == "TE")
        return DomainType::Te;
    else if (s == "AOxA")
        return DomainType::AOxA;
}

template<>
std::string nrps::toString(Configuration c)
{
    switch (c) {
        case Configuration::L:
            return "L";
            break;
        case Configuration::D:
            return "D";
            break;
        default:
            return "L";
            break;
    }
}

template<>
Configuration NRPSDESIGNER_EXPORT nrps::fromString(const std::string &s)
{
    if (s == "L")
        return Configuration::L;
    else if (s == "D")
        return Configuration::D;
}

template<>
std::string nrps::toString(DomainTPosition p)
{
    switch (p) {
        case DomainTPosition::BeforeC:
            return "C";
            break;
        case DomainTPosition::BeforeE:
            return "E";
            break;
        default:
            return std::string();
            break;
    }
}

template<>
DomainTPosition NRPSDESIGNER_EXPORT nrps::fromString(const std::string &s)
{
    if (s == "BeforeC")
        return DomainTPosition::BeforeC;
    else if (s == "BeforeE")
        return DomainTPosition::BeforeE;
}

template<>
std::string nrps::toString(OriginSourceType o)
{
    switch (o) {
        case OriginSourceType::Species:
            return "Species";
            break;
        case OriginSourceType::DNA:
            return "DNA";
            break;
        default:
            return std::string();
            break;
    }
}

template<>
OriginSourceType NRPSDESIGNER_EXPORT nrps::fromString(const std::string &s)
{
    if (s == "Species")
        return OriginSourceType::Species;
    else if (s == "DNA")
        return OriginSourceType::DNA;
}
