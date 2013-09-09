#ifndef NRPSDESIGNER_GLOBAL_ENUMS_H
#define NRPSDESIGNER_GLOBAL_ENUMS_H

#include "nrpsdesigner_export.h"
#include <string>

namespace nrps
{
    enum class NRPSDESIGNER_EXPORT DomainType {A, AC, C, T, E, Te, Tailoring};
    enum class NRPSDESIGNER_EXPORT Configuration {None, L, D};
    enum class NRPSDESIGNER_EXPORT DomainTPosition {BeforeC, BeforeE};
    enum class NRPSDESIGNER_EXPORT OriginSourceType {Species, DNA};

    template<class T>
    std::string toString(T);
    template<>
    std::string NRPSDESIGNER_EXPORT toString(DomainType);
    template<>
    std::string NRPSDESIGNER_EXPORT toString(Configuration);
    template<>
    std::string NRPSDESIGNER_EXPORT toString(DomainTPosition);
    template<>
    std::string NRPSDESIGNER_EXPORT toString(OriginSourceType);

    template<class T>
    T fromString(const std::string&);
    template<>
    DomainType NRPSDESIGNER_EXPORT fromString(const std::string&);
    template<>
    Configuration NRPSDESIGNER_EXPORT fromString(const std::string&);
    template<>
    DomainTPosition NRPSDESIGNER_EXPORT fromString(const std::string&);
    template<>
    OriginSourceType NRPSDESIGNER_EXPORT fromString(const std::string&);
}

#endif
