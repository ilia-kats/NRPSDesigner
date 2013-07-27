#ifndef NRPSDESIGNER_NRPS_H

#include "nrpsdesigner_export.h"
#include "domain.h"

#include <vector>

namespace nrps
{
class NRPSDESIGNER_EXPORT Nrps : public std::vector<std::shared_ptr<Domain>>
{
public:
    Nrps(const std::vector<std::shared_ptr<std::vector<std::shared_ptr<Domain>>>>&);
};
}

#endif
