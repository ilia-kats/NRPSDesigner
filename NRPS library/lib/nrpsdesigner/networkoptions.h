#ifndef NRPSDESIGNER_NETWORKOPTIONS_H
#define NRPSDESIGNER_NETWORKOPTIONS_H

#include "nrpsdesigner_export.h"

#include <boost/program_options/options_description.hpp>

namespace nrps
{
class NetworkOptions
{
public:
    ~NetworkOptions();

    boost::program_options::options_description options();

    long timeout() const;
    void setTimeout(long);

    static NetworkOptions* getInstance();

private:
    NetworkOptions();

    long m_timeout;

    static NetworkOptions* s_instance;
};
}

#endif
