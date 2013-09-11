#include "networkoptions.h"

using namespace nrps;
namespace po = boost::program_options;

NetworkOptions *NetworkOptions::s_instance = nullptr;

NetworkOptions::NetworkOptions()
: m_timeout(0)
{}

NetworkOptions::~NetworkOptions()
{
    s_instance = nullptr;
}

NetworkOptions* NetworkOptions::getInstance()
{
    if (s_instance == nullptr) {
        s_instance = new NetworkOptions();
    }
    return s_instance;
}

po::options_description NetworkOptions::options()
{
    po::options_description options("Networking options");
    options.add_options()("timeout", po::value<long>(&m_timeout)->default_value(30), "Timeout value for connection attempts.");
    return options;
}

long NetworkOptions::timeout() const
{
    return m_timeout;
}

void NetworkOptions::setTimeout(long timeout)
{
    m_timeout = timeout;
}
