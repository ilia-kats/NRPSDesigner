#include <monomer.h>
#include <abstractdatabaseconnector.h>
#include <nrps.h>
#include <nrpsbuilder.h>

#include <curl/curl.h>

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

#include <boost/program_options.hpp>

using namespace nrps;
namespace po = boost::program_options;

size_t write_fun (char *ptr, size_t size, size_t nmemb, void *userdata)
{
    ((std::string*)userdata)->append(ptr, size * nmemb);
    return size * nmemb;
}

int main(int argc, char *argv[])
{
    std::string monomersopt;
    std::string outfile;
    po::options_description options("General options");
    options.add_options()("help,h", "print this help")
                         ("monomers,m", po::value<std::string>(&monomersopt)->required(), "Comma-separated list of monomer ids for the NRP.")
                         ("outfile,o", po::value<std::string>(&outfile)->default_value("-"), "Output file. Use - for stdout.");
    auto dbConn = AbstractDatabaseConnector::getInstance();

    po::options_description alloptions;
    alloptions.add(options).add(dbConn->options());

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(alloptions).run(), vm);

    if (vm.count("help")) {
        std::cout << alloptions;
        delete dbConn;
        return 0;
    }
    po::notify(vm);

    dbConn->initialize();

    std::vector<Monomer> monomers;
    char *monomersch = (char*)std::malloc((monomersopt.length() + 1) * sizeof(char));
    std::strncpy(monomersch, monomersopt.c_str(), monomersopt.length());
    monomersch[monomersopt.length()] = '\0';
    char *ptr = std::strtok(monomersch, ",");
    while (ptr != nullptr) {
        monomers.push_back(dbConn->getMonomer(std::atoi(ptr)));
        ptr = std::strtok(nullptr, ",");
    }
    std::free(monomersch);

    Nrps nrps = NrpsBuilder().build(monomers);

    if (outfile == "-")
        nrps.toXml(std::cout);
    else
        nrps.toXml(outfile);
    
    delete dbConn;
    return 0;
}
