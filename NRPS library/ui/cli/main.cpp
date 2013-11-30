#include <monomer.h>
#include <abstractdatabaseconnector.h>
#include <nrps.h>
#include <exceptions.h>
#include <nrpsbuilder.h>
#include <taxonbuilder.h>
#include <nrpslibrary.h>
#include <networkoptions.h>

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
    std::string outsbol;
    bool indTag = false;
    std::string inputnrps;
    po::options_description options("General options");
    options.add_options()("help,h", "print this help")
                         ("monomers,m", po::value<std::string>(&monomersopt), "Comma-separated list of monomer ids for the NRP.")
#ifdef WITH_INTERNAL_XML
                         ("outfile,o", po::value<std::string>(&outfile), "Path to output file in internal XML format. Use - for stdout.")
#endif
#ifdef WITH_SBOL
                         ("outsbol,s", po::value<std::string>(&outsbol), "Path to output file in SBOL format. Use - for stdout.")
#endif
                         ("indigoidine-tag,t", po::bool_switch(&indTag), "Append an indigoidine tag to the NRP.")
                         ("input-nrps,n", po::value<std::string>(&inputnrps), "Path to input NRPS file for library generation. Use - for stdin.");
    auto dbConn = AbstractDatabaseConnector::getInstance();

    TaxonBuilder *tb = TaxonBuilder::getInstance();
    po::options_description alloptions;
    alloptions.add(options).add(dbConn->options()).add(NetworkOptions::getInstance()->options()).add(tb->options());

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(alloptions).run(), vm);
    po::store(po::parse_environment(tb->options(), [&tb](std::string s){return tb->mapEnvironment(s);}), vm);

    if (argc == 1 || vm.count("help") || !vm.count("monomers") && !vm.count("input-nrps")) {
        std::cout << alloptions;
        delete dbConn;
        return 0;
    }
    po::notify(vm);

    try {
        dbConn->initialize();

        if (monomersopt.size()) {
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

            Nrps nrps = NrpsBuilder().build(monomers, indTag);
#ifdef WITH_INTERNAL_XML
            if (!outfile.empty()) {
                if (outfile == "-")
                    nrps.toXml(std::cout);
                else
                    nrps.toXml(outfile);
            }
#endif
#ifdef WITH_SBOL
            if (!outsbol.empty()) {
                if (outsbol == "-")
                    nrps.toSbol(std::cout);
                else
                    nrps.toSbol(outsbol);
            }
#endif
        } else {
            NrpsLibrary lib;
            if (inputnrps == "-")
                lib.fromFile(std::cin);
            else
                lib.fromFile(inputnrps);
            lib.makeLibrary();
#ifdef WITH_INTERNAL_XML
            if (!outfile.empty()) {
                if (outfile == "-")
                    lib.toXml(std::cout);
                else
                    lib.toXml(outfile);
            }
#endif
#ifdef WITH_SBOL
            if (!outsbol.empty()) {
                if (outsbol == "-")
                    lib.toSbol(std::cout);
                else
                    lib.toSbol(outsbol);
            }
#endif
        }
        delete dbConn;
        return 0;
    } catch (const NCBITaxonomyError &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    } catch(const TaxonomyDumpError &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 2;
    } catch (const NetworkError &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 3;
    } catch (const DatabaseError &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 4;
    } catch (const std::logic_error &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 5;
    } catch (const std::system_error &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 6;
    } catch (const std::exception &e) {
        delete dbConn;
        std::cerr << "ERROR: " << e.what() << std::endl;
        return -1;
    }
}
