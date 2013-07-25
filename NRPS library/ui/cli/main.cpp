#include <nrp.h>
#include <generalpathway.h>
#include <abstractdatabaseconnector.h>
#include <taxon.h>

#include <iostream>

#include <curl/curl.h>
#include <string>

#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>

using namespace nrps;

size_t write_fun (char *ptr, size_t size, size_t nmemb, void *userdata)
{
    ((std::string*)userdata)->append(ptr, size * nmemb);
    return size * nmemb;
}

int main(int argc, char *argv[])
{
    Nrp nrp(argv[1]);

    std::cout << "Type: " << (nrp.type() == Nrp::Type::Linear ? "linear" : "circular") << std::endl;
    for (const Monomer &monomer : nrp) {
        std::cout << std::endl
        << "ID: " << monomer.id() << std::endl
        << "name: " << monomer.name() << std::endl
        << "configuration:" << (monomer.configuration() == Configuration::L ? "L" : "D") << std::endl
        << "modification:" << std::endl;
        if (monomer.hasModification(Monomer::Modification::Nmethyl))
            std::cout << "\tN-methylation" << std::endl;
    }

    GeneralPathway pathway(nrp);
    for (const std::shared_ptr<AbstractDomainType> &ptr : pathway) {
        std::cout << static_cast<int>(ptr->type()) << std::endl;
    }
    auto dbConn = AbstractDatabaseConnector::getInstance();
    dbConn->initialize();
    auto possiblePathway = dbConn->getPotentialDomains(pathway);
    for (const auto &position : possiblePathway) {
        if (!position->empty()) {
            for (const auto &domain : *position) {
                std::cout << domain->domainId() << "\t";
            }
        }
        std::cout << std::endl;
    }

    Taxon taxon(562);
    char time[100];
    std::time_t t = std::chrono::system_clock::to_time_t(taxon.pubDate());
    std::strftime(time, 100, "%F %T", std::localtime(&t));
    std::cout << "scientific name: " << taxon.scientificName() << std::endl
    << "rank: " << static_cast<int>(taxon.rank()) << std::endl
    << "pubDate: " << time << std::endl << "lineage: ";
    for (const auto &t : taxon.lineage()) {
        std::cout << std::endl << "\t" << "id: " << t->id() << " name: " << t->scientificName();
    }
    std::cout << std::endl;

    /*char url[] = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=562&rettype=xml";
    std::string xml;
    CURLcode init = curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL *handle = curl_easy_init();
    CURLcode ret = curl_easy_setopt(handle, CURLOPT_URL, url);
    ret = curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, write_fun);
    ret = curl_easy_setopt(handle, CURLOPT_WRITEDATA, &xml);

    ret = curl_easy_perform(handle);
    curl_easy_cleanup(handle);
    curl_global_cleanup();
    std::cout << xml << std::endl;*/

    try {
        sql::Driver *driver = get_driver_instance();
        sql::Connection *con = driver->connect("tcp://127.0.0.1:3306", "root", "");
        con->setSchema("nrps_designer");
        sql::Statement *stmt = con->createStatement();
        sql::ResultSet *res = stmt->executeQuery("SELECT * FROM `main`");
        sql::ResultSetMetaData *resMeta = res->getMetaData();
        for (unsigned int i = 1; i <= resMeta->getColumnCount(); ++i)
            std::cout << resMeta->getColumnName(i) << "\t";
        std::cout << std::endl;
        while (res->next()) {
            for (unsigned int i = 1; i <= resMeta->getColumnCount(); ++i)
                std::cout << res->getString(i) << "\t";
            std::cout << std::endl;
        }
        con->close();
        delete res;
        delete stmt;
        delete con;
    } catch (sql::SQLException &e) {
        std::cout << "# ERR: SQLException in " << __FILE__;
        std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << std::endl;
        std::cout << "# ERR: " << e.what();
        std::cout << " (MySQL error code: " << e.getErrorCode();
        std::cout << ", SQLState: " << e.getSQLState() << " )" << std::endl;
    }
}
