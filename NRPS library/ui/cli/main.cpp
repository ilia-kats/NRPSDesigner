#include <monomer.h>
#include <abstractdatabaseconnector.h>
#include <taxon.h>
#include <nrps.h>
#include <nrpsbuilder.h>
#include <origin.h>

#include <iostream>

#include <curl/curl.h>
#include <string>
#include <cstring>
#include <cstdlib>

#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>

using namespace nrps;

size_t write_fun (char *ptr, size_t size, size_t nmemb, void *userdata)
{
    ((std::string*)userdata)->append(ptr, size * nmemb);
    return size * nmemb;
}

int main(int argc, char *argv[])
{
    auto dbConn = AbstractDatabaseConnector::getInstance();
    dbConn->initialize();

    std::vector<Monomer> monomers;
    char *ptr = std::strtok(argv[1], ",");
    while (ptr != nullptr) {
        monomers.push_back(dbConn->getMonomer(std::atoi(ptr)));
        ptr = std::strtok(nullptr, ",");
    }

    for (const Monomer &m : monomers) {
        std::cout << m.name() << std::endl
        << "\tid: " << m.id() << std::endl
        << "\tparentId: " << m.parentId() << std::endl
        << "\tenantiomerId:" << m.enantiomerId() << std::endl
        << "\tconfiguration:" << (int)m.configuration() << std::endl
        << "\tmodifications:" << std::endl;
        for (const uint32_t &mod : m.modifications()) {
            std::cout << "\t\t" << mod << std::endl;
        }
    }

    try {
        uint32_t id = 8;
        Origin *ori = Origin::makeOrigin(id);
        dbConn->fillOrigin(ori);
        while (ori != nullptr) {
            std::cout << "Origin: id: " << ori->id() << std::endl
            << "\tsourceType: " << toString(ori->sourceType()) << std::endl
            << "\tsource: " << ori->source() << std::endl
            << "\tspecies: " << ori->species() << std::endl
            << "\tdesc: " << ori->description() << std::endl
            << "\ttaxid: " << ori->taxId() << std::endl;
            ori = ori->parent();
        }

    } catch (sql::SQLException &e) {
        std::cout << "# ERR: SQLException in " << __FILE__;
        std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << std::endl;
        std::cout << "# ERR: " << e.what();
        std::cout << " (MySQL error code: " << e.getErrorCode();
        std::cout << ", SQLState: " << e.getSQLState() << " )" << std::endl;
    }
    Nrps nrps = NrpsBuilder().build(monomers);
    nrps.toXml(std::cout);

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

    Taxon tac(54914);
    auto d = taxon - tac;
    std::cout << "diff1: " << (int)d[0] << " ; diff2: " << (int)d[1] << std::endl;

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

//     try {
//         sql::Driver *driver = get_driver_instance();
//         sql::Connection *con = driver->connect("tcp://127.0.0.1:3306", "root", "");
//         con->setSchema("nrps_designer");
//         //sql::PreparedStatement *stmt = con->prepareStatement("SELECT @id:=16, @enant_id:=17, domain.id AS did, substr_xref.substrate_id AS sid, origin.id AS orid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_domain_substrateSpecificity AS substr_xref, databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype) ON (domain.id = substr_xref.id AND domain.cds_id = cds.id AND cds.origin_id = origin.id AND domain.domainType_id = dtype.id) WHERE (0 OR substr_xref.substrate_id = @id OR substr_xref.substrate_id = @enant_id) AND dtype.name = 'A' AND domain.chirality = 'N' ORDER BY IF(substr_xref.substrate_id = @id, 0, 1) ASC;");
//         //stmt->setUInt(1, 16);
//         //stmt->setUInt(2, 17);
//         //stmt->setString(2, std::string("A"));
//         sql::Statement *stmt = con->createStatement();
//         /*sql::ResultSet *res = */stmt->execute("SET @id:=16; SET @enant_id:=17; SELECT domain.id AS did, substr_xref.substrate_id AS sid, origin.id AS orid FROM databaseInput_domain AS domain INNER JOIN (databaseInput_domain_substrateSpecificity AS substr_xref, databaseInput_cds AS cds, databaseInput_origin AS origin, databaseInput_type AS dtype) ON (domain.id = substr_xref.id AND domain.cds_id = cds.id AND cds.origin_id = origin.id AND domain.domainType_id = dtype.id) WHERE (0 OR substr_xref.substrate_id = @id OR substr_xref.substrate_id = @enant_id) AND dtype.name = 'A' AND domain.chirality = 'N' ORDER BY IF(substr_xref.substrate_id = @id, 0, 1) ASC;");
//         /*sql::ResultSetMetaData *resMeta = res->getMetaData();
//         for (unsigned int i = 1; i <= resMeta->getColumnCount(); ++i)
//             std::cout << resMeta->getColumnName(i) << "\t";
//         std::cout << std::endl;
//         while (res->next()) {
//             for (unsigned int i = 1; i <= resMeta->getColumnCount(); ++i)
//                 std::cout << res->getString(i) << "\t";
//             std::cout << std::endl;
//         }*/
//         con->close();
//         //delete res;
//         delete stmt;
//         delete con;
//     } catch (sql::SQLException &e) {
//         std::cout << "# ERR: SQLException in " << __FILE__;
//         std::cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << std::endl;
//         std::cout << "# ERR: " << e.what();
//         std::cout << " (MySQL error code: " << e.getErrorCode();
//         std::cout << ", SQLState: " << e.getSQLState() << " )" << std::endl;
//     }
}
