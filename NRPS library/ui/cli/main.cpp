#include <nrp.h>
#include <generalpathway.h>

#include <iostream>

#include <curl/curl.h>
#include <string>

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
        if (monomer.modifications() & static_cast<Monomer::modification_type>(Monomer::Modification::Nmethyl))
            std::cout << "\tN-methylation" << std::endl;
    }

    GeneralPathway pathway(nrp);
    for (const std::shared_ptr<AbstractDomainType> &ptr : pathway) {
        std::cout << static_cast<int>(ptr->type()) << std::endl;
    }

    char url[] = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=562&rettype=xml";
    std::string xml;
    CURLcode init = curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL *handle = curl_easy_init();
    CURLcode ret = curl_easy_setopt(handle, CURLOPT_URL, url);
    ret = curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, write_fun);
    ret = curl_easy_setopt(handle, CURLOPT_WRITEDATA, &xml);

    ret = curl_easy_perform(handle);
    curl_easy_cleanup(handle);
    curl_global_cleanup();
    std::cout << xml << std::endl;
}
