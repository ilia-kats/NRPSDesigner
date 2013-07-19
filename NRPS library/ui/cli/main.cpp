#include <nrp.h>

#include <iostream>

using namespace nrps;
int main(int argc, char *argv[])
{
    Nrp nrp(argv[1]);
    
    std::cout << "Type: " << (nrp.type() == Nrp::Type::Linear ? "linear" : "circular") << std::endl;
    for (const Monomer &monomer : nrp) {
        std::cout << std::endl
        << "ID: " << monomer.id() << std::endl
        << "name: " << monomer.name() << std::endl
        << "configuration:" << (monomer.configuration() == Monomer::Configuration::L ? "L" : "D") << std::endl
        << "modification:" << std::endl;
        if (monomer.modifications() & static_cast<Monomer::modification_type>(Monomer::Modification::Nmethyl))
            std::cout << "\tN-methylation" << std::endl;
    }
}
