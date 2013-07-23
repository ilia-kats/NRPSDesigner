#include "abstractdatabaseconnector.h"
#include "mysqldatabaseconnector.h"

using namespace nrps;

AbstractDatabaseConnector *AbstractDatabaseConnector::s_instance = nullptr;

AbstractDatabaseConnector* AbstractDatabaseConnector::getInstance()
{
    if (s_instance == nullptr) {
        s_instance = new MySQLDatabaseConnector();
    }
    return s_instance;
}
