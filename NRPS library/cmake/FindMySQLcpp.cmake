# - Try to find MySQL C++
# Find the MySQL includes and client library
# This module defines
#  MYSQLCPP_INCLUDE_DIR, where to find mysql.h
#  MYSQLCPP_LIBRARIES, the libraries needed to use MySQL.
#  MYSQLCPP_LIB_DIR, path to the MYSQL_LIBRARIES
#  MYSQLCPP_FOUND, If false, do not try to use MySQL.

# Copyright (c) 2013 Ilia Kats, based on FindMYSQL.cmake, (c) 2006-2008, Jarosław Staniek <staniek@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying kde/COPYING-CMAKE-SCRIPTS file.

if(WIN32)
    find_path(MYSQLCPP_INCLUDE_DIR cppconn/driver.h
        PATHS
        $ENV{MYSQL_INCLUDE_DIR}
        $ENV{MYSQL_DIR}/include
        $ENV{ProgramFiles}/MySQL/*/include
        $ENV{SystemDrive}/MySQL/*/include
   )
else(WIN32)
    find_path(MYSQLCPP_INCLUDE_DIR cppconn/driver.h
        PATHS
        $ENV{MYSQL_INCLUDE_DIR}
        $ENV{MYSQL_DIR}/include
        /usr/local/mysql/include
        /opt/mysql/mysql/include
        PATH_SUFFIXES
        mysql
   )
endif(WIN32)

if(WIN32)
    set(MYSQLCPP_LIB_PATHS
        $ENV{MYSQL_DIR}/lib/opt
        $ENV{MYSQL_DIR}/client/release
        $ENV{ProgramFiles}/MySQL/*/lib/opt
        $ENV{SystemDrive}/MySQL/*/lib/opt
   )
   find_library(MYSQLCPP_LIBRARIES NAMES mysqlcppconn
        PATHS
        ${MYSQLCPP_LIB_PATHS}
   )
else(WIN32)
    set(MYSQLCPP_LIB_PATHS
        $ENV{MYSQL_DIR}/libmysql_r/.libs
        $ENV{MYSQL_DIR}/lib
        $ENV{MYSQL_DIR}/lib/mysql
        /usr/local/mysql/lib
        /opt/mysql/mysql/lib
        PATH_SUFFIXES
        mysql
   )
   find_library(MYSQLCPP_LIBRARIES NAMES mysqlcppconn
        PATHS
        ${MYSQLCPP_LIB_PATHS}
   )
endif(WIN32)

if(MYSQLCPP_LIBRARIES)
    get_filename_component(MYSQLCPP_LIB_DIR ${MYSQLCPP_LIBRARIES} PATH)
endif(MYSQLCPP_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MySQLCpp  DEFAULT_MSG  MYSQLCPP_LIBRARIES MYSQLCPP_INCLUDE_DIR)

mark_as_advanced(MYSQLCPP_INCLUDE_DIR MYSQLCPP_LIBRARIES )
