# - Try to find MySQL / MySQL Embedded library
# Find the MySQL includes and client library
# This module defines
#  MYSQL_INCLUDE_DIR, where to find mysql.h
#  MYSQL_LIBRARIES, the libraries needed to use MySQL.
#  MYSQL_LIB_DIR, path to the MYSQL_LIBRARIES
#  MYSQL_EMBEDDED_LIBRARIES, the libraries needed to use MySQL Embedded.
#  MYSQL_EMBEDDED_LIB_DIR, path to the MYSQL_EMBEDDED_LIBRARIES
#  MYSQL_FOUND, If false, do not try to use MySQL.
#  MYSQL_EMBEDDED_FOUND, If false, do not try to use MySQL Embedded.

# Copyright (c) 2013 Ilia Kats, based on FindMYSQL.cmake, (c) 2006-2008, Jarosław Staniek <staniek@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying kde/COPYING-CMAKE-SCRIPTS file.

if(WIN32)
    find_path(MYSQLCPP_INCLUDE_DIR mysql_connection.h
        PATHS
        $ENV{MYSQL_INCLUDE_DIR}
        $ENV{MYSQL_DIR}/include
        $ENV{ProgramFiles}/MySQL/*/include
        $ENV{SystemDrive}/MySQL/*/include
   )
else(WIN32)
    find_path(MYSQLCPP_INCLUDE_DIR mysql_connection.h
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
        ${MYSQL_LIB_PATHS}
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
        ${MYSQL_LIB_PATHS}
   )
endif(WIN32)

if(MYSQLCPP_LIBRARIES)
    get_filename_component(MYSQLCPP_LIB_DIR ${MYSQLCPP_LIBRARIES} PATH)
endif(MYSQLCPP_LIBRARIES)

if(MYSQLCPP_EMBEDDED_LIBRARIES)
    get_filename_component(MYSQLCPP_EMBEDDED_LIB_DIR ${MYSQLCPP_EMBEDDED_LIBRARIES} PATH)
endif(MYSQLCPP_EMBEDDED_LIBRARIES)

if(MYSQLCPP_INCLUDE_DIR AND MYSQLCPP_LIBRARIES)
    set(MYSQLCPP_FOUND TRUE)
    message(STATUS "Found MySQL C++ connector: ${MYSQLCPP_INCLUDE_DIR}, ${MYSQLCPP_LIBRARIES}")
else(MYSQLCPP_INCLUDE_DIR AND MYSQLCPP_LIBRARIES)
    set(MYSQLCPP_FOUND FALSE)
    message(STATUS "MySQL C++ connector not found.")
endif(MYSQLCPP_INCLUDE_DIR AND MYSQLCPP_LIBRARIES)

if(MYSQLCPP_INCLUDE_DIR AND MYSQLCPP_EMBEDDED_LIBRARIES AND HAVE_MYSQLCPP_OPT_EMBEDDED_CONNECTION)
    set(MYSQLCPP_EMBEDDED_FOUND TRUE)
    message(STATUS "Found MySQL Embedded: ${MYSQLCPP_INCLUDE_DIR}, ${MYSQLCPP_EMBEDDED_LIBRARIES}")
else(MYSQLCPP_INCLUDE_DIR AND MYSQLCPP_EMBEDDED_LIBRARIES AND HAVE_MYSQLCPP_OPT_EMBEDDED_CONNECTION)
    set(MYSQLCPP_EMBEDDED_FOUND FALSE)
    message(STATUS "MySQL Embedded not found.")
endif(MYSQLCPP_INCLUDE_DIR AND MYSQLCPP_EMBEDDED_LIBRARIES AND HAVE_MYSQLCPP_OPT_EMBEDDED_CONNECTION)

mark_as_advanced(MYSQLCPP_INCLUDE_DIR MYSQLCPP_LIBRARIES MYSQLCPP_EMBEDDED_LIBRARIES)
