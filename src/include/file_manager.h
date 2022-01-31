/// @file      file_manager.h
/// @brief     Functions to manage a list of directories to search for input files
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef file_manager_H
#define file_manager_H 1

#include <string>

#define OPTSYN_PATH_NAME "OPTSYN_PATH"
//#define AUX_DIRECTORY_NAME "aux/"
#define AUX_DIRECTORY_NAME ""

void append_search_directory(std::string dir);

void prepend_search_directory(std::string dir);

std::string find_file(std::string name);

std::string search_directories();

std::string first_directory();

#endif
