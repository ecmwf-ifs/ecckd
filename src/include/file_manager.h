/// @file      file_manager.h
/// @brief     Functions to manage a list of directories to search for input files
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright 2015-2016 European Centre for Medium Range Weather Forecasts
/// @license   Apache License Version 2 or ESA Software Community License Type 1
///            (see the NOTICE.md file for details)

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
