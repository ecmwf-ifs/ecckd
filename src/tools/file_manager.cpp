/// @file      file_manager.cpp
/// @brief     Implements functions to manage a list of directories to search for input files
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#include <cstdlib>
#include <cstdio>
#include "file_manager.h"
#include <iostream>
#include <deque>
typedef std::deque<std::string> directory_list;

static
directory_list
initialize_directory_list()
{
  directory_list list;
  char* path_cstr = getenv(OPTSYN_PATH_NAME);
  if (path_cstr) {
    std::string path = path_cstr;
    while (path.size() > 0) {
      size_t start = 0;
      size_t end = path.find_first_of(':');
      if (end == std::string::npos) {
	// No colon: just one directory
	list.push_back(path);
	break;
      }
      else {
	if (end > start) {
	  list.push_back(path.substr(start, end));
	}
	path = path.substr(end+1, std::string::npos);
      }
    }
  }
  else {
    list.push_front(".");
  }
  return list;
}


directory_list directory_list_ = initialize_directory_list();


void
append_search_directory(std::string path)
{
  //  directory_list_.push_back(dir);
  while (path.size() > 0) {
    size_t start = 0;
    size_t end = path.find_first_of(':');
    if (end == std::string::npos) {
      // No colon: just one directory
      directory_list_.push_back(path);
      break;
    }
    else {
      if (end > start) {
	directory_list_.push_back(path.substr(start, end));
      }
      path = path.substr(end+1, std::string::npos);
    }
  }
}

void
prepend_search_directory(std::string path)
{
  directory_list tmp_list;
  //  directory_list_.push_front(dir);
  while (path.size() > 0) {
    size_t start = 0;
    size_t end = path.find_first_of(':');
    if (end == std::string::npos) {
      // No colon: just one directory
      tmp_list.push_back(path);
      break;
    }
    else {
      if (end > start) {
	tmp_list.push_back(path.substr(start, end));
      }
      path = path.substr(end+1, std::string::npos);
    }
  }
  for (directory_list::reverse_iterator it = tmp_list.rbegin();
       it != tmp_list.rend(); it++) {
    directory_list_.push_front(*it);
  }
}

std::string
find_file(std::string name)
{
  // If file starts with "/", treat it as a full path
  if (name[0] == '/') {
    FILE* file = fopen(name.c_str(), "r");
    if (file) {
      // File opened OK: close and return the full path
      fclose(file);
      return name;
    }
  }
  else {
    // Search through directory list instead
    for (directory_list::iterator it = directory_list_.begin();
	 it != directory_list_.end(); it++) {
      std::string full_path = *it + "/" + name;
      FILE* file = fopen(full_path.c_str(), "r");
      if (file) {
	// File opened OK: close and return the full path 
	fclose(file);
	return full_path;
      }
    }
  }
  // File not found: return an empty string
  std::string empty_string;
  return empty_string;
}

std::string
search_directories()
{
  std::string dirs;
  for (directory_list::iterator it = directory_list_.begin();
       it != directory_list_.end(); it++) {
    dirs += *it + ":";
  }
  dirs.erase(dirs.size()-1);
  return dirs;
}

std::string
first_directory()
{
  return *(directory_list_.begin());
}
