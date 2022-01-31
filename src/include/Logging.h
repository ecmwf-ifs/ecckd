/// @file      Logging.h
/// @brief     Functions and classes to handle logging messages
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef Logging_H
#define Logging_H 1

#include <string>
#include <sstream>
#include <cstring>
//#include <iostream>
//#include <fstream>
#include <typeinfo>

#include "Error.h"

// ECSIM-conformant prefixes to the different levels of logging
extern const char* log_prefix_str[];

/* enum LogLevel { */
/*   LOG_LEVEL_ERROR = 0, */
/*   LOG_LEVEL_WARNING, */
/*   LOG_LEVEL_INFO, */
/*   LOG_LEVEL_PROGRESS, */
/*   LOG_LEVEL_DETAIL, */
/*   LOG_LEVEL_DEBUG */
/* }; */

enum LogLevel {
  LOG_LEVEL_DEBUG = 0,
  LOG_LEVEL_INFO,
  LOG_LEVEL_PROGRESS,
  LOG_LEVEL_WARNING,
  LOG_LEVEL_ERROR,
  LOG_LEVEL_DETAIL=5,
};

// A single object of this type should be created at the start of the
// main function - its destructor ensures cleanup of outputs
// struct LogDemon {
//   LogDemon(const char* processor_name, const char* processor_version,
// 	      int log_level = 3);
//   ~LogDemon();
// };


struct LogStream {
  LogStream(LogLevel log_level)
    : source_file_(0), source_func_(0), source_line_(0),
      log_level_(log_level), last_level_(LOG_LEVEL_DEBUG), 
      is_file_(false), line_complete_(true), use_gmv_(false),
      print_prefix_(false), active_line_(false), verbose_(false) { }
  ~LogStream();

  void flush();

  // Data
  std::ofstream file_;
  const char* source_file_;
  const char* source_func_;
  int source_line_;
  LogLevel log_level_;
  LogLevel last_level_;
  bool is_file_;
  bool line_complete_;
  bool use_gmv_;
  bool print_prefix_;
  bool active_line_;
  bool verbose_; 
  std::stringstream ss_;
};

void set_log_file(const std::string& file_name);
void set_log_level(const std::string& log_level_string);
int set_file_log_level(const std::string& log_level_string);
void set_print_log_prefix(bool b = true);
void set_gmv_logging(bool b = true, int job_order_log_level = LOG_LEVEL_INFO);
void set_logging_verbose(bool b = false);

template <typename T>
LogStream& operator<<(LogStream& s, T const & val) {
  if (s.active_line_) {
    s.ss_ << val;
  }
  return s;
}

// If we see "\n" then flush to the stream immediately
inline
LogStream& operator<<(LogStream& s, const char val[]) {
  if (s.active_line_) {
    s.ss_ << val;
    int len = std::strlen(val);
    if (val[len-1] == '\n') {
      s.flush();
    }
  }
  return s;
}

LogStream& log_stream(LogLevel log_level, const char* source_file,
		      const char* source_func, int source_line);
void flush_stream(LogLevel log_level);

#define ERROR log_stream(LOG_LEVEL_ERROR, __FILE__, __FUNCTION__, __LINE__)
#define ENDERROR flush_stream(LOG_LEVEL_ERROR)

#define THROW(object) flush_stream(LOG_LEVEL_ERROR);	      \
  if (trace_exceptions_) { throw DebugPreserveStackTrace(); } \
  else { throw object; }

#define STAR_LINE "**********************************************************************\n"

#define WARNING log_stream(LOG_LEVEL_WARNING, __FILE__, __FUNCTION__, __LINE__)
#define ENDWARNING flush_stream(LOG_LEVEL_WARNING)
#define LOG log_stream(LOG_LEVEL_INFO,        __FILE__, __FUNCTION__, __LINE__)
#define DETAIL log_stream(LOG_LEVEL_DETAIL,   __FILE__, __FUNCTION__, __LINE__)
#define DEBUG log_stream(LOG_LEVEL_DEBUG,     __FILE__, __FUNCTION__, __LINE__)



inline
int num_digits(int x) {  
  return (x < 10 ? 1 :   
	  (x < 100 ? 2 :   
	   (x < 1000 ? 3 :   
	    (x < 10000 ? 4 :   
	     (x < 100000 ? 5 :   
	      (x < 1000000 ? 6 :   
	       (x < 10000000 ? 7 :  
		(x < 100000000 ? 8 :  
		 (x < 1000000000 ? 9 :  
		  10)))))))));  
}

inline
std::string pad_right(int x, int ndigits) {
  std::stringstream s;
  s << x;
  for (int i = num_digits(x); i < ndigits; ++i) {
    s << " ";
  }
  return s.str();
}

// Print ECSIM progress bar
void progress_bar(int i, int n);


#endif
