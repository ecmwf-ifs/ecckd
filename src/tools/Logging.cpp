/// @file      Logging.cpp
/// @brief     Implements the functions and classes to handle logging messages
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_GMVECSL
#include "LogCAPI.h"
#define MAX_STRING_LEN 1000
#endif

#define STAR_LINE "**********************************************************************\n"

#include "Logging.h"

LogStream active_log_stream_(LOG_LEVEL_INFO);

//ECSIM-conformant prefixes to the different levels of logging
 // const char* log_prefix_str[] = {"Error   |",
 // 				"Warning |",
 // 				"Info    |",
 // 				"Progress|",
 // 				"Debug   |",
 // 				"Debug   |"};

const char* log_prefix_str[] = {"Debug   |",
				"Info    |",
				"Progress|",
				"Warning |",
				"Error   |",
				"Detail  |"};


LogStream::~LogStream() {
  if (is_file_) file_.close();
}

void set_log_file(const std::string& file_name) {
  active_log_stream_.file_.open(file_name.c_str());
  active_log_stream_.is_file_ = true;
}

void set_log_level(const std::string& log_level_string)
{
  if (log_level_string == "error") {
    active_log_stream_.log_level_ = LOG_LEVEL_ERROR;
  }
  else if (log_level_string == "warning") {
    active_log_stream_.log_level_ = LOG_LEVEL_WARNING;
  }
  else if (log_level_string == "info") {
    active_log_stream_.log_level_ = LOG_LEVEL_INFO;
  }
  else if (log_level_string == "progress") {
    active_log_stream_.log_level_ = LOG_LEVEL_PROGRESS;
  }
  else if (log_level_string == "detail") {
    active_log_stream_.log_level_ = LOG_LEVEL_DETAIL;
  }
  else if (log_level_string == "debug") {
    active_log_stream_.log_level_ = LOG_LEVEL_DEBUG;
  }
  else {
    ERROR << "Invalid log level string \"" << log_level_string
	  << "\"; must be one of error, warning, info, progress, debug";
    THROW(PARAMETER_ERROR);
  }
}

void
set_print_log_prefix(bool b)
{
  active_log_stream_.print_prefix_ = b;
}

void
set_logging_verbose(bool b)
{
  active_log_stream_.verbose_ = b;
}

#ifdef HAVE_GMVECSL
void
set_gmv_logging(bool b, int job_order_log_level)
{
  active_log_stream_.use_gmv_ = b;
  active_log_stream_.log_level_ = static_cast<LogLevel>(job_order_log_level);
}
int set_file_log_level(const std::string& log_level_string)
{
  int file_log_level;
  if (log_level_string == "error") {
    file_log_level = LOG_LEVEL_ERROR;
  }
  else if (log_level_string == "warning") {
    file_log_level = LOG_LEVEL_WARNING;
  }
  else if (log_level_string == "info") {
    file_log_level = LOG_LEVEL_INFO;
  }
  else if (log_level_string == "progress") {
    file_log_level = LOG_LEVEL_PROGRESS;
  }
  else if (log_level_string == "detail") {
    file_log_level = LOG_LEVEL_DETAIL;
  }
  else if (log_level_string == "debug") {
    file_log_level = LOG_LEVEL_DEBUG;
  }
  else {
    ERROR << "Invalid log level string for log file \"" << log_level_string
	  << "\"; must be one of error, warning, info, progress, debug";
    THROW(PARAMETER_ERROR);
  }
  return file_log_level;
}
#else
void
set_gmv_logging(bool b)
{
  ERROR << "GMV logging is not available";
  THROW(WRITE_ERROR);
}

#endif


//   char str[MAX_STRING_LEN];
//   std::string line;
//   do {
//     std::getline(ss_, line, MAX_STRING_LEN);
//     int len = line.size();
//     if (len > 0) {
//       if (log_level_ > LOG_LEVEL_INFO) {
// 	LOG_DEBUG(line.c_str());
//       }
//       else if (log_level_ == LOG_LEVEL_INFO) {
// 	LOG_INFO(line.c_str());
//       }
//       else if (log_level_ == LOG_LEVEL_WARNING) {
// 	LOG_WARNING(line.c_str());
//       }
//       else {
// 	LOG_ERROR(line.c_str());
//       }
//     }
//   } while (len > 0);
// }

void
LogStream::flush()
{
#ifdef HAVE_GMVECSL
  if (use_gmv_) {
    // if (last_level_ == LOG_LEVEL_ERROR) {
    //   std::stringstream s;
    //   s << "*** Error at line " << source_line_ << " of " << source_file_ << ":";
    //   log_error(source_file_, source_func_, source_line_, s.str().c_str());
    // }
    std::string str;
    while (getline(ss_, str)) {
      if (last_level_ == LOG_LEVEL_DEBUG) {
	log_debug(source_file_, source_func_, source_line_, str.c_str());
      }
      else if (last_level_ == LOG_LEVEL_INFO) {
	log_info(source_file_, source_func_, source_line_, str.c_str());
      }
      else if (last_level_ == LOG_LEVEL_DETAIL && active_log_stream_.verbose_) {
	log_info(source_file_, source_func_, source_line_, str.c_str());
      }
      else if (last_level_ == LOG_LEVEL_WARNING) {
	log_warning(source_file_, source_func_, source_line_, str.c_str());
      }
      else if (last_level_ == LOG_LEVEL_ERROR) {
	log_error(source_file_, source_func_, source_line_, str.c_str());
      }
    }
  } 
  else
#endif
    {
      std::ostream& stream = is_file_ ? file_ : std::cout;
      const char* prefix = "";
      if (print_prefix_) {
	prefix = log_prefix_str[static_cast<int>(last_level_)];
      }
      else if (last_level_ == LOG_LEVEL_ERROR) {
	prefix = "*** ";
      }

      // Opening info
      std::string str;
      if (last_level_ == LOG_LEVEL_ERROR) {
	getline(ss_, str);
	if (!str.empty()) {
	  stream << prefix << STAR_LINE;
	  stream << prefix << "*** Error at line " << source_line_
		 << " of " << source_file_ << ": " << str << "\n";
	  while (getline(ss_, str)) {
	    stream << prefix << "*** " << str << "\n";
	  }
	  // Closing info
	  stream << prefix << STAR_LINE;
	}
      }
      else {
	if (last_level_ == LOG_LEVEL_WARNING) {
	  getline(ss_, str);
	  if (!str.empty()) {
	    stream << prefix << "Warning at line " << source_line_
		   << " of " << source_file_ << ": " << str << "\n";
	  }
	}

	while (getline(ss_, str)) {
	  stream << prefix << str << "\n";
	}
      }

      if (!is_file_) {
	std::cout.flush();
      }
    }
    // Clear the string stream
  ss_.str(std::string());
  ss_.clear();
  //  last_level_ = LOG_LEVEL_DEBUG;
}

LogStream&
log_stream(LogLevel log_level, const char* source_file, 
	   const char* source_func, int source_line)
{
  if (log_level == LOG_LEVEL_DETAIL && active_log_stream_.verbose_){
    log_level = LOG_LEVEL_INFO;
  }else if (log_level == LOG_LEVEL_DETAIL && !active_log_stream_.verbose_){
    log_level = LOG_LEVEL_DEBUG;
  }
  if (static_cast<int>(log_level)
      >= static_cast<int>(active_log_stream_.log_level_)) {
    if (log_level != active_log_stream_.last_level_) {
      active_log_stream_.flush();
    }
    active_log_stream_.source_file_ = source_file;
    active_log_stream_.source_func_ = source_func;
    active_log_stream_.source_line_ = source_line;
    active_log_stream_.last_level_  = log_level;
    active_log_stream_.active_line_ = true;
  }
  else {
    active_log_stream_.active_line_ = false;
  }
  return active_log_stream_;
}

void
flush_stream(LogLevel log_level)
{
  /*
  if (static_cast<int>(log_level)
      <= static_cast<int>(active_log_stream_.log_level_)) {
    active_log_stream_.flush();
  }
  */
  active_log_stream_.flush();
}

// Print ECSIM progress bar
void progress_bar(int i, int n) {
#ifdef HAVE_GMVECSL
  if (active_log_stream_.use_gmv_) {
    int prog = (100*i)/n;
    log_progress(prog);
  }
  else 
#endif
    if (active_log_stream_.print_prefix_
	&& active_log_stream_.log_level_ >= LOG_LEVEL_INFO) {
      int prog = (100*i)/n;
      log_stream(LOG_LEVEL_PROGRESS,"","",0) << prog << "\n";
    }
}
