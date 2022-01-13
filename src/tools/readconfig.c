/* readconfig.c -- Functions to read configuration information from a file

   Copyright (C) 2003 The Universtiy of Reading
   Author: Robin Hogan <r.j.hogan@reading.ac.uk>

   This software is licensed under the terms of the Apache Licence
   Version 2.0 which can be obtained at
   http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
//#include <stdio.h>

#include "readconfig.h"

#define COMMENT_CHAR '#'
#define SECTION_CHAR '.'
#define ESCAPE_CHAR '\\'
#define REFERENCE_CHAR '$'
#define FILE_STRING ""
#define BEGIN_COMMAND "\\begin"
#define END_COMMAND "\\end"
#define INCLUDE_COMMAND "\\include"
// Uncomment in order for values in output of rc_sprint to be wrapped
// in either quotation marks or curly brackets
//#define SPRINT_CLASSIC 1
#ifdef SPRINT_CLASSIC
  // Classic behaviour writes each entry on a new line with no whitespace suppression
  #define SPRINT_EQUALS ' '
  //#define SPRINT_COMPRESS_WHITESPACE 1
  #define SPRINT_SEPARATOR "\n"
  #define SPRINT_STRING_BEGIN '"'
  #define SPRINT_STRING_END '"'
  #define SPRINT_NUMBER_BEGIN '{'
  #define SPRINT_NUMBER_END '}'
  //#define SPRINT_WRAP 1
#else
  // Alternatively put everything on a new line
  #define SPRINT_EQUALS '='
  #define SPRINT_COMPRESS_WHITESPACE
  #define SPRINT_SEPARATOR "; "
  #define SPRINT_STRING_BEGIN '{'
  #define SPRINT_STRING_END '}'
  #define SPRINT_NUMBER_BEGIN '{'
  #define SPRINT_NUMBER_END '}'
//#define SPRINT_WRAP 1
#endif

#define REPLACE_VALUE(dest, new) \
  if (dest) { free(dest); }	 \
  if (new) { dest = new; }	 \
  else { dest = strdup("1"); }

static char *memory_error 
= "Error allocating memory for configuration information\n";

#define ERROR(message)							\
  if (err_file) fprintf(err_file, "*** Error at line %d of %s: %s\n",	\
			__LINE__, __FILE__, message)

/* Move file pointer to the start of the next line, returning '\n' on
   success or EOF if the file ended first. */
static
int
__rc_skip_line(FILE *file)
{
  char c;
  do {
    c = fgetc(file);
  } while (c != '\n' && c != EOF);
  return c;
}

/* Move the file pointer over whitespace, returning the first
   non-whitespace character found, or '\n' if the line ended first, or
   EOF if the file ended first. */
static
int
__rc_skip_whitespace(FILE *file)
{
  int c;
  do {
    c = fgetc(file);
  } while (c <= ' ' && c != '\n' && c != EOF);
  return c;
}

/* Remove trailing whitespace from a string */
static
void
__rc_remove_trailing_whitespace(char* value)
{
  if (value) {
    char *ch = value + strlen(value) - 1;
    while (ch >= value && *ch <= ' ') {
      *ch = '\0';
      ch--;
    }
  }
}

static
void
__rc_copy_compress_whitespace(char* dest, const char* src)
{
  const char* s = src;
  char* d = dest;
  int space_required = 0;
  int in_word = 0;
  while(*s != '\0') {
    if (*s > ' ') {
      /* Non-whitespace found */
      if (!in_word && space_required) {
	*d = ' ';
	++d;
	space_required = 0;
	in_word = 1;
      }
      *d = *s;
      ++d;
      in_word = 1;
    }
    else if (in_word) {
      /* Whitespace found after word */
      space_required = 1;
      in_word = 0;
    }
    ++s;
  }
  *d = '\0';
}

/* Find param in data, returning an element of the rc_data list, or
   NULL if not present. Note that the search is case insensitive. */
static
rc_data *
__rc_find_no_section(rc_data *data, const char *param)
{
  if (data->param) {
    if (strcasecmp(param, data->param) == 0) {
      return data;
    }
    else if (data->next) {
      return __rc_find_no_section(data->next, param);
    }
  }
  return NULL;
}

static
rc_data *
__rc_find(rc_data *data, const char *param)
{
  if (!data->section_reqd) {
    return __rc_find_no_section(data, param);
  }
  else if (data->param) {
    int section_len = strlen(data->section_reqd);
    rc_data *current = data;
    while (current->next) { // Is this right?
      if (strncasecmp(data->section_reqd, current->param,
		      section_len) == 0) {
	/* We are in the right section... */
	if (current->param[section_len] == SECTION_CHAR
	    && strcasecmp(current->param+section_len+1,
			  param) == 0) {
	  /* Found a match */
	  return current;
	}
      }
      current = current->next;
    }
  }
  return NULL;
}

static
int
__rc_register_table(rc_data *data, const char *params, const char *value,
		    FILE *err_file)
{
  rc_data **data_list = NULL;
  rc_data *data_end = NULL;

  const char *c = params;
  const char *param_tail = params;
  char *param;
  char *section_reqd = data->section_reqd;
  int section_len = 0;
  int list_len = 0;
  int ilist = 0;
  int old_value_len = 0;
  int new_value_len = 0;

  if (section_reqd) {
    section_len = strlen(section_reqd);
  }

  /* Count how many parameters, allocating each one */
  while (*param_tail != '\0') {
    rc_data *data_current;
    int param_length = 0;
    int m = 0, n = 0;

    /* Skip whitespace */
    while (*param_tail && *param_tail <= ' ') {
      param_tail++;
    }
    if (*param_tail == '\0') {
      break;
    }

    /* Found the start of a name: find the end but also check if this
       is a vector or matrix and read in the dimensions */
    c = param_tail;
    while (*c > ' ' && *c != '[') {
      c++;
    }
    param_length = c - param_tail;
    if (*c == '[') {
      c++;
      if (*c) {
	char *newc;
	int val = strtol(c, &newc, 10);
	if (val > 0) {
	  m = val;
	}
	if (newc && newc > c) {
	  c = newc;
	  if (*c == ']') {
	    c++;
	    if (*c == '[') {
	      c++;
	      if (*c) {
		int val = strtol(c, &newc, 10);
		if (val > 0) {
		  n = val;
		}
	      }
	    }
	  }
	}
      }
      /* Skip over remaining characters */
      while (*c && *c > ' ') {
	c++;
      }
    }
    
    if (section_len) {
      /* Add section to param in the form "section.param" */
      param = malloc(section_len + param_length + 2);
      if (!param) {
	ERROR(memory_error);
	return 0;
      }
      strncpy(param, section_reqd, section_len);
      param[section_len] = SECTION_CHAR;
      strncpy(param+section_len+1, param_tail, param_length);
      param[section_len + param_length + 1] = '\0';
    }
    else {
      //      param = strdup(param_tail);
      param = malloc(param_length + 1);
      if (!param) {
	ERROR(memory_error);
	return 0;
      }
      strncpy(param, param_tail, param_length);
      param[param_length] = '\0';      
    }

    /* Search for param in existing data */
    data_current = __rc_find_no_section(data, param);
    if (data_current) {
      /* Existing element found: delete current value */
      if (data_current->value) {
	free(data_current->value);
	data_current->value = NULL;
      }
      data_current->m = m;
      data_current->n = n;
      /* Don't need storage of param */
      free(param);
      param = NULL;
    }
    else {
      /* Add a new list element */
      if (!data_end) {
	/* Find the end of the list */
	data_end = data;
	while (data_end->next) {
	  data_end = data_end->next;
	}
      }
      /* Create a new data element */
      data_end->next = malloc(sizeof(rc_data));
      if (!data_end->next) {
	ERROR(memory_error);
	return 0;
      }
      else {
	data_end->next->param = NULL;
	data_end->next->value = NULL;
	data_end->next->section_reqd = NULL;
	data_end->next->next = NULL;
	data_end->next->m = 0;
	data_end->next->n = 0;
      }
      data_current = data_end;
      data_end = data_end->next;
      data_current->param = param;
      data_current->m = m;
      data_current->n = n;
    }

    /* Add to local list */
    data_list = realloc(data_list, sizeof(*data_list)*(++list_len));
    if (!data_list) {
      ERROR(memory_error);
      return 0;
    }
    data_list[list_len-1] = data_current;
    
    /* Prepare for next cycle */
    param_tail = c;
  } // end of while loop
    
  if (!list_len) {
    /* No param names found in brackets */
    if (err_file) {
      fprintf(err_file, "No param names found in brackets of table assignment\n");
    }
    return 0;
  }
  
  /* Read through "value" */
  while (*value != '\0') {
    /* Skip whitespace */
    while (*value <= ' ' && *value != '\0') {
      value++;
    }
    if (*value == '\0') {
      /* Found terminating character: quit */
      break;
    }

    /* Found the start of a value: find the end */
    c = value;
    if (*c == '{') {
      c++;
      /* Curly braces: find the terminating brace */
      while (*c && *c != '}') {
	c++;
      }
      if (*c) {
	c++;
      }
    }
    else if (*c == '"') {
      c++;
      /* Quotes: find the terminating quote */
      while (*c && *c != '"') {
	c++;
      }
      if (*c) {
	c++;
      }
    }
    else {
      while (*c > ' ') {
	c++;
      }
    }
    /* Found the end: add to value */
    if (data_list[ilist]->value) {
      old_value_len = strlen(data_list[ilist]->value);
    }
    else {
      old_value_len = 0;
    }
    new_value_len = c-value;
    //    data_list[ilist]->m = data_list[ilist]->n = 0;
    data_list[ilist]->value = realloc(data_list[ilist]->value,
					 old_value_len + new_value_len + 2);
    if (!data_list[ilist]->value) {
      ERROR(memory_error);
      return 0;
    }
    data_list[ilist]->value[old_value_len] = ' ';
    strncpy(data_list[ilist]->value + old_value_len + 1, value, new_value_len);
    data_list[ilist]->value[old_value_len + 1 + new_value_len] = '\0';

    /* Increase counter */
    if (++ilist >= list_len) {
      ilist = 0;
    }
    /* Move value to next whitespace */
    value = c;
  }

  /* Clear allocated memory */
  free(data_list);

  return 1;
}


/* Add a param-value pair to an existing rc_data structure,
   overwriting an existing param with the same name (case
   insensitive). Note that param and value are the actual pointers
   added to the structure so they should not be freed after calling
   this function. Return 1 on success, 0 if memory allocation of the
   next element in the list failed. */
static
int
__rc_register(rc_data *data, char *param, char *value)
{
  char *c = param;
  char *section_reqd = data->section_reqd;
  int vector_dim_offset = 0;
  int m = 0, n = 0;
  int section_len = 0;

  if (section_reqd) {
    section_len = strlen(section_reqd);
  }

  /* Check if this is a vector and read in the dimensions */
  while (*c) {
    if (*c == '[') {
      vector_dim_offset = c - param;
      *c = '\0';
      break;
    }
    c++;
  }
  if (vector_dim_offset) {
    c++;
    if (*c) {
      char *newc;
      int val = strtol(c, &newc, 10);
      if (val > 0) {
	m = val;
      }
      if (newc && newc > c) {
	c = newc;
	while (*c && *c != '[') {
	  c++;
	}
	c++;
	if (*c) {
	  int val = strtol(c, &newc, 10);
	  if (val > 0) {
	    n = val;
	  }
	}
      }
    }
  }

  /* Check whether "value" is actually a reference to another
     variable, and if so, substitute it */
  if (value && value[0] == REFERENCE_CHAR) {
    rc_data *subst_data = __rc_find_no_section(data, value+1);
    if (subst_data) {
      char* tmp_value = strdup(subst_data->value);
      REPLACE_VALUE(value, tmp_value);
    }
    else {
      // Replacement value not found: the variable will not be
      // assigned. This means that any existing value will still be
      // present.
      //      return 1;
    }
  }

  /* Loop through data structure, looking for existing parameter */
  while (data->param) {
    /* If a "section" is specified, then the parameter will be stored
       as section.param, so we need to test both parts */
    if (section_len) {
      if (strncasecmp(section_reqd, data->param, section_len) == 0
	  && data->param[section_len] == SECTION_CHAR
	  && strcasecmp(data->param+section_len, param) == 0) {
	REPLACE_VALUE(data->value, value);
        free(param);
	return 1;
      }
    }
    else if (strcasecmp(param, data->param) == 0) {
      REPLACE_VALUE(data->value, value);
      free(param);
      return 1;
    }
    data = data->next;
  }

  if (section_len) {
    data->param = malloc(sizeof(char)*(strlen(param)
					  +section_len+2));
    if (!data->param) {
      return 0;
    }
    strncpy(data->param, section_reqd, section_len);
    data->param[section_len] = SECTION_CHAR;
    strcpy(data->param+section_len+1, param);
    free(param);
  }
  else {
    data->param = param;
  }
  REPLACE_VALUE(data->value, value);
  data->section_reqd = NULL;
  data->m = m;
  data->n = n;
  data->next = malloc(sizeof(rc_data));
  if (!data->next) {
    return 0;
  }
  else {
    data->next->param = data->next->value = NULL;
    data->next->section_reqd = NULL;
    data->next->next = NULL;
    data->next->m = data->next->n = 0;
    return 1;
  }
}


/* Free an rc_data structure created with rc_read() - rc_data is a
   singly-linked list and its contents are recursively detelted. */
void
rc_clear(rc_data *data)
{
  if (data->param) {
    free(data->param);
  }
  if (data->value) {
    free(data->value);
  }
  if (data->section_reqd) {
    free(data->section_reqd);
  }
  if (data->next) {
    rc_clear(data->next);
  }
  free(data);
}

/* Append contents of file to an existing data structure */
int
rc_append(rc_data *data, const char *file_name, FILE *err_file)
{
  FILE *file;
  int c;

  if (!data) {
    if (err_file) {
      fprintf(err_file, "No data specified in call to rc_append()\n");
    }
    return 0;
  }

  if (!file_name) {
    if (err_file) {
      fprintf(err_file, "No file name specified in call to rc_append()\n");
    }
    return 0;
  }

  file = fopen(file_name, "r");

  if (!file) {
    if (err_file) {
      fprintf(err_file, "Error opening %s\n", file_name);
    }
    return 0;
  }

  while (!feof(file)) {
    char *param = NULL, *value = NULL;
    int param_length = 0, value_length = 0;
    /* Skip whitespace */
    c = __rc_skip_whitespace(file);
    if (c == EOF) {
      break;
    }
    else if (c == COMMENT_CHAR) {
      __rc_skip_line(file);
      continue;
    }
    else if (c == '\n') {
      continue;
    }
    else if (c == '(') {
      /* Table of objects */
      while (c != ')') {
	if (c == EOF) {
	  if (err_file) {
	    fprintf(err_file,
		    "File ended before table column names finished: \"");
	    fprintf(err_file, param);
	    fprintf(err_file, "\"\n");
	    fclose(file);
	    free(param);
	    return 0;
	  }
	}
	param = realloc(param, (++param_length)+1);
	if (!param) {
	  ERROR(memory_error);
	  fclose(file);
	  return 0;
	}
	param[param_length-1] = c;
	c = fgetc(file);
      }
      /*      ungetc(c, file);*/
      param[param_length] = '\0';
    }
    else {
      /* Read param name */
      while (c > ' ' && c != COMMENT_CHAR && c != EOF) {
	param = realloc(param, (++param_length)+1);
	if (!param) {
	  ERROR(memory_error);
	  fclose(file);
	  return 0;
	}
	param[param_length-1] = c;
	c = fgetc(file);
      }
      ungetc(c, file);
      param[param_length] = '\0';
    }

    /* Skip whitespace */
    c = __rc_skip_whitespace(file);
    if (c == COMMENT_CHAR) {
      __rc_skip_line(file);
    }
    else if (c != '\n') {
      /* Read value */
      if ((c == '\'') || (c == '"')) {
	int quote = c;
	c = fgetc(file);
	while (c != EOF && c != quote) {
	  value = realloc(value, (++value_length)+1);
	  if (!value) {
	    ERROR(memory_error);
	    fclose(file);
	    return 0;
	  }
	  value[value_length-1] = c;
	  value[value_length] = '\0';
	  c = fgetc(file);
	}
      }
      else if (c == '{') {
	c = fgetc(file);
	while (c != EOF && c != '}') {
	  if (c == COMMENT_CHAR) {
	    __rc_skip_line(file);
	  }
	  else {
	    value = realloc(value, (++value_length)+1);
	    if (!value) {
	      ERROR(memory_error);
	      fclose(file);
	      return 0;
	    }
	    value[value_length-1] = c;
	    value[value_length] = '\0';
	  }
	  c = fgetc(file);
	}
      }
      else {
	while (c != EOF && c != '\n') {
	  if (c == COMMENT_CHAR) {
	    __rc_skip_line(file);
	    break;
	  }
	  else if (c != '\r') {
	    value = realloc(value, (++value_length)+1);
	    if (!value) {
	      ERROR(memory_error);
	      fclose(file);
	      return 0;
	    }
	    value[value_length-1] = c;
	    value[value_length] = '\0';
	  }
	  c = fgetc(file);
	}
	/* Remove trailing whitespace */
	__rc_remove_trailing_whitespace(value);
      }
    }

    /* Check for special commands */
    if (param[0] == ESCAPE_CHAR) {
      if (strcasecmp(param,BEGIN_COMMAND) == 0) {
	/* \begin command starts a new section */
	if (data->section_reqd) {
	  int len = strlen(data->section_reqd);
	  char *section_reqd
	    = realloc(data->section_reqd, len+value_length+2);
	  if (!section_reqd) {
	    ERROR(memory_error);
	    fclose(file);
	    return 0;
	  }
	  data->section_reqd = section_reqd;
	  section_reqd[len] = SECTION_CHAR;
	  strcpy(section_reqd+len+1, value);
	  free(param);
	  free(value);
	  param = NULL;
	  value = NULL;
	}
	else {
	  data->section_reqd = value;
          free(param);
          param = NULL;
	}
      }
      else if (strcasecmp(param,END_COMMAND) == 0) {
	/* \end command ends a section */
	if (data->section_reqd) {
	  int len = strlen(data->section_reqd);
	  char* last_section = data->section_reqd+len;
	  /* Find section separator */
	  do {
	    last_section--;
	  } while (last_section[0] != SECTION_CHAR
		   && last_section > data->section_reqd);
	  if (last_section[0] == SECTION_CHAR) {
	    last_section++;
	  }
	  /* If \end has an argument, check it agrees */
	  if (value) {
	    if (strcasecmp(last_section, value) != 0) {
	      if (err_file) {
		fprintf(err_file, "\"%s %s\" ended by \"%s %s\"\n",
			BEGIN_COMMAND, last_section,
			END_COMMAND, value);
	      }
	      fclose(file);
	      free(value);
	      free(param);
	      return 0;
	    }
	    free(value);
	    value = NULL;
	  }
	  if (last_section == data->section_reqd) {
	    /* Only one section in force: remove */
	    free(data->section_reqd);
	    data->section_reqd = NULL;
	  }
	  else {
	    /* More than one section in force: remove last one */
	    last_section[-1] = '\0';
	  }
	}
	else {
	  if (err_file) {
	    fprintf(err_file, "\"%s\" with no \"%s\"\n",
		    END_COMMAND, BEGIN_COMMAND);
	  }
	  if (value) {
	    free(value);
	  }
	  fclose(file);
	  return 0;
	}
      }
      else if (strcasecmp(param,INCLUDE_COMMAND) == 0) {
	/* \include config data from another file */
	if (!value) {
	  if (err_file) {
	    fprintf(err_file, "%s does not specify a file\n",
		    INCLUDE_COMMAND);
	  }
	  free(param);
	  fclose(file);
	  return 0;
	}
	if (value[0] != '/') {
	  /* Interpret value as a filename relative to the directory
	     of the current file; find root directory from the current
	     file */
	  int i = strlen(file_name)-1;
	  while (file_name[i] != '/' && i > 0) {
	    i--;
	  }
	  if (file_name[i] == '/') {
	    /* Need to combine directory and filename */
	    char *new_value = malloc(sizeof(char)*(i+2+value_length));
	    if (!new_value) {
	      ERROR(memory_error);
	      free(param);
	      free(value);
	      fclose(file);
	      return 0;
	    }
	    strncpy(new_value, file_name, i+1);
	    strcpy(new_value+i+1, value);
	    free(value);
	    value = new_value;
	  }
	}
	/* Check we are not reading the same file leading to infinite
	   recursion */
	if (strcmp(value, file_name) == 0) {
	  if (err_file) {
	    fprintf(err_file, "%s attempts to %s itself\n", file_name,
		    INCLUDE_COMMAND);
	  }
	  free(param);
	  free(value);
	  fclose(file);
	  return 0;
	}
	/* Read data from other file */
	if (rc_append(data, value, err_file) == 0) {
	  free(value);
	  free(param);
	  fclose(file);
	  return 0;
	}
	free(value);
	free(param);
	value = NULL;
	param = NULL;
      }
      else {
	/* A command found that we don't understand: ignore... */
	free(param);
	param = NULL;
	if (value) {
	  free(value);
	  value = NULL;
	}
      }
      if (param) {
        free(param);
        param = NULL;
      }
    }
    /* The first character is not '\', so this is a normal line of
       data */
    else if (param[0] == '(') {
      /* This is a table of numbers */
      if (!__rc_register_table(data, param+1, value, err_file)) {
	if (err_file) {
	  fprintf(err_file, "Error assigning table of values\n");
	}
	fclose(file);
	return 0;
      }
      free(param);
      free(value);
    }
    else if (!__rc_register(data, param, value)) {
      /* Register result */
      ERROR(memory_error);
      fclose(file);
      return 0;
    }
  }
  fclose(file);

  return 1;
}

/* Add a param-value pair to an existing rc_data structure,
   overwriting an existing param with the same name (case
   insensitive). Returns 1 on success and 0 on failure. */
int
rc_register(rc_data *data, const char *param, const char *value)
{
  char *newparam = strdup(param);
  if (value) {
    char *newvalue = strdup(value);
    return __rc_register(data, newparam, newvalue);
  }
  else {
    return __rc_register(data, newparam, NULL);
  }
}

/* Search the command-line arguments for param=value pairs and -param
   arguments and add them to the rc_data structure. */
int
rc_register_args(rc_data *data, int argc, const char **argv)
{
  int i;

  for (i = 1; i < argc; i++) {
    /* Find an "=" sign in argument i */
    if (argv[i][0] == '-' && argv[i][1]) {
      char *param = strdup(argv[i]+1);
      if (!param) {
	return 0;
      }
      if (!__rc_register(data, param, NULL)) {
	return 0;
      }
    }
    else {
      const char *c = argv[i];
      while (*c != '\0') {
	if (*c == '=') {
	  /* Found one */
	  int param_length = c-argv[i];
	  char *param = malloc(param_length+1);
	  char *value;
	  if (c[1] == REFERENCE_CHAR) {
	    /* Found a reference to another variable: substitute */
	    rc_data *subst_data = __rc_find_no_section(data, c+2);
	    if (subst_data) {
	      value = strdup(subst_data->value);
	    }
	    else {
	      /* Referenced variable not found: skip */
	      free(param);
	      continue;
	    }
	  }
	  else {
	    value = strdup(c+1);
	  }
	  if (!value || !param) {
	    return 0;
	  }
	  strncpy(param, argv[i], param_length);
	  param[param_length] = '\0';
	  if (!__rc_register(data, param, value)) {
	    return 0;
	  }
	}
	c++;
      }
    }
  }
  return 1;
}

/* Search the command-line arguments for files: ones that are *not*
   param=value pairs or -param arguments, and add them to the rc_data
   structure as FILE0, FILE1 etc. */
int
rc_register_files(rc_data *data, int argc, const char **argv)
{
  int i;
  int nfiles = 0;

  {
    char *param = malloc(2+strlen(FILE_STRING));
    char *value = strdup(argv[0]);
    sprintf(param, "%s0", FILE_STRING,0);

    /* Register argv[0], the name of the executable */
    if (!__rc_register(data, param, value)) {
      return 0;
    }
  }

  nfiles++;

  for (i = 1; i < argc; i++) {
    /* Check arg doesn't start with "-" */
    if (argv[i][0] == '-' && argv[i][1]) {
      continue;
    }
    else {
      /* Look for an "=" sign in argument i */
      char found_file = 1;
      const char *c = argv[i];
      while (*c != '\0') {
	if (*c == '=') {
	  found_file = 0;
	  break;
	}
	c++;
      }
      if (found_file) {
	int str_len = 10+strlen(FILE_STRING);
	char *param = malloc(str_len);
	char *value = strdup(argv[i]);
	snprintf(param, str_len, "%s%d", FILE_STRING, nfiles);
	nfiles++;
	if (!__rc_register(data, param, value)) {
	  return 0;
	}
      }
    }
  }
  return 1;
}


/* Read configuration information from file called file_name and
   return a pointer to the rc_data structure, or NULL if an error
   occurred.  If file_name is NULL then an empty rc_data structure is
   returned. If err_file is not NULL, errors messages will be written
   to err_file. */
rc_data *
rc_read(const char *file_name, FILE *err_file)
{
  rc_data *data;
  int status;

  /* Initialize an empty data structure */
  data = malloc(sizeof(rc_data));
  if (!data) {
    ERROR(memory_error);
    return NULL;
  }

  data->param = NULL;
  data->value = NULL;
  data->next = NULL;
  data->section_reqd = NULL;
  data->m = data->n = 0;

  if (!file_name) {
    /* Return empty data structure */
    return data;
  }

  status = rc_append(data, file_name, err_file);
  if (status) {
    if (data->section_reqd) {
      if (err_file) {
	fprintf(err_file, "Section \"%s\" unterminated by %s\n",
		data->section_reqd, END_COMMAND);
      }
      free(data->section_reqd);
      data->section_reqd = NULL;
      return NULL;
    }
    else {
      return data;
    }
  }
  else {
    rc_clear(data);
    return NULL;
  }
}


/* Find the first config file on the command line, returning an index
   to the argv array, and 0 if none is found. Basically the first
   argument that contains no "=" sign and doesn't start with "-" is
   returned. However, after a "--" argument, a filename begining with
   "-" could in principle be returned. */
int
rc_get_file(int argc, const char **argv)
{
  int i;
  int ignore_hyphen = 0;
  for (i = 1; i < argc; i++) {
    /* Find an "=" sign in argument i */
    const char *c = argv[i];

    if (!ignore_hyphen && c[0] == '-') {
      if (strcmp(c, "--")) {
	ignore_hyphen = 1;
      }
      continue;
    }

    while (*c != '\0') {
      if (*c == '=') {
	break;
      }
      c++;
    }
    if (*c != '=') {
      if(strstr(argv[i], ".cfg") != NULL){
	return i;
      }
    }
  }
  return 0;
}

/* Print contents of rc_data structure to file. */
void
rc_print(rc_data *data, FILE *file)
{
  char start_quote = '"', end_quote = '"';
  while (data->param) {
    fprintf(file, "%s", data->param);
    if (data->m > 0) {
      fprintf(file, "[%d]", data->m);
      if (data->n > 0) {
	fprintf(file, "[%d]", data->n);
      }
      start_quote = '{'; end_quote = '}';
    }
    if (data->value) {
      fprintf(file, " %c%s%c\n", start_quote, data->value, end_quote);
    }
    else {
      fprintf(file, " (no value)\n");
    }
    if (!(data = data->next)) {
      break;
    }   
  }
}

/* Return the contents of an rc_data structure as a string. Returns
   NULL on failure. Free with rc_free(). */
char *
rc_sprint(rc_data *data)
{
  int length = 0;
  char *out = NULL;
  if (!data) {
    return NULL;
  }

  while (data->param) {
    int param_length = strlen(data->param);
    out = realloc(out, length+param_length+strlen(SPRINT_SEPARATOR)+3);
    if (!out) {
      return NULL;
    }
    /* Add a separator between param-value pairs */
    if (length > 0) {
      strcpy(out+length, SPRINT_SEPARATOR);
      length += strlen(SPRINT_SEPARATOR);
    }
    /* Write the parameter name */
    strcpy(out+length, data->param);
    length += (param_length+1);

    if (data->value) {
      int value_length = strlen(data->value);
      char *val = data->value;
      int matsize_length = 0;
      if (data->m > 0 || data->n > 0) {
	matsize_length = 20;
      }
      out = realloc(out, length+value_length+4+matsize_length);
      if (!out) {
	return NULL;
      }
      // Write matrix dimensions if appropriate
      if (matsize_length > 0) {
	sprintf(out+length-1, "[%d][%d]", data->m, data->n);
	length += strlen(out+length-1);
      }
       
      out[length-1] = SPRINT_EQUALS;
      while (*val && *val <= ' ') {
	val++;
      }

#ifdef SPRINT_WRAP
      if (*val == '.' || (*val >= '0' && *val <= '9')) {
	/* Wrap numbers */
	out[length] = SPRINT_NUMBER_BEGIN;
	strcpy(out+length+1, data->value);
	length += (value_length+3);
	out[length-2] = SPRINT_NUMBER_END;
      }
      else {
	/* Wrap strings */
	out[length] = SPRINT_STRING_BEGIN;
	strcpy(out+length+1, data->value);
	length += (value_length+3);
	out[length-2] = SPRINT_STRING_END;
      }
#else
      /* Check whether value needs to be wrapped in quotes or
	 curlies */
      {
	int found_single_quote = 0;
	int found_double_quote = 0;
	int found_closing_curly = 0;
	int found_newline = 0;
	int found_separator = 0;
	char *c = data->value;
	while (*c != '\0') {
	  if (*c == '}') {
	    found_closing_curly = 1;
	  }
	  else if (*c == '\'') {
	    found_single_quote = 1;
	  }
	  else if (*c == '"') {
	    found_double_quote = 1;
	  }
	  else if (*c == '\n') {
	    found_newline = 1;
	  }
	  else if (*c <= ' ') {
	    found_separator = 1;
	  }
	  c++;
	}
#ifdef SPRINT_COMPRESS_WHITESPACE
	if (found_newline || found_separator) {
#else
	if (found_newline) {
#endif
	  /* Some kind of wrapper required */
	  if (!found_double_quote) {
	    out[length] = SPRINT_STRING_BEGIN;
#ifdef SPRINT_COMPRESS_WHITESPACE
	    __rc_copy_compress_whitespace(out+length+1, data->value);
	    length += strlen(out+length+1)+3;
#else
	    strcpy(out+length+1, data->value);
	    length += (value_length+3);
#endif
	    out[length-2] = SPRINT_STRING_END;
	  }
	  else if (!found_closing_curly) {
	    out[length] = SPRINT_NUMBER_BEGIN;
#ifdef SPRINT_COMPRESS_WHITESPACE
	    __rc_copy_compress_whitespace(out+length+1, data->value);
	    length += strlen(out+length+1)+3;
#else
	    strcpy(out+length+1, data->value);
	    length += (value_length+3);
#endif
	    out[length-2] = SPRINT_NUMBER_END;
	  }
	  else {
	    out[length] = SPRINT_NUMBER_BEGIN;
#ifdef SPRINT_COMPRESS_WHITESPACE
	    __rc_copy_compress_whitespace(out+length+1, data->value);
	    length += strlen(out+length+1)+3;
#else
	    strcpy(out+length+1, data->value);
	    length += (value_length+3);
#endif
	    out[length-2] = SPRINT_NUMBER_END;
	  }
	}
	else {
	  /* No wrapper required */
	  strcpy(out+length, data->value);
	  length += (value_length+1); 
	}
      }
#endif
    }

    /*    out[length-1] = '\n'; */
    --length;
    out[length] = '\0';
    
    if (!(data = data->next)) {
      break;
    }
  } /* Loop over param-value pairs */
  return out;
}

/* Return 1 if param exists in data, 0 otherwise. */
int
rc_exists(rc_data *data, const char *param)
{
  return (__rc_find(data, param) != NULL);  
}

/* Interpret the value associated with param as a boolean, returning 1
   if true and 0 if false. Note that 0 will be returned if param is
   not present, or if it exists and is "0", "no" or "false" (or any
   case insensitive variants). Any other scenario will result in 1
   being returned */
int
rc_get_boolean(rc_data *data, const char *param)
{
  data = __rc_find(data, param);
  if (!data) {
    return 0;
  }
  else if (!data->value) {
    return 1;
  }
  else if (strncasecmp(data->value, "false", 5) == 0
	   || strncasecmp(data->value, "no", 2) == 0) {
    return 0;
  }
  else {
    char *endptr = data->value;
    double val = strtod(data->value, &endptr);
    if (data->value == endptr || val != 0.0) {
      return 1;
    }
    else {
      return 0;
    }
  }
}

/* Interpret the value associated with param as an integer and return
   it. If successful *status is set to 1, but if either param is not
   present or the associated value cannot be interpreted as an
   integer, *status is set to 0. */
int
rc_get_int(rc_data *data, const char *param, int *status)
{
  data = __rc_find(data, param);
  if (!data || !data->value) {
    *status = 0;
    return 0;
  }
  else {
    char *endptr;
    long val = strtol(data->value, &endptr, 10);
    if (data->value == endptr) {
      *status = 0;
      return 0;
    }
    else {
      *status = 1;
      return val;
    }
  }
}

/* If param exists in data, set *value to the value associated with
   param and return 1, otherwise leave *value untouched and return
   0. */
int
rc_assign_int(rc_data *data, const char *param, int *value)
{
  int status;
  int val = rc_get_int(data, param, &status);
  if (status) {
    *value = val;
  }
  return status;
}

rc_real
rc_get_real(rc_data *data, const char *param, int *status)
{
  data = __rc_find(data, param);
  if (!data || !data->value) {
    *status = 0;
    return 0.0;
  }
  else {
    char *endptr;
    double val = strtod(data->value, &endptr);
    if (data->value == endptr) {
      *status = 0;
      return 0.0;
    }
    else {
      *status = 1;
      return val;
    }
  }
}

rc_real
rc_get_real_element(rc_data *data, const char *param, int index, int *status)
{
  char *str = rc_get_substring(data, param, index);
  if (!str) {
    *status = 0;
    return 0.0;
  }
  else {
    char *endptr;
    double val = strtod(str, &endptr);
    if (str == endptr) {
      *status = 0;
      free(str);
      return 0.0;
    }
    else {
      *status = 1;
      free(str);
      return val;
    }
  }
}

/* As rc_assign_int() but with a real */
int
rc_assign_real(rc_data *data, const char *param, rc_real *value)
{
  int status;
  rc_real val = rc_get_real(data, param, &status);
  if (status) {
    *value = val;
  }
  return status;
}

int
rc_assign_real_element(rc_data *data, const char *param, int index,
		       rc_real *value)
{
  int status;
  rc_real val = rc_get_real_element(data, param, index, &status);
  if (status) {
    *value = val;
  }
  return status;
}

/* If param exists in data and can be interpretted as 1 or more
   reals, assign *value to a vector of reals, returning the number
   found. If fewer than min_length are present, the vector will be
   padded with the last valid real found up to min_length. In this
   case the number returned is the number of valid reals found, not
   min_length. If no reals are found or if there is an error
   allocating memory, *value is untouched and 0 is returned. */
int
rc_assign_real_vector(rc_data *data, const char *param,
		      rc_real **value, int min_length)
{
  int length;
  rc_real *val = rc_get_real_vector(data, param, &length);
  if (!val) {
    return 0;
  }
  else if (length < min_length) {
    int i;
    val = realloc(val, min_length*sizeof(rc_real));
    if (!val) {
      return 0;
    }
    for (i = length; i < min_length; i++) {
      val[i] = val[length-1];
    }
  }
  *value = val;
  return length;
}

/* As rc_assign_real_vector() except that if no reals are found then
   a vector of length "min_length" is returned containing
   default_value. This function always returns min_length, except if
   there was an error allocating memory in which case 0 is
   returned. */
int
rc_assign_real_vector_default(rc_data *data, const char *param,
	      rc_real **value, int min_length, rc_real default_value)
{
  int length;
  rc_real *val = rc_get_real_vector(data, param, &length);
  if (!val) {
    int i;
    val = malloc(min_length*sizeof(rc_real));
    if (!val) {
      return 0;
    }
    for (i = 0; i < min_length; i++) {
      val[i] = default_value;
    }
  }
  else if (length < min_length) {
    int i;
    val = realloc(val, min_length*sizeof(rc_real));
    if (!val) {
      return 0;
    }
    for (i = length; i < min_length; i++) {
      val[i] = val[length-1];
    }
  }
  *value = val;
  return min_length;
}

/* Return the value associated with param as a string. NULL is
   returned on failure or memory allocation error. The string should
   be deallocated with rc_free().  */
char *
rc_get_string(rc_data *data, const char *param)
{
  data = __rc_find(data, param);
  if (!data || !data->value) {
    return NULL;
  }
  else {
    char *value = strdup(data->value);
    /* Remove trailing whitespace */
    __rc_remove_trailing_whitespace(value);
    return value;
  }
}

/* Interpret the string as a white-space separated list of strings,
   and return the ith string as a strdup-type copy */
char *
rc_substring(char *string, int i)
{
  int count = 0;
  char *c = string;
  char *start = c;
  char *end;

  while (*c) {
    /* Skip whitespace */
    while (*c && *c <= ' ') {
      c++;
    }
    if (*c == '\0') {
      return NULL;
    }
    
    /* Found the start of a value: find the end */
    if (*c == '{') {
      c++;
      start = c;
      /* Curly braces: find the terminating brace */
      while (*c && *c != '}') {
	c++;
      }
      end = c;
      if (*c) {
	c++;
      }
    }
    else if (*c == '"') {
      c++;
      start = c;
      /* Quotes: find the terminating quote */
      while (*c && *c != '"') {
	c++;
      }
      end = c;
      if (*c) {
	c++;
      }
    }
    else {
      /* Normal text */
      start = c;
      while (*c > ' ') {
	c++;
      }
      end = c;
    }
    
    if (count == i) {
      /* Found the required substring: copy and return */
      char terminator = *end;
      *end = '\0';
      char *substr = strdup(start);
      *end = terminator;
      return substr;
    }
    
    count++;
  }
  return NULL;
}

/* Interpret the string as a white-space separated list of substrings,
   and return the number of substrings */
int
rc_count_substrings(char *string)
{
  int count = 0;
  char *c = string;

  while (*c) {
    /* Skip whitespace */
    while (*c && *c <= ' ') {
      c++;
    }
    if (*c == '\0') {
      return 0;
    }
    
    /* Found the start of a value: find the end */
    if (*c == '{') {
      c++;
      /* Curly braces: find the terminating brace */
      while (*c && *c != '}') {
	c++;
      }
      if (*c) {
	c++;
      }
    }
    else if (*c == '"') {
      c++;
      /* Quotes: find the terminating quote */
      while (*c && *c != '"') {
	c++;
      }
      if (*c) {
	c++;
      }
    }
    else {
      /* Normal text */
      while (*c > ' ') {
	c++;
      }
    }
    count++;
  }
  return count;
}

/* Interpret the "value" corresponding to "param" as a white-space
   separated list of strings, and return the ith string */
char *
rc_get_substring(rc_data *data, const char *param, int i)
{
  data = __rc_find(data, param);
  if (!data || !data->value) {
    return NULL;
  }
  else {
    char *substr = rc_substring(data->value, i);
    return substr;
  } 
}

/* Interpret the "value" corresponding to "param" as a white-space
   separated list of substrings, and return the number of
   substrings */
int
rc_num_substrings(rc_data *data, const char *param)
{
  data = __rc_find(data, param);
  if (!data || !data->value) {
    return 0;
  }
  else {
    return rc_count_substrings(data->value);
  } 
}

/* Interpret the value associated with param as a string of
   whitespace-separated substrings, and return the number of
   substrings.  Additionally, if the size of the variable has been
   specified with square brackets [m,n] then return the values m and
   n. */
int
rc_size(rc_data *data, const char *param, int *m, int *n)
{
  data = __rc_find(data, param);
  if (!data || !data->value) {
    *m = *n = 0;
    return 0;
  }
  else {
    *m = data->m;
    *n = data->n;
    return rc_count_substrings(data->value);
  }
}

/* Free dynamically allocated data returned by rc_get_string(),
   rc_assign_string(), rc_assign_real_vector(),
   rc_assign_real_vector_default(), rc_get_int_vector() and
   rc_get_real_vector(). */
void
rc_free(void *string)
{
  if (string)
    free(string);
}

/* Free dynamically allocated matrix returned by
   rc_get_real_matrix() */
void
rc_free_matrix(void **matrix)
{
  if (matrix) {
    if (matrix[0]) {
      free(matrix[0]);
    }
    free(matrix);
  }
}


/* If param exists in data, set *value to a pointer to a copy of the
   string associated with param and return 1, otherwise leave value
   untouched and return 0. Memory allocation error also results in 0
   being returned. The string should be freed with rc_free().  */
int
rc_assign_string(rc_data *data, const char *param, char **value)
{
  char *str = rc_get_string(data, param);
  if (str) {
    if (*str) {
      /* A genuine non-empty string */
      *value = str;
      return 1;
    }
    else {
      /* An empty string */
      free(str);
    }
  }
  return 0;
}

/* If param exists, interpret the associated value as integers and
   return a pointer to an int vector or NULL on failure. *length
   contains the number of integers assigned, 0 if none or memory
   allocation error. The vector should be freed with rc_free(). */
int *
rc_get_int_vector(rc_data *data, const char *param, int *length)
{
  int *out = NULL;
  data = __rc_find(data, param);
  *length = 0;

  if (!data || !data->value) {
    return NULL;
  }
  else {
    char *endptr;
    char *c = data->value;
    while (*c) {
      long val = strtol(c, &endptr, 10);
      if (endptr != c) {
	out = realloc(out, (++*length)*sizeof(int));
	if (!out) {
	  *length = 0;
	  return NULL;
	}
	out[*length-1] = val;
	c = endptr;
      }
      else {
	break;
      }
    }
    return out;
  }
}

/* As rc_get_real_vector() but with reals. */
rc_real *
rc_get_real_vector(rc_data *data, const char *param, int *length)
{
  rc_real *out = NULL;
  data = __rc_find(data, param);
  *length = 0;

  if (!data || !data->value) {
    return NULL;
  }
  else {
    char *endptr;
    char *c = data->value;
    while (*c) {
      double val = strtod(c, &endptr);
      if (endptr != c) {
	out = realloc(out, (++*length)*sizeof(rc_real));
	if (!out) {
	  *length = 0;
	  return NULL;
	}
	out[*length-1] = val;
	c = endptr;
      }
      else {
	break;
      }
    }
    return out;
  }
}

/* Return 1 if param exists and is a 2D matrix, 0 otherwise */
int
rc_is_matrix(rc_data *data, const char *param)
{
  data = __rc_find(data, param);
  if (data && data->m > 0 && data->n > 0) {
    return 1;
  }
  else {
    return 0;
  }
}

/* If param exists, interpret the associated value as reals in a 2D
   matrix and return a pointer to an array of pointers to each row of
   the matrix, or NULL on failure. *m and *n contain the number of
   rows and columns, respectively, 0 if none or memory allocation
   error. The matrix should be freed with rc_free_matrix(). */
rc_real **
rc_get_real_matrix(rc_data *data, const char *param, int *m, int *n)
{
  int length, i;
  rc_real *val, **pval;

  val = rc_get_real_vector(data, param, &length);
  if (!val) {
    return NULL;
  }

  data = __rc_find(data, param);
  if (!data) {
    rc_free(val);
    return NULL;
  }

  if (data->m * data->n != length) {
    rc_free(val);
    fprintf(stderr, "Error reading \"%s\": should have %d*%d=%d elements but found %d\n",
	    param, data->m, data->n, data->m*data->n, length);
    return NULL;
  }

  pval = malloc(length*sizeof(rc_real*));
  if (!pval) {
    rc_free(val);
    return NULL;
  }
  for (i = 0; i < data->m; i++) {
    pval[i] = val + i*data->n;
  }
  *m = data->m;
  *n = data->n;
  return pval;
}


/* As rc_get_real_matrix, but with ints */
int **
rc_get_int_matrix(rc_data *data, const char *param, int *m, int *n)
{
  int length, i;
  int *val, **pval;

  val = rc_get_int_vector(data, param, &length);
  if (!val) {
    return NULL;
  }

  data = __rc_find(data, param);
  if (!data) {
    rc_free(val);
    return NULL;
  }

  if (data->m * data->n != length) {
    rc_free(val);
    fprintf(stderr, "Error reading \"%s\": should have %d*%d=%d elements but found %d\n",
	    param, data->m, data->n, data->m*data->n, length);
    return NULL;
  }

  pval = malloc(length*sizeof(int*));
  if (!pval) {
    rc_free(val);
    return NULL;
  }
  for (i = 0; i < data->m; i++) {
    pval[i] = val + i*data->n;
  }
  *m = data->m;
  *n = data->n;
  return pval;
}



int
rc_set_section(rc_data *data, const char *section)
{
  if (data->section_reqd) {
    free(data->section_reqd);
    data->section_reqd = NULL;
  }
  if (section != NULL) {
    data->section_reqd = strdup(section);
  }
  return 1;
}

int
rc_append_section(rc_data *data, const char *subsection)
{
  if (data->section_reqd) {
    int section_len = strlen(data->section_reqd);
    int subsection_len = strlen(subsection);
    data->section_reqd = realloc(data->section_reqd,
				    subsection_len + subsection_len + 2);
    if (!data->section_reqd) {
      return 0;
    }
    data->section_reqd[section_len] = SECTION_CHAR;
    strcpy(data->section_reqd+section_len+1, subsection);
  }
  else {
    data->section_reqd = strdup(subsection);
    if (!data->section_reqd) {
      return 0;
    }
  }
  return 1;
}

int
rc_pop_section(rc_data *data)
{
  if (!data->section_reqd) {
    /* No section to "pop"! */
    return 0;
  }
  else {
    int section_len = strlen(data->section_reqd);
    char *c = data->section_reqd + section_len - 1;
    while (c != data->section_reqd && *c != SECTION_CHAR) {
      c--;
    }
    if (*c == SECTION_CHAR) {
      *c = '\0';
    }
    else {
      free(data->section_reqd);
      data->section_reqd = NULL;
    }
    return 1;
  }
}

rc_data*
rc_find_next_member(rc_data *data, const char *param)
{
  int len = strlen(param);
  while (data) {
    if (data->param
	&& strncasecmp(data->param, param, len) == 0
	&& data->param[len] == '.') {
      break;
    }
    data = data->next;
  }
  return data;
}
 
rc_data*
rc_find_next_simple_member(rc_data *data, const char *param)
{
  data = rc_find_next_member(data, param);
  if (data
      && data->param
      && !strchr(data->param+strlen(param)+1, '.')) {
    return data;
  }
  else {
    if (data && data->next) {
      return rc_find_next_simple_member(data->next, param);
    }
    else {
      return NULL;
    }
  }
} 
