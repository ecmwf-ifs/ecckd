/* readconfig.h -- Functions to read configuration information from a file

   Copyright (C) 2003-2005 The University of Reading
   Author: Robin Hogan <r.j.hogan@reading.ac.uk>

   This software is licensed under the terms of the Apache Licence
   Version 2.0 which can be obtained at
   http://www.apache.org/licenses/LICENSE-2.0.
*/

#ifndef _READCONFIG_H
#define _READCONFIG_H 1

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

typedef double rc_real;


/* Configuration information is stored as a singly linked list as a
   set of param-value pairs. */
typedef struct __rc_data rc_data;
struct __rc_data {
  char *param;
  char *value;
  char *section_reqd;
  int m, n;
  rc_data *next;
};


/* READING FUNCTIONS */

/* Read configuration information from file called file_name and
   return a pointer to the rc_data structure, or NULL if an error
   occurred.  If err_file is not NULL, errors messages will be written
   to err_file. */
rc_data *rc_read(const char *file_name, FILE *err_file);

/* Append contents of file to an existing data structure, returning 1
   on success or 0 on failure */
int rc_append(rc_data *data, const char *file_name, FILE *err_file);

/* Free an rc_data structure created with rc_read() - rc_data is a
   singly-linked list and its contents are recursively detelted. */
void rc_clear(rc_data *data);

/* Add a param-value pair to an existing rc_data structure,
   overwriting an existing param with the same name (case
   insensitive). Returns 1 on success and 0 on failure. */
int rc_register(rc_data *data, const char *param, const char *value);

/* Search the command-line arguments for files: ones that are *not*
   param=value pairs or -param arguments, and add them to the rc_data
   structure as FILE0, FILE1 etc. */
int rc_register_files(rc_data *data, int argc, const char **argv);

/* Search the command-line arguments for param=value pairs and -param
   arguments and add them to the rc_data structure. Returns 1 on
   success and 0 on failure. */
int rc_register_args(rc_data *data, int argc, const char **argv);

/* Find the first config file on the command line, returning an index
   to the argv array, and 0 if none is found. Basically the first
   argument that contains no "=" sign and doesn't start with "-" is
   returned. */
int rc_get_file(int argc, const char **argv);

/* Print contents of rc_data structure to file. */
void rc_print(rc_data *data, FILE *file);

/* Return the contents of an rc_data structure as a string. Returns
   NULL on failure. Free with rc_free(). */
char *rc_sprint(rc_data *data);

/* Return 1 if param exists in data, 0 otherwise. */
int rc_exists(rc_data *data, const char *param);

/* Interpret the value associated with param as a boolean, returning 1
   if true and 0 if false. Note that 0 will be returned if param is
   not present, or if it exists and is "0", "no" or "false" (or any
   case insensitive variants). Any other scenario will result in 1
   being returned */
int rc_get_boolean(rc_data *data, const char *param);

/* Interpret the value associated with param as an integer and return
   it. If successful *status is set to 1, but if either param is not
   present or the associated value cannot be interpreted as an
   integer, *status is set to 0. */
int rc_get_int(rc_data *data, const char *param, int *status);

/* Interpret the value associated with param as a real. */
rc_real rc_get_real(rc_data *data, const char *param, int *status);

/* Return the value associated with param as a string. NULL is
   returned on failure or memory allocation error. The string should
   be deallocated with rc_free().  */
char *rc_get_string(rc_data *data, const char *param);

/* Interpret the value associated with param as a string of
   whitespace-separated substrings, and return the ith one. The
   string should be deallocated with rc_free(). */
char *rc_get_substring(rc_data *data, const char *param, int i);

/* Interpret the value associated with param as a string of
   whitespace-separated substrings, and return the number of
   substrings, or 0 if the value is either all whitespace or if param
   does not exist. */
int rc_num_substrings(rc_data *data, const char *param);

/* Interpret the value associated with param as a string of
   whitespace-separated substrings, and return the number of
   substrings.  Additionally, if the size of the variable has been
   specified with square brackets [m,n] then return the values m and
   n. */
int rc_size(rc_data *data, const char *param, int *m, int *n);

/* As rc_get_substring, but extracting the substring from a
   specified string. */
char *rc_substring(char *string, int i);

/* If param exists, interpret the associated value as integers and
   return a pointer to an int vector or NULL on failure. *length
   contains the number of integers assigned, 0 if none or memory
   allocation error. The vector should be freed with rc_free(). */
int *rc_get_int_vector(rc_data *data, const char *param, int *length);

/* As rc_get_int_vector() but with reals. */
rc_real *rc_get_real_vector(rc_data *data, const char *param, int *length);

/* Get the specified element of a real vector; otherwise as rc_get_real() */
rc_real rc_get_real_element(rc_data *data, const char *param, int index,
			      int *status);

/* Return 1 if param exists and is a 2D matrix, 0 otherwise */
int rc_is_matrix(rc_data *data, const char *param);

/* If param exists, interpret the associated value as reals in a 2D
   matrix and return a pointer to an array of pointers to each row of
   the matrix, or NULL on failure. *m and *n contain the number of
   rows and columns, respectively, 0 if none or memory allocation
   error. The matrix should be freed with rc_free_matrix(). */
rc_real **rc_get_real_matrix(rc_data *data, const char *param, int *m, int *n);

  /* As rc_get_real_matrix but with ints. */
int **rc_get_int_matrix(rc_data *data, const char *param, int *m, int *n);

/* Free dynamically allocated data returned by rc_get_string(),
   rc_assign_string(), rc_assign_real_vector(),
   rc_assign_real_vector_default(), rc_get_int_vector() and
   rc_get_real_vector(). */
void rc_free(void *string);

/* Free dynamically allocated matrix returned by
   rc_get_real_matrix() */
void rc_free_matrix(void **matrix);

/* If param exists in data, set *value to the value associated with
   param and return 1, otherwise leave *value untouched and return
   0. */
int rc_assign_int(rc_data *data, const char *param, int *value);

/* As rc_assign_int() but with a real */
int rc_assign_real(rc_data *data, const char *param, rc_real *value);

int rc_assign_real_element(rc_data *data, const char *param, int index, rc_real *value);

/* If param exists in data, set *value to a pointer to a copy of the
   string associated with param and return 1, otherwise leave value
   untouched and return 0. Memory allocation error also results in 0
   being returned. The string should be freed with rc_free().  */
int rc_assign_string(rc_data *data, const char *param, char **value);

/* If param exists in data and can be interpretted as 1 or more
   reals, assign *value to a vector of reals, returning the number
   found. If fewer than min_length are present, the vector will be
   padded with the last valid real found up to min_length. In this
   case the number returned is the number of valid reals found, not
   min_length. If no reals are found or if there is an error
   allocating memory, *value is untouched and 0 is returned. */
int rc_assign_real_vector(rc_data *data, const char *param,
			  rc_real **value, int min_length);

/* As rc_assign_real_vector() except that if no reals are found then
   a vector of length "min_length" is returned containing
   default_value. This function always returns min_length, except if
   there was an error allocating memory in which case 0 is
   returned. */
int rc_assign_real_vector_default(rc_data *data, const char *param,
	  rc_real **value, int min_length, rc_real default_value);

int rc_set_section(rc_data *data, const char *section);
int rc_append_section(rc_data *data, const char *section);
int rc_pop_section(rc_data *data);

/* Return pointer to next parameter-value pair where parameter is of
   the form param.member, or NULL if no member is found. The "simple"
   version does not return members that themselves have sub-members. */
rc_data* rc_find_next_member(rc_data *data, const char *param); 
rc_data* rc_find_next_simple_member(rc_data *data, const char *param); 


/* WRITING FUNCTIONS */

int rc_write_comment(FILE *file, const char *comment);
int rc_write_string(FILE *file, const char *param,
		    const char *value);
int rc_write_ints(FILE *file, const char *param, int n,
		  const int *values);
int rc_write_reals(FILE *file, const char *param, int n,
		   const rc_real *values);
int rc_write_matrix(FILE *file, const char *param, int m, int n,
		    const rc_real **values);

#ifdef __cplusplus
}
#endif

#endif
