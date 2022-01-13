/* writeconfig.c -- Functions to write data to a configuration file

   Copyright (C) 2006 The University of Reading
   Author: Robin Hogan <r.j.hogan@reading.ac.uk>

   This software is licensed under the terms of the Apache Licence
   Version 2.0 which can be obtained at
   http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "readconfig.h"

int
rc_write_comment(FILE *file, const char *comment)
{
  return fprintf(file, "# %s\n", comment);
}

int
rc_write_string(FILE *file, const char *param, const char *value)
{
  return fprintf(file, "%s %s\n", param, value);
}

int
rc_write_ints(FILE *file, const char *param, int n, const int *values)
{
  int i;
  fprintf(file, "%s", param);
  for (i = 0; i < n; i++) {
    fprintf(file, " %d", values[i]);
  }
  return fprintf(file, "\n");
}

int
rc_write_reals(FILE *file, const char *param, int n, const rc_real *values)
{
  int i;
  fprintf(file, "%s", param);
  for (i = 0; i < n; i++) {
    fprintf(file, " %g", values[i]);
  }
  return fprintf(file, "\n");
}

int
rc_write_matrix(FILE *file, const char *param, int m, int n,
		const rc_real **values)
{
  int i, j;
  fprintf(file, "%s[%d][%d] {\n", param, m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      fprintf(file, " %g", values[i][j]);
    }
    fprintf(file, "\n");
  }
  return fprintf(file, "}\n");  
}

