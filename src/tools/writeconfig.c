/* writeconfig.c -- Functions to write data to a configuration file

   Copyright (C) 2006 Robin Hogan <r.j.hogan@reading.ac.uk>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

