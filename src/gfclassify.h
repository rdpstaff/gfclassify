/*
 *    Copyright (C) 2012 Michigan State University
 *
 *    This file is part of GFClassify.
 *
 *    GFClassify is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    GFClassify is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with GFClassify.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _SCORE_HH
#define _SCORE_HH

#include  "icm.hh"

typedef struct score_opts_t {
  char** model_paths;
  char* infile_name;
  int num_models;
  bool use_bg;
  bool check_reverse;
  double threshold;
} score_opts;

static void Parse_Command_Line(int argc, char * argv [], score_opts* opts);
static int Read_String(FILE * fp, char * & s, long int & s_size, char * & tag, long int & tag_size);
static void Usage(char * command);
static char* get_model_name(const char* );
double score(char* seq, int len, ICM_t* models, int num_models, double* scores);
void rev_comp(char* seq, int len) ;
char comp(char c);


#endif
