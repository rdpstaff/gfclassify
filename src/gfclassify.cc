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

// GFClassify core program, heavily modified version of
// score-fixed.cc from Glimmer 3.02 (command line parser
// and fasta parser are the only code blocks present
// in the original, and have been modified)

#include "gfclassify.h"
#include <stdio.h>
#include <libgen.h>
#include <float.h>

//**ALD  Gets rid of make undefined reference error
int Unused = Filter ('a');

//Report a fatal error to stderr and
//exit with a non-zero exit code
void fatal_error(const char* msg, ...) {
  va_list argp;
  va_start(argp, msg);
  vfprintf(stderr, msg, argp);
  va_end(argp);
  exit(-1);
}

int main(int argc, char * argv []) {
  //First parse the command line arguments
  score_opts opts;
  Parse_Command_Line(argc, argv, &opts);

  //Variable declorations and allocates 
  char* string = NULL, * tag = NULL;       //Current sequence storage
  long int  string_size = 0, tag_size = 0; //buffer sizes
  int  string_num = 0;                     //Sequence count

  //Model allocations
  ICM_t* models = (ICM_t*)SAFE_MALLOC(sizeof(ICM_t) * opts.num_models);
  char** model_names = (char**)SAFE_MALLOC(sizeof(char*) * opts.num_models);
  double* scores_f = (double*)SAFE_MALLOC(sizeof(double) * opts.num_models);
  double* scores_r = (double*)SAFE_MALLOC(sizeof(double) * opts.num_models);

  //Write the header out
  printf("#seqid\tdescription\torientation");
  //Read the models in and finish writing the header
  for(int index = 0;index < opts.num_models;index++) {
    models[index].Read(opts.model_paths[index]); //Read the model
    model_names[index] =  get_model_name(opts.model_paths[index]); //Store the name
    printf("\t%s log odds", model_names[index]);  //Write the model column header
  }  
  printf("\tlabel\n");

  //Figure out (and open if required) the input source
  FILE* in = NULL;
  if(strcmp(opts.infile_name, "-") == 0) {
    in = stdin;
  } else {
    in = fopen(opts.infile_name, "r");
    if(in == NULL) {
      fatal_error("Couldn't open input file %s\n", opts.infile_name);
    }
  }

  //Read sequences from the input (the last four are ouput variables)
  while (Read_String (in, string, string_size, tag, tag_size)) {
    int len;
    char dir = '+'; //Default direction is +
    
    string_num ++;
    len = strlen (string); //How long is our sequence?
    double* scores;

    //Figure out the header (if there is one)
    char* header = NULL;
    for (int i = 0; i < tag_size;i++) {
      if(tag[i] == ' ') {
	//Pointer play, split the one string in to two by putting a \0 in it and keeping a pointer to the next part
	header = tag + i + 1; 
	tag[i] = '\0';
	break;
      }
    }

    //Score the forward sequence against the model (storing the scores in scores_f), returns the highest score
    double f_max = score(string, len, models, opts.num_models, scores_f);
    scores = scores_f;

    if(opts.check_reverse) {
      //Do the reverse checking
      rev_comp(string, len);
      
      double r_max = score(string, len, models, opts.num_models, scores_r);

      //Assume the best direction is the one that gives
      //the single best log odds (not highest log likelihood)
      if(r_max > f_max) {
	dir = '-';
	scores = scores_r;
      }
    }

    //Print out the scores for the sequence
    printf("%s\t%s\t%c", tag, header, dir);
    for(int index = 0;index < opts.num_models;index++) {
      printf("\t%10.4f", scores[index]);
    }
    
    //If we've got more than one model we need to figure out a label
    if(opts.num_models > 1) {
      int c1 = 0, c2 = 1;
      //So we need to figure out the top two models to
      //find the likelihood
      for(int index = 1;index < opts.num_models;index++) {
	if(scores[index] > scores[c1]) {
	  c2 = c1;
	  c1 = index;
	} 
      }
      
      int c;
      //The likelihood ratio can be computed in
      //log space by subtraction, test the log
      //likelihood to see if it is above the reporting
      //threshold, class -1 means rejected
      if(scores[c1] - scores[c2] > opts.threshold) {
	c = c1;
      } else {
	c = -1;
      }

      //Find the actual string label from the index
      char* label = NULL;
      if(c > -1) {
	label = model_names[c];
      } else {
	label = "rejected";
      }

      printf("\t%s", label);
    }
    printf("\n");
  }

  //Release all the memory
  free(scores_f);
  free(scores_r);
  free(models);
  for(int index = 0;index < opts.num_models;index++) {
    free(model_names[index]);
  }
  free(model_names);
  
  return 0;
}

//Score a sequence (seq, with length len) against the num_models input models
//store the resulting log odds in scores (assumed allocated with enough space
//to hold num_models scores)
double score(char* seq, int len, ICM_t* models, int num_models, double* scores) {
  double max_score = -DBL_MAX;

  for(int index = 0;index < num_models;index++) {
    int m_len = models[index].Get_Periodicity(); //Be sure we're supplying the correct frame
    //Run the actual Glimmer ICM scoring function
    scores[index] = models[index].Score_String(seq, len, (index % m_len));
    if(scores[index] > max_score) {
      max_score = scores[index];
    }
  }

  return max_score;
}

//Helper function to strip the path and extension
//from an input file name
static char* get_model_name(const char* model_path) {
  char* ret = (char*)SAFE_MALLOC(sizeof(char) * strlen(model_path));
  strcpy(ret, model_path);
  ret = basename(ret);

  char* tmp = (char*)SAFE_MALLOC(sizeof(char) * strlen(ret));
  strcpy(tmp, ret);
  //  free(ret);
  ret = tmp;

  for(int index = strlen(ret) - 1;index >= 0;index--) {
    if(ret[index] == '.') {
      ret[index] = '\0';
      break;
    }
  }

  return ret;
}

static void Parse_Command_Line(int argc, char * argv [], score_opts* opts) {

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .
  bool  errflg = false;
  int  ch;

  opts->threshold = 0;
  opts->use_bg = false;
  opts->check_reverse = false;

  optarg = NULL;

  while  (! errflg
	  && ((ch = getopt (argc, argv, "ht:bc")) != EOF)) {
    switch  (ch) 
      {
      case 'h':
	errflg = true;
	break;
      case 't':
	opts->threshold = strtof(optarg, NULL);
	break;
      case 'b':
	opts->use_bg = true;
	break;
      case 'c':
	opts->check_reverse = true;
	break;
      case '?':
	fprintf (stderr, "Unrecognized option -%c\n", optopt);
      default:
	errflg = TRUE;
      }
  }

  if (optind > argc - 2) {
    fprintf(stderr, "You must supply a sequence file and at least one icm to score queries against\n");
    errflg = true;
  }

  if (errflg) {
    Usage(argv [0]);
    exit(EXIT_FAILURE);
  }

  opts->infile_name = argv[optind++];
  opts->num_models = argc - optind;
  opts->model_paths = argv + optind;

  return;
}

void rev_comp(char* seq, int len) {
  int start = 0, end = len;

  for(;start < end;start++, end--) {
    char tmp = seq[start];
    seq[start] = comp(seq[end]);
    seq[end] = comp(tmp);
  }

  if(start == end) {
    seq[start] = comp(seq[start]);
  }
}

char comp(char c) {
  switch(c) {
  case 'A': case 'a':
    return 't';
  case 'G': case 'g':
    return 'c';
  case 'C': case 'c':
    return 'g';
  case 'T': case 't': case 'U': case 'u':
    return 'a';
  default:
    return 'n';
  }
}


int  Read_String
    (FILE * fp, char * & s, long int & s_size, char * & tag,
     long int & tag_size)

//  Read next string from  fp  (assuming FASTA format) into  s [0 .. ]
//  which has  s_size  characters.  Allocate extra memory if needed
//  and adjust  s_size  accordingly.  Return  TRUE  if successful,  FALSE
//  otherwise (e.g., EOF).  Put FASTA header line into  tag [0 .. ]
//  (and adjust  tag_size  if needed).

{
  int  ch, ct;

  while  ((ch = fgetc (fp)) != EOF && ch != '>')
    ;

  if  (ch == EOF)
    return  FALSE;

  ct = 0;
  while  ((ch = fgetc (fp)) != EOF && ch != '\n' && isspace (ch))
    ;

  if  (ch == EOF)
    return  FALSE;
  if  (ch != '\n' && ! isspace (ch))
    ungetc (ch, fp);
  while  ((ch = fgetc (fp)) != EOF && ch != '\n') {
    if  (ct >= tag_size - 1) {
      tag_size += INCR_SIZE;
      tag = (char *) Safe_realloc (tag, tag_size);
    }
    if(isspace(ch)) {
      ch = int(' ');
    }
    tag [ct ++] = char (ch);
  }
  tag [ct ++] = '\0';
  
  ct = 0;
  while  ((ch = fgetc (fp)) != EOF && ch != '>') {
    if  (isspace (ch))
      continue;

    ch = tolower(ch);
    if(ch != 'a' && ch != 'c' && ch != 'g' && ch != 't' && ch != 'u' && ch != 'n') {
      fatal_error("Unable to parse sequence %s\n", tag);
    }

    if  (ct >= s_size - 1)
      {
	s_size += INCR_SIZE;
	s = (char *) Safe_realloc (s, s_size);
      }
    s [ct ++] = char (ch);
  }
  s [ct ++] = '\0';
  
  if  (ch == '>')
    ungetc (ch, fp);
  
  return  TRUE;
}



static void Usage(char * command) {
//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  fprintf (stderr,
           "USAGE:  %s [options] <seq_file> <model1> [model2 ...] \n"
           "\n"
           "Read sequences from  seq_file (- for stdin) and score them using\n"
           "model(s) in file(s)  <model1> [model2 ...] .  Output scores to  stdout.\n"
           "\n"
           "Options:\n"
           " -h        Print this message\n"
	   " -t        Threshold for classification\n"
           " -c        Check the reverse complement of the sequences\n"
           "\n",
           command);

  return;
}



