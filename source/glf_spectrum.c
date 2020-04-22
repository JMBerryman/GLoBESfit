/* GLoBESfit -- GLoBES fitting tools
*  (C) 2019-2020 The GLoBESfit Team
*
* GLoBESfit is mainly intended for academic purposes. Proper
* credit must be given if you use GLoBESfit or parts of it. Please
* read the section 'Credit' in the README file.
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* (C) 2005, 2007, 2019, 2020 Patrick Huber, J.M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <argp.h>

#include <globes/globes.h>   /* GLoBES library */

#include "glf_types.h"
#include "glf_spectrum_chi.h"
#include "glf_precomputed_probabilities.h"
#include "glf_probability.h"

#define VERSION "1.0"
#define MAXNUM 12 /* largest possible number of oscillation parameters */
#define YES 1
#define NO 0

double glf_systematic[10];

const char *argp_program_version =
"GLoBESfit_spectra "VERSION"\n(C) 2007, 2019, 2020 J. M. Berryman, P. Huber \n"
"This is free software see the source for copying conditions. There is NO\n"
"warranty; not even for MERCHANTABILITY or"
" FITNESS FOR A PARTICULAR PURPOSE.";
const char *argp_program_bug_address = "globesfit-g@vt.edu";

/* Program documentation. */
static char doc[] ="Data fitting for SBL experiments with GLoBES";

/* A description of the arguments we accept. */
static char args_doc[] = "(No Inputs Allowed!)";

/* The options we understand. */
static struct argp_option options[] ={
  {"resolution",'r',  "NUMBER", 0,  
   "Number of points to be used in each direction",3},
  {"output",'o', "FILE", 0,
   "Output to FILE instead of standard output",3},
  {"parameters",'p',"PARAMETERS",0,
   "Set of oscillation parameters\n   for which the data are computed\n"
   "		'th12,th13,th23,th24,dm41,Eres,th14,th34,dm21,dm31'",0},
  {"xrange",'x',"RANGE",0,"Range in log(sin^2 2theta13), Format 'a,b'",3},
  {"yrange",'y',"RANGE",0,"Range in log(Delta m_{41}^2/eV^2), Format 'a,b'",3},
  {"marginalize",'m',"DIRECTION",0,"Marginalize with respect to parameter number DIRECTION",0},
  {"Rates",'R',0,0,"Show Raw Event Rates",0},
  {"Project",'P',"DIRECTION",0,"Project onto DIRECTION",0}, 
  {"chi",'c',0,0,"Compute chi^2 for parameters given by -p",0},
  {"Initial-step",'I',"STEPSIZE",0,"Initial step size for minimizer",0},
  {"block",'b',"BLOCK", 0, "Experiments that are **NOT** \n to be included in present analysis; \n see documentation for definitions",0},
  {"Theta13", 'T', 0, 0, "Execute a scan over theta_13 as well as sterile parameters"},
  {"Calibrate", 'C', 0, 0, "Execute a scan over theta_13 and Dm_{31}^2"},
  { 0 } 
};

/* Used by `main' to communicate with `parse_opt'. */
struct arguments
{
/*  int exp; */ /* number of experiments */
  char* args[32];                /* up to 32 experiments */
  int resolution,rate;
  int marg[MAXNUM]; 
  int parg[MAXNUM];
  int barg[16]; /* Array for blocks of experiments being disregarded */
  double step;
  int margf,pargf,bargf;
  char *output_file;
  char *params,*xrange,*yrange;
  int chi, doth13, calibrate;
};

static void parse_definition(const char *in)
{
  const char *delim="=";
  double val;
  char *token=NULL;
  char *lhs=NULL;
  char *rhs=NULL;
  char *inc=NULL;
  size_t length=0;

  inc=strdup(in);
  if(inc==NULL) return;
  token=strtok(inc,delim);
  if(token!=NULL) 
    {
     lhs=strdup(token);
     length++;
    }
  token=strtok(NULL,delim); 
  if(token!=NULL) 
    {
     rhs=strdup(token);
     length++;
    }
  
  if(length!=2) {
    free(lhs);
    free(rhs);
    fprintf(stderr,"ERROR: Definition is not of form 'DEFINITION=VALUE'\n");
    return;}
  val=atof(rhs);
  glbDefineAEDLVariable(lhs,val);
  free(lhs);
  free(rhs);
  return;
}

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the INPUT argument from `argp_parse', which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;
  
  int witch;

  switch (key)
    {
    case 'r':
      arguments->resolution = atoi(arg);
      break; 

    case 'R':
      arguments->rate = YES;
      break; 

    case 'T':
      arguments->doth13 = YES;
      break;

    case 'C':
      arguments->calibrate = YES;
      break;

    case 'c':
      arguments->chi = YES;
      break;
    
    case 'I':
      arguments->step = atof(arg);
      break;
      
    case 'o':
      arguments->output_file = arg;
      break;

    case 'm':
      witch=atoi(arg);
      if(witch<MAXNUM&&witch>-1)
	{
	  arguments->marg[witch]=GLB_FREE;
	  arguments->margf=YES;
	}
      break;

    case 'P':
      witch=atoi(arg);
      if(witch<MAXNUM&&witch>-1)
	{
	  arguments->parg[witch]=GLB_FIXED;
	  arguments->pargf=YES;
	}
      break;

    case 'b':
      witch=atoi(arg);
      if(witch<16&&witch>-1)
	{
	  arguments->barg[witch]=0;
	  arguments->bargf=YES;
	}
      break;

    case 'p':
      arguments->params = arg;
      break;

    case 'x':
      arguments->xrange = arg;
      break;

    case 'y':
      arguments->yrange = arg;
      break;
/*
  This is where the arguments are limited; keep an eye out here...
  Commenting out everything that's looking to load in experiments

    case ARGP_KEY_ARG:
      if (state->arg_num > 32) */
	/* Too many arguments. */
/*	argp_usage (state);
     
      arguments->args[state->arg_num]=arg;
      
      break;
      
    case ARGP_KEY_END:
      if (state->arg_num < 1) */
	/* Not enough arguments. */
/*	argp_usage (state);
      break;
      
    default:
      return ARGP_ERR_UNKNOWN;
*/    }
/*  arguments->exp=state->arg_num; */
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

int main(int argc, char *argv[])
{
  double progress;
  int count=0;
  int s,s2,s3,s4,s5,s6,s7,s8;
  int i,rv;
  FILE *output=NULL;
  int resolution;

  int NUMP=10;

  double osc[MAXNUM];
  double xrange[]={-4,0};
  double yrange[]={-2,2};
  struct arguments arguments;

  arguments.args[0]="-";
  arguments.resolution=51;
  arguments.rate=NO;
  arguments.output_file=NULL;

  arguments.params=NULL;
  arguments.xrange=NULL;
  arguments.yrange=NULL;
/*  arguments.exp=0; */

  for(i=0;i<MAXNUM;i++)  arguments.marg[i]=GLB_FIXED;
  for(i=0;i<3;i++)  arguments.parg[i]=GLB_FIXED;
  for(i=3;i<MAXNUM;i++)  arguments.parg[i]=GLB_FREE;
  for(i=0;i<16;i++)  arguments.barg[i]=1;

  arguments.margf=NO;
  arguments.pargf=NO;
  arguments.bargf=NO;

  arguments.chi=NO;
  arguments.doth13=NO;
  arguments.calibrate=NO;
  arguments.step=0.05;
  /* parsing the command line */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);  
  
  

  if(arguments.yrange!=NULL){
    rv=sscanf(arguments.yrange,"%lf , %lf",&yrange[0],&yrange[1]);
    if(rv!=2){
      fprintf(stderr,"%s: FATAL: Wrong format range\n ",argv[0]);
      exit(1);
      }
    }
 
  if(arguments.xrange!=NULL){
      rv=sscanf(arguments.xrange,"%lf , %lf",&xrange[0],&xrange[1]);
      if(rv!=2){
	  fprintf(stderr,"%s: FATAL: Wrong format range\n ",argv[0]);
	  exit(1);
	}
    } 

  if(arguments.output_file==NULL)  output=stdout;
  else {
    output=fopen(arguments.output_file,"w");
    if(output==NULL) {
      fprintf(stderr,"FATAL: %s: Could not open output file '%s'\n",argv[0],
	      arguments.output_file);
      exit(1);
    }
  }

  resolution=arguments.resolution;

  /* Initialize libglobes */
  glbInit(argv[0]);

  glbSetInitialStep(arguments.step);

#ifndef  DEBUG

/*************************************************
*                                                *
*        DEFINING THE OSCILLATION ENGINES        *
*                                                *
**************************************************/

  /**********
  *  DANSS  *
  ***********/


glbDefineOscEngine(NUMP, &glf_two_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
			"DANSS_upper", &DANSS_up_s);

glbDefineOscEngine(NUMP, &glf_two_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DANSS_lower", (void *) &DANSS_down_s);
  

  /*************
  *  DAYA BAY  *
  **************/

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH1_AD1", (void *) &DB_EH1_AD1_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH1_AD2", (void *) &DB_EH1_AD2_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
		   	&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH2_AD3", (void *) &DB_EH2_AD3_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH2_AD8", (void *) &DB_EH2_AD8_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD4", (void *) &DB_EH3_AD4_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD5", (void *) &DB_EH3_AD5_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD6", (void *) &DB_EH3_AD6_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DayaBay_EH3_AD7", (void *) &DB_EH3_AD7_s);

  /*********
  *  NEOS  *
  **********/

 glbDefineOscEngine(NUMP, &glf_two_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "NEOS_eng", (void *) &NEOS_s);

  /*****************
  *  DOUBLE CHOOZ  *
  ******************/

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DC_ND", (void *) &DC_ND_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "DC_FD", (void *) &DC_FD_s);

  /**********
  *  BUGEY  *
  ***********/

glbDefineOscEngine(NUMP, &glf_two_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "bugey_near", (void *) &bugey_15_s);

glbDefineOscEngine(NUMP, &glf_two_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "bugey_far", (void *) &bugey_40_s);

  /*********
  *  RENO  *
  **********/

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "RENO_ND", (void *) &RENO_ND_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&glf_get_oscillation_parameters,
			&glf_set_oscillation_parameters,
		   "RENO_FD", (void *) &RENO_FD_s);

#endif

/***********************************************
*                                              *
*        DEFINING CHI-SQUARED FUNCTIONS        *
*                                              *
************************************************/

  int YesNo[10];
  int q;
  for(q=0; q<10; q++){
    YesNo[q]=1;
  }

/* 
  This object tells us which experiments to use or not use; we will flip its
  elements to 0 depending on what has been input at the command line (with a lot of
  if statements...)
*/

  if(arguments.bargf==YES){
    if(arguments.barg[0]==0){
      YesNo[0] = 0;
    }
    if(arguments.barg[1]==0){
      YesNo[1] = 0;
    }
    if(arguments.barg[2]==0){
      YesNo[2] = 0;
    }
    if(arguments.barg[3]==0){
      YesNo[3] = 0;
    }
    if(arguments.barg[4]==0){
      YesNo[4] = 0;
    }
    if(arguments.barg[5]==0){
      YesNo[5] = 0;
    }
    if(arguments.barg[6]==0){
      YesNo[6] = 0;
    }
    if(arguments.barg[7]==0){
      YesNo[7] = 0;
    }
  }

/* A check to make sure that YesNo is being handled properly...
  for(q=0; q<10; q++){
    fprintf(stdout, "%d", YesNo[q]);
  }
  fprintf(stdout, "\n");
*/

  glbDefineChiFunction(&glf_spectra_chi, 4, "spectra-chisq", (void *)YesNo);

/*************************************
*                                    *
*        THE MEAT OF THE CODE        *
*                                    *
**************************************/


  for(i=0;i<glbGetNumOfOscParams();i++) osc[i]=0; 

/* Loading in experiment files... */

  s=glbInitExperiment("glb/spectra.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s2=glbInitExperiment("glb/spectra2.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s3=glbInitExperiment("glb/spectra3.glb",&glb_experiment_list[0],
	&glb_num_of_exps);

  s4=glbInitExperiment("glb/spectra4.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s5=glbInitExperiment("glb/spectra5.glb",&glb_experiment_list[0],
	&glb_num_of_exps);
  s6=glbInitExperiment("glb/spectra6.glb",&glb_experiment_list[0],
	&glb_num_of_exps);

  /* Testing for failure */

  if(s<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s2<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s3<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}

  if(s4<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s5<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}
  if(s6<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",argv[0]);exit(1);}

  /* Lexing the oscillation parameters as given by -p='1,2,3...' */
  if(arguments.params!=NULL){
     char *p=arguments.params;
      for (i=0; i < glbGetNumOfOscParams(); i++){
          if (p != NULL)
            rv=sscanf(p, "%lf", &osc[i]);
          else
            rv=0;
          if (rv < 1){
              fprintf(stderr,"%s: FATAL: Wrong format for oscillation"
                      " parameters\n ",argv[0]);
              exit(1);
            }
          if ((p=strchr(p, ',')) != NULL)
            p++;    /* Go to character after the delimiter */
        }	
    }
  
  /* Initialize parameter vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params starting_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_params minimum = glbAllocParams();

  glb_projection proj=glbAllocProjection();

  for(i=0;i<glbGetNumOfOscParams();i++){glbSetOscParams(true_values,osc[i],i);}
  for(i=0;i<glbGetNumOfOscParams();i++){glbSetOscParams(input_errors,0,i);}
  for(i=0;i<glbGetNumOfOscParams();i++){glbSetProjectionFlag(proj,arguments.marg[i],i);}

  glbSetDensityProjectionFlag(proj,GLB_FIXED,GLB_ALL);
  glbCopyParams(true_values,test_values);  

  glbSetOscParams(test_values, asin(sqrt(0.084))/2.0, GLF_THETA_13);
  glbSetOscParams(test_values, 0.0, GLF_THETA_14);
  glbSetOscParams(test_values, 1.0, GLF_DELTA_M);
  glbSetOscParams(test_values, 2.525e-3, GLF_DELTA_ATM);
  glbSetOscParams(test_values, 0.590, GLF_THETA_12);
  glbSetOscParams(test_values, 7.39e-5, GLF_DELTA_SOLAR);
  
  glbCopyParams(test_values,starting_values);
 
  glbSetDensityParams(input_errors,0,GLB_ALL);
  glbSetInputErrors(input_errors);
  glbSetStartingValues(starting_values);

  /* The simulated data are computed */
  glbSetOscillationParameters(true_values);
  glbSetRates();
  double res;

  double thetheta;
  /*switching off the atmospheric prior */ 
  glbSetDensityParams(input_errors,0,GLB_ALL);
  glbSetInputErrors(input_errors);
  double x,y,res2;

/*************************************************************************

  Initializing energy-scale uncertainty systematics...                 
  Also initializing the minimization over these nuisance parameters... 

*************************************************************************/

  double *scales;

  scales = glbGetSysErrorsListPtr(0, 0, GLB_ON);

  double SetSys[4];
  for (i=0; i<4; i++){
    SetSys[i] = 0.5*scales[i];
  }
  glbSetSysStartingValuesList(0, 0, GLB_ON, &SetSys);

/* 
   We define this object now so that we can read out the ending values of
  the nuisance parameters...
*/

double *EndNuisances;

/*************************************************************************/
  
  if(arguments.margf==YES)
    {
      glbSetProjection(proj);
      res=glbChiNP(true_values,minimum,GLB_ALL);
      fprintf(output,"%f\n",res);
      glbPrintParams(output,minimum);
    }
  else if (arguments.pargf==YES)
    {
      for(i=0;i<glbGetNumOfOscParams();i++){glbSetProjectionFlag(proj,arguments.parg[i],i);}
      glbSetDensityProjectionFlag(proj,GLB_FIXED,GLB_ALL);
      glbSetProjection(proj);
      res=glbChiNP(true_values,minimum,GLB_ALL);
      fprintf(output,"%f\n",res);
      glbPrintParams(output,minimum);
    }

  else if(arguments.chi==YES){
      res=glbChiSys(true_values,GLB_ALL,GLB_ALL);
      fprintf(output,"%f\n",res);
    }

/***************************************************************

  This block spits out the raw GLoBES rates...

***************************************************************/

  else if(arguments.rate==YES){
    double *rate; 
    int rule;
    int bin;

    glbChiSys(test_values,GLB_ALL,GLB_ALL);

    for(i=23;i<29;i++){
      for(rule=0;rule<glbGetNumberOfRules(i);rule++){
        rate= glbGetSignalFitRatePtr(i,rule);

	for(bin=0;bin<glbGetNumberOfBins(i);bin++){
	  fprintf(output,"%f\t",100000.*rate[bin]);
        }
        fprintf(output,"\n");
      }
    }
  }

/***************************************************************

  This block of code executes a scan over theta13 as well as theta14 and Dm241

***************************************************************/

  else if(arguments.doth13==YES){

     resolution=resolution-1;

     double s22t13min, s22t13max, s22t13step, s22t13, th13;
     s22t13min = 0.0;
     s22t13max = 0.2;
     s22t13step = 0.0025;

     for (s22t13 = s22t13min; s22t13<s22t13max+ s22t13step; s22t13 += s22t13step){
      th13=asin(sqrt(s22t13))/2.0;
      glbSetOscParams(test_values, th13, GLF_THETA_13);

      for(x=xrange[0];x<=xrange[1]+(xrange[1]-xrange[0])/resolution/2;
	  x=x+(xrange[1]-xrange[0])/resolution)
	for(y=yrange[0];y<=yrange[1]+(yrange[1]-yrange[0])/resolution/2;
	    y=y+(yrange[1]-yrange[0])/resolution){

	    /* Set vector of test values */
	    if(x<0) thetheta=asin(sqrt(pow(10,x)))/2.0;
	    else thetheta=M_PI/4.0;
	    
#ifdef DEBUG
	    glbSetOscParams(test_values,thetheta,0);
	    glbSetOscParams(test_values,pow(10,y),1);
#else 
	    glbSetOscParams(test_values,thetheta,GLF_THETA_14);
	    glbSetOscParams(test_values,pow(10,y),GLF_DELTA_M);
#endif

	    /* Compute Chi^2 for all loaded experiments and all rules */
	    res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
	    fprintf(output,"%f %f %f %f\n", s22t13,x,y,res);
	  }
     }
    }

/***************************************************************

  This block of code executes a scan over theta13 and Dm_{31}^2

***************************************************************/

  else if(arguments.calibrate == YES){
    double s22t13min, s22t13max, s22t13step, s22t13, th13;
    s22t13min = 0.05;
    s22t13max = 0.2;
    s22t13step = 0.001;

    double Dm31min, Dm31max, Dm31step, Dm31;
    Dm31min = 1.0;
    Dm31max = 4.0;
    Dm31step = 0.01;

    for (s22t13 = s22t13min; s22t13 <= s22t13max; s22t13 += s22t13step){
      th13=asin(sqrt(s22t13))/2.0;
      glbSetOscParams(test_values, th13, GLF_THETA_13);
	
      for (Dm31 = Dm31min; Dm31 <= Dm31max; Dm31 += Dm31step){
        glbSetOscParams(test_values, Dm31*(1.0e-3), GLF_DELTA_ATM);
	    
        /* Compute Chi^2 for all loaded experiments and all rules */
        res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        fprintf(output,"%f %f %f\n", s22t13, Dm31, res);
        }
      }
    }

/***************************************************************

  This block of code executes a scan over theta14 and Dm241;
  it is the default option for this program

***************************************************************/

  else{
     resolution=resolution-1;
      for(x=xrange[0];x<=xrange[1]+(xrange[1]-xrange[0])/resolution/2;
	  x=x+(xrange[1]-xrange[0])/resolution)
	for(y=yrange[0];y<=yrange[1]+(yrange[1]-yrange[0])/resolution/2;
	    y=y+(yrange[1]-yrange[0])/resolution){

	    /* Set vector of test values */
	    if(x<0) thetheta=asin(sqrt(pow(10,x)))/2.0;
	    else thetheta=M_PI/4.0;
	    
#ifdef DEBUG
	    glbSetOscParams(test_values,thetheta,0);
	    glbSetOscParams(test_values,pow(10,y),1);
#else 

	    glbSetOscParams(test_values,thetheta,GLF_THETA_14);
	    glbSetOscParams(test_values,pow(10,y),GLF_DELTA_M);
#endif

	    /* Compute Chi^2 for all loaded experiments and all rules */
	    res=glbChiSys(test_values,GLB_ALL,GLB_ALL);

	    /* Get nuisance parameters; write... */
            EndNuisances = glbGetSysStartingValuesListPtr(0, 0, GLB_ON);
	    fprintf(output,"%f %f %f\n",
               x,y,res);

	  }
    }

  /* Destroy parameter vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(starting_values);
  glbFreeParams(input_errors); 
  glbFreeParams(minimum); 
  fclose(output);
  exit(0);
}
