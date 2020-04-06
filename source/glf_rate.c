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

/* (C) 2005, 2007, 2019 Patrick Huber, J. M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <argp.h>

#include <globes/globes.h>   /* GLoBES library */

#include "glf_types.h"
#include "glf_rate_chi.h"
#include "glf_probability.h"
#include "glf_precomputed_probabilities.h"

#define VERSION "1.0"
#define MAXNUM 12 /* largest possible number of oscillation parameters */
#define GT GITVERSION

double theta_14;
double delta_m; /* Delta m_{41}^2 */
double delta_atm; /* Delta m_{31}^2 */
double theta_13;
double theta_12;
double delta_solar;

/* the infamous exponent to E_res^N */
int N;

double systematic[10];

const char *argp_program_version =
"GLoBESfit_rate "VERSION"\n(C) 2007, 2019, 2020 J. M. Berryman, P. Huber \n"
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
  {"Project",'P',"DIRECTION",0,"Project onto DIRECTION",0}, 
  {"chi",'c',0,0,"Compute chi^2 for parameters given by -p",0},
  {"Initial-step",'I',"STEPSIZE",0,"Initial step size for minimizer",0},
  {"Systematics",'S',0,0,"Systematic Uncertainties on Reactor Fluxes (Default is OFF)",0},
  {"block",'b',"BLOCK", 0, "Experiments that are **NOT** \n to be included in present analysis; \n see documentation for definitions.",0},
  {"Unfix", 'u',0,0, "Toggle on unfixing the U235, Pu239 fluxes (Default is OFF)" ,0},
  {"SM", 'M',0,0, "Toggle whether to use summation method flux (Default is OFF)" ,0},
  {"HKSS", 'H',0,0, "Toggle whether to use the HKSS \n correction to the HM fluxes (Default is OFF)" ,0},
  {"Write Nuisance", 'w',0,0, "Toggle whether to write nuisance parameters (Default is OFF)" ,0},
  { 0 } 
};

/* Used by `main' to communicate with `parse_opt'. */
struct arguments
{
/*  int exp; */ /* number of experiments */
  char* args[32];                /* up to 32 experiments */
  int resolution;
  int marg[MAXNUM]; 
  int parg[MAXNUM];
  int barg[16]; /* Array for blocks of experiments being disregarded */
  double step;
  int margf,pargf,Systematics,bargf, Unfix, SumMethod, HKSS, write;
  char *output_file;
  char *params,*xrange,*yrange;
  int chi;
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

    case 'S':
      arguments->Systematics = YES;
      break; 

    case 'c':
      arguments->chi = YES;
      break;

    case 'w':
      arguments->write = YES;
      break;

    case 'u':
      arguments->Unfix = YES;
      break; 

    case 'M':
      arguments->SumMethod = YES;
      break; 

    case 'H':
      arguments->HKSS = YES;
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
int main(int argc, char *argv[])
{

  double progress;
  int count=0;
  int s,s2,s3,i,rv;
  FILE *output=NULL;
  int resolution;

  int NUMP=6;

  double *osc=NULL;
  double xrange[]={-4,0};
  double yrange[]={-2,2};
  struct arguments arguments;

  arguments.args[0]="-";
  arguments.resolution=51;
  arguments.output_file=NULL;
  arguments.params=NULL;
  arguments.xrange=NULL;
  arguments.yrange=NULL;

  for(i=0;i<MAXNUM;i++)  arguments.marg[i]=GLB_FIXED;
  for(i=0;i<3;i++)  arguments.parg[i]=GLB_FIXED;
  for(i=3;i<MAXNUM;i++)  arguments.parg[i]=GLB_FREE;
  for(i=0;i<16;i++)  arguments.barg[i]=1;

  arguments.margf=NO;
  arguments.pargf=NO;
  arguments.bargf=NO;

  arguments.chi=NO;
  arguments.step=0.05;

  arguments.Systematics=NO;
  arguments.Unfix=NO;
  arguments.SumMethod=NO;
  arguments.HKSS=NO;
  arguments.write=NO;

  /* parsing the command line */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);  
  
  N=1;

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

/* Oscillation engine for Bugey/Rovno91...remade for "rate" calculation */
	    
double bugey_wide[1] = {0.0015};
		    
glbDefineOscEngine(NUMP, &combo_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "bugey_rate"
		    ,(void *) &bugey_wide);

/* Oscillation engine for Gosgen...NEED ACTUAL GEOMETRY*/
	    
double gosgen_wide[1] = {0.001};
		    
glbDefineOscEngine(NUMP, &combo_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "gosgen_rate"
		    ,(void *) & gosgen_wide);

/* Oscillation engine for Krasnoyarsk...NEED ACTUAL GEOMETRY*/

double krasnoyarsk_wide[1] = {0.001};
		    
glbDefineOscEngine(NUMP, &combo_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "krasnoyarsk_rate"
		    ,(void *) &krasnoyarsk_wide);

/* Oscillation engine for Savannah River...NEED ACTUAL GEOMETRY*/

double SRP_wide[1] = {0.001};
		    
glbDefineOscEngine(NUMP, &combo_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "SRP_rate"
		    ,(void *) &SRP_wide);

/* Oscillation engine for Rovno88...NEED ACTUAL GEOMETRY*/

double Rovno88_wide[1] = {0.0005};
		    
glbDefineOscEngine(NUMP, &combo_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "rovno88_rate"
		    ,(void *) &Rovno88_wide);

/* Oscillation engine for Nucifer...NEED ACTUAL GEOMETRY*/

double Nucifer_wide[1] = {0.0015};
		    
glbDefineOscEngine(NUMP, &combo_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "nucifer_rate"
		    ,(void *) &Nucifer_wide);

/* Standard oscillation engine */

/*
  A fictitious width, to regulate high-mass oscillations at medium- and
  long-baseline reactor experiments...
*/

double FakeWide[1] = {0.002};
		    
glbDefineOscEngine(NUMP, &standard_probability_matrix,
		    &combo_get_oscillation_parameters,
		    &combo_set_oscillation_parameters,
		    "standard"
		    ,(void *) &FakeWide);

  /*************
  *  DAYA BAY  *
  **************/
 
glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&combo_get_oscillation_parameters,
			&combo_set_oscillation_parameters,
		   "DayaBay_EH1_AD1", (void *) &DB_EH1_AD1_s);
 

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&combo_get_oscillation_parameters,
			&combo_set_oscillation_parameters,
		   "DayaBay_EH1_AD2", (void *) &DB_EH1_AD2_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&combo_get_oscillation_parameters,
			&combo_set_oscillation_parameters,
		   "DayaBay_EH2_AD3", (void *) &DB_EH2_AD3_s);

glbDefineOscEngine(NUMP, &glf_four_state_probability_matrix,
			&combo_get_oscillation_parameters,
			&combo_set_oscillation_parameters,
		   "DayaBay_EH2_AD8", (void *) &DB_EH2_AD8_s);


#endif

/***********************************************
*                                              *
*        DEFINING CHI-SQUARED FUNCTIONS        *
*                                              *
************************************************/

  int YesNo[26];
  int q;
  for(q=0; q<26; q++){
    YesNo[q]=1;
  }

/* 
  This object tells us which experiments to use or not use; we will flip its
  elements to 0 depending on what has been input at the command line (with a lot of
  if statements...)

  There's got to be a smarter way to do this...I'll think on it...
*/

  if(arguments.bargf==YES){
    if(arguments.barg[0]==0){
      YesNo[0] = 0;
      YesNo[1] = 0;
    }
    if(arguments.barg[1]==0){
      YesNo[2] = 0;
      YesNo[3] = 0;
      YesNo[4] = 0;
    }
    if(arguments.barg[2]==0){
      YesNo[5] = 0;
      YesNo[6] = 0;
      YesNo[7] = 0;
      YesNo[8] = 0;
    }
    if(arguments.barg[3]==0){
      YesNo[9] = 0;
      YesNo[10] = 0;
    }
    if(arguments.barg[4]==0){
      YesNo[11] = 0;
    }
    if(arguments.barg[5]==0){
      YesNo[12] = 0;
    }
    if(arguments.barg[6]==0){
      YesNo[13] = 0;
    }
    if(arguments.barg[7]==0){
      YesNo[14] = 0;
    }
    if(arguments.barg[8]==0){
      YesNo[15] = 0;
      YesNo[16] = 0;
      YesNo[17] = 0;
      YesNo[18] = 0;
      YesNo[19] = 0;
    }
    if(arguments.barg[9]==0){
      YesNo[20] = 0;
    }
    if(arguments.barg[10]==0){
      YesNo[21] = 0;
    }
    if(arguments.barg[11]==0){
      YesNo[22] = 0;
    }
    if(arguments.barg[12]==0){
      YesNo[23] = 0;
    }
    if(arguments.barg[13]==0){
      YesNo[24] = 0;
    }
    if(arguments.barg[14]==0){
      YesNo[25] = 0;
    }
  }

/* A check to make sure that YesNo is being handled properly...
  for(q=0; q<26; q++){
    fprintf(stdout, "%d", YesNo[q]);
  }
  fprintf(stdout, "\n");
*/

/*
  This array contains the U235 & Pu239 rescaling factors, in that order, followed by
  the elements of the array "YesNo"
*/

  double Ratios[28];
  Ratios[0] = 1.0;
  Ratios[1] = 1.0;
  for(q=2; q<28; q++){
    Ratios[q]=1.0*YesNo[q-2];
  }

/*
  for(q=0; q<28; q++){
    fprintf(stdout, "%f", Ratios[q]);
  }
  fprintf(stdout, "\n");
*/

  glbDefineChiFunction(&combo_rate_chi, 4, "rate-combo", (void *)YesNo);
  glbDefineChiFunction(&combo_rate_chi_unfix, 2, "rate-combo-unfix", (void *)Ratios);
  glbDefineChiFunction(&combo_rate_chi_nosys, 0, "rate-combo-no-sys", (void *)YesNo);

  glbDefineChiFunction(&combo_rate_chi_SM, 4, "rate-combo-SM", (void *)YesNo);
  glbDefineChiFunction(&combo_rate_chi_HKSS, 4, "rate-combo-HKSS", (void *)YesNo);

/*************************************
*                                    *
*        THE MEAT OF THE CODE        *
*                                    *
**************************************/

 /* memory for osc */
  osc = (double *) glb_malloc(sizeof(*osc) * glbGetNumOfOscParams());
  for(i=0;i<glbGetNumOfOscParams();i++) osc[i]=0; /* FIXME glb_calloc
						     is missing */

/*
  IMPORTANT REMARK: There are five distinct analyses that can be performed using
  this code, as things currently stand, according to the presence of certain flags.
  The flags, and what they do are as follows:

  -H: Uses the HKSS corrections to the HM fluxes, taken from 1908.08302

  -M: Uses fluxes from summation method calculation of 1904.09358 to scan over
	sterile neutrino parameter space; no systematics are included.

  -U: Unfix the overall HM fluxes for U235 and Pu239 and scan over these

  -S: Scan over sterile neutrino parameter space and *include* systematics on the
	Huber-Mueller fluxes

  <none of these>: Scan over sterile neutrino parameter space and *ignore* 
	systematics on the Huber-Mueller fluxes

  These are listed in order of precedence, i.e., -M and -S causes the code to run
  in "-M" mode
*/

  if (arguments.HKSS==YES){
  /* Sterile oscillations using HKSS FLUXES */
    s=glbInitExperiment("glb/rate_combo_HKSS.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s2=glbInitExperiment("glb/rate_combo_HKSS2.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s3=glbInitExperiment("glb/rate_combo_HKSS3.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
  }
  else if (arguments.SumMethod==YES){
  /* Sterile oscillations using SUMMATION METHOD FLUXES */
    s=glbInitExperiment("glb/rate_combo_SM.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s2=glbInitExperiment("glb/rate_combo_SM2.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s3=glbInitExperiment("glb/rate_combo_SM3.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
  }
  else if (arguments.Unfix==YES){
  /* Unfix U235 and Pu239 fluxes; NO OSCILLATIONS HERE */
    s=glbInitExperiment("glb/rate_combo_unfix.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s2=glbInitExperiment("glb/rate_combo2.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s3=glbInitExperiment("glb/rate_combo3.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
  }
  else if (arguments.Systematics==YES){
  /* Sterile oscillations WITH flux prediction systematics */
    s=glbInitExperiment("glb/rate_combo.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s2=glbInitExperiment("glb/rate_combo2.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s3=glbInitExperiment("glb/rate_combo3.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
  }
  else{
  /* Sterile oscillations WITHOUT flux prediction systematics */
    s=glbInitExperiment("glb/rate_combo_no_sys.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s2=glbInitExperiment("glb/rate_combo2.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
    s3=glbInitExperiment("glb/rate_combo3.glb",&glb_experiment_list[0],
		&glb_num_of_exps);
  }

  /* Testing for failure */
  if(s<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error (1)\n",argv[0]);exit(1);}
  if(s2<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error (2)\n",argv[0]);exit(1);}
  if(s3<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error (3)\n",argv[0]);exit(1);}

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

  glbSetOscParams(true_values, asin(sqrt(0.084))/2.0, MY_THETA_13);
  glbSetOscParams(true_values, 0.0, MY_THETA_14);
  glbSetOscParams(true_values, 1.0, MY_DELTA_M);
  glbSetOscParams(true_values, 2.525e-3, MY_DELTA_ATM);
  glbSetOscParams(true_values, 0.590, MY_THETA_12);
  glbSetOscParams(true_values, 7.39e-5, MY_DELTA_SOLAR);

  glbSetDensityProjectionFlag(proj,GLB_FIXED,GLB_ALL);
  glbCopyParams(true_values,test_values);  

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
  else if (arguments.Unfix==YES){
    for(x=0.6; x<=1.6; x=x+0.01){
      for(y=0.6; y<=1.4; y=y+0.01){
	glbSetOscParams(test_values,0.0,MY_THETA_14);
	glbSetOscParams(test_values,0.0,MY_DELTA_M);

        Ratios[0]=x;
        Ratios[1]=y;

        res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        fprintf(output,"%f %f %f\n",x,y,res);
      }
    }
  }

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

	    glbSetOscParams(test_values,thetheta,MY_THETA_14);
	    glbSetOscParams(test_values,pow(10,y),MY_DELTA_M);
#endif
	    count++;
	   
	    progress=count/pow(resolution+1,2)*100;
	    /* Compute Chi^2 for all loaded experiments and all rules */
	    res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
            if (arguments.write == YES){ fprintf(output,"%f %f %f %f %f %f %f \n",
			x,y,res, systematic[0], systematic[1], 
			systematic[2], systematic[3]);}
            else{fprintf(output,"%f %f %f\n",x,y,res); }
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
