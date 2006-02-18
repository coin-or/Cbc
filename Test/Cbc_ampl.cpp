/****************************************************************
Copyright (C) 1997-2000 Lucent Technologies
Modifications for Coin -  Copyright (C) 2006, International Business Machines Corporation and others.
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/
#ifdef CBC_AMPL
#include "getstub.h"
#include "Cbc_ampl.h"
#include "unistd.h"
#include <string>
/* so decodePhrase and clpCheck can access */
static ampl_info * saveInfo=NULL;
// Set to 1 if algorithm found
static char algFound[20]="";
static char*
checkPhrase(Option_Info *oi, keyword *kw, char *v)
{
  if (strlen(v))
    printf("string %s\n",v);
  // Say algorithm found
  strcpy(algFound,kw->desc);;
  return v;
}
static char*
checkPhrase2(Option_Info *oi, keyword *kw, char *v)
{
  if (strlen(v))
    printf("string %s\n",v);
  // put out keyword
  saveInfo->arguments=(char **) realloc(saveInfo->arguments,(saveInfo->numberArguments+1)*sizeof(char *));
  saveInfo->arguments[saveInfo->numberArguments++]=strdup(kw->desc);
  return v;
}
static fint
decodePhrase(char * phrase,ftnlen length)
{
  char * blank = strchr(phrase,' ');
  if (blank) {
    /* split arguments */
    *blank='\0';
    saveInfo->arguments=(char **) realloc(saveInfo->arguments,(saveInfo->numberArguments+2)*sizeof(char *));
    saveInfo->arguments[saveInfo->numberArguments++]=strdup(phrase);
    *blank=' ';
    phrase=blank+1; /* move on */
    saveInfo->arguments[saveInfo->numberArguments++]=strdup(phrase);
  } else {
    saveInfo->arguments=(char **) realloc(saveInfo->arguments,(saveInfo->numberArguments+1)*sizeof(char *));
    saveInfo->arguments[saveInfo->numberArguments++]=strdup(phrase);
  }
  return 0;
}
static char xxxxxx[20];
#define VP (char*)
 static keyword keywds[] = { /* must be sorted */
	{ "barrier",	checkPhrase,		(char *) xxxxxx ,"-barrier" },
	{ "dual",	checkPhrase,		(char *) xxxxxx , "-dualsimplex"},
	{ "help",	checkPhrase2,		(char *) xxxxxx , "-?"},
	{ "initial",	checkPhrase,		(char *) xxxxxx , "-initialsolve"},
	{ "max",	checkPhrase2,		(char *) xxxxxx , "-maximize"},
	{ "maximize",	checkPhrase2,		(char *) xxxxxx , "-maximize"},
	{ "primal",	checkPhrase,		(char *) xxxxxx , "-primalsimplex"},
	{ "quit",	checkPhrase,		(char *) xxxxxx , "-quit"},
	{ "wantsol",	WS_val,		NULL, "write .sol file (without -AMPL)" }
	};
static Option_Info Oinfo = {"cbc", "Cbc 1.01", "cbc_options", keywds, nkeywds, 0, "",
				0,decodePhrase,0,0,0, 20060130 };
// strdup used to avoid g++ compiler warning
 static SufDecl
suftab[] = {
#if 0
	{ "current", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "current", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "direction", 0, ASL_Sufkind_var },
	{ "down", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "down", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "priority", 0, ASL_Sufkind_var },
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sos", 0, ASL_Sufkind_var },
	{ "sos", 0, ASL_Sufkind_con },
	{ "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
#endif
	{ strdup("sstatus"), 0, ASL_Sufkind_var, 0 },
	{ strdup("sstatus"), 0, ASL_Sufkind_con, 0 }
#if 0
	{ "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_outonly},
	{ "up", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "up", 0, ASL_Sufkind_var | ASL_Sufkind_outonly }
#endif
	};
#include "float.h"
#include "limits.h"
static ASL *asl=NULL;
static FILE *nl=NULL;

static void
stat_map(int *stat, int n, int *map, int mx, const char *what)
{
  int bad, i, i1=0, j, j1=0;
  static char badfmt[] = "Coin driver: %s[%d] = %d\n";
  
  for(i = bad = 0; i < n; i++) {
    if ((j = stat[i]) >= 0 && j <= mx)
      stat[i] = map[j];
    else {
      stat[i] = 0;
      i1 = i;
      j1 = j;
      if (!bad++)
        fprintf(Stderr, badfmt, what, i, j);
    }
  }
  if (bad > 1) {
    if (bad == 2)
      fprintf(Stderr, badfmt, what, i1, j1);
    else
      fprintf(Stderr,
              "Coin driver: %d messages about bad %s values suppressed.\n",
              bad-1, what);
  }
}
int
readAmpl(ampl_info * info, int argc, char **argv)
{
  char *stub;
  ograd *og;
  int i;
  SufDesc *csd;
  SufDesc *rsd;
  /*bool *basis, *lower;*/
  /*double *LU, *c, lb, objadj, *rshift, *shift, t, ub, *x, *x0, *x1;*/
  char * environment = getenv("cbc_options");
  char tempBuffer[20];
  double * obj;
  double * columnLower;
  double * columnUpper;
  double * rowLower;
  double * rowUpper;
  char ** saveArgv=argv;
  int saveArgc = argc;
  memset(info,0,sizeof(ampl_info));
  /* save so can be accessed by decodePhrase */
  saveInfo = info;
  info->numberArguments=0;
  info->arguments=(char **) malloc(2*sizeof(char *));
  info->arguments[info->numberArguments++]=strdup("ampl");
  info->arguments[info->numberArguments++]=strdup("cbc");
  asl = ASL_alloc(ASL_read_f);
  stub = getstub(&argv, &Oinfo);
  if (!stub)
    usage_ASL(&Oinfo, 1);
  nl = jac0dim(stub, 0);
  /*void * specialOrderedInfo = sos_add(nl,0);*/
  suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
  
  /* set A_vals to get the constraints column-wise (malloc so can be freed) */
  A_vals = (double *) malloc(nzc*sizeof(double));
  if (!A_vals) {
    printf("no memory\n");
    return 1;
  }
  /* say we want primal solution */
  want_xpi0=1;
  /* for basis info */
  info->columnStatus = (int *) malloc(n_var*sizeof(int));
  info->rowStatus = (int *) malloc(n_con*sizeof(int));
  csd = suf_iput("sstatus", ASL_Sufkind_var, info->columnStatus);
  rsd = suf_iput("sstatus", ASL_Sufkind_con, info->rowStatus);
  /* read linear model*/
  f_read(nl,0);

  /*sos_finish(&specialOrderedInfo, 0, &j, 0, 0, 0, 0, 0);*/
  Oinfo.uinfo = tempBuffer;
  if (getopts(argv, &Oinfo))
    return 1;
  /* objective*/
  obj = (double *) malloc(n_var*sizeof(double));
  for (i=0;i<n_var;i++)
    obj[i]=0.0;;
  if (n_obj) {
    for (og = Ograd[0];og;og = og->next)
      obj[og->varno] = og->coef;
  }
  if (objtype[0])
    info->direction=-1.0;
  else
    info->direction=1.0;
  info->offset=objconst(0);
  /* Column bounds*/
  columnLower = (double *) malloc(n_var*sizeof(double));
  columnUpper = (double *) malloc(n_var*sizeof(double));
#define COIN_DBL_MAX DBL_MAX
  for (i=0;i<n_var;i++) {
    columnLower[i]=LUv[2*i];
    if (columnLower[i]<= negInfinity)
      columnLower[i]=-COIN_DBL_MAX;
    columnUpper[i]=LUv[2*i+1];
    if (columnUpper[i]>= Infinity)
      columnUpper[i]=COIN_DBL_MAX;
  }
  /* Row bounds*/
  rowLower = (double *) malloc(n_con*sizeof(double));
  rowUpper = (double *) malloc(n_con*sizeof(double));
  for (i=0;i<n_con;i++) {
    rowLower[i]=LUrhs[2*i];
    if (rowLower[i]<= negInfinity)
      rowLower[i]=-COIN_DBL_MAX;
    rowUpper[i]=LUrhs[2*i+1];
    if (rowUpper[i]>= Infinity)
      rowUpper[i]=COIN_DBL_MAX;
  }
  info->numberRows=n_con;
  info->numberColumns=n_var;
  info->numberElements=nzc;;
  info->numberBinary=nbv;
  info->numberIntegers=niv;
  info->objective=obj;
  info->rowLower=rowLower;
  info->rowUpper=rowUpper;
  info->columnLower=columnLower;
  info->columnUpper=columnUpper;
  info->starts=A_colstarts;
  /*A_colstarts=NULL;*/
  info->rows=A_rownos;
  /*A_rownos=NULL;*/
  info->elements=A_vals;
  /*A_vals=NULL;*/
  info->primalSolution=NULL;
  /* put in primalSolution if exists */
  if (X0) {
    info->primalSolution=(double *) malloc(n_var*sizeof(double));
    memcpy(info->primalSolution,X0,n_var*sizeof(double));
  }
  info->dualSolution=NULL;
  if ((!(niv+nbv)&&(csd->kind & ASL_Sufkind_input))
      ||(rsd->kind & ASL_Sufkind_input)) {
    /* convert status - need info on map */
    static int map[] = {1, 3, 1, 1, 2, 1, 1};
    stat_map(info->columnStatus, n_var, map, 6, "incoming columnStatus");
    stat_map(info->rowStatus, n_con, map, 6, "incoming rowStatus");
  } else {
    /* all slack basis */
    free(info->rowStatus);
    info->rowStatus=NULL;
    free(info->columnStatus);
    info->columnStatus=NULL;
  }
  /* add -solve - unless something there already
   - also check for sleep=yes */
  {
    int found=0;
    int foundLog=0;
    int foundSleep=0;
    const char * something[]={"solve","branch","duals","primals"};
    for (i=0;i<info->numberArguments;i++) {
      unsigned int j;
      const char * argument = info->arguments[i];
      for (j=0;j<sizeof(something)/sizeof(char *);j++) {
        const char * check = something[j];
        if (!strncmp(argument,check,sizeof(check))) {
          found=(int)(j+1);
        } else if (!strncmp(argument,"log",3)) {
          foundLog=1;
        } else if (!strncmp(argument,"sleep",5)) {
          foundSleep=1;
        }
      }
    }
    if (foundLog) {
      /* print options etc */
      for (i=0;i<saveArgc;i++)
        printf("%s ",saveArgv[i]);
      printf("\n");
      if (environment)
        printf("env %s\n",environment);
      /*printf("%d rows %d columns %d elements\n",n_con,n_var,nzc);*/
    }
    if (!found) {
      if (!strlen(algFound)) {
        info->arguments=(char **) realloc(info->arguments,(info->numberArguments+1)*sizeof(char *));
        info->arguments[info->numberArguments++]=strdup("-solve");
      } else {
        // use algorithm from keyword
        info->arguments=(char **) realloc(info->arguments,(info->numberArguments+1)*sizeof(char *));
        info->arguments[info->numberArguments++]=strdup(algFound);
      }
    }
    if (foundSleep) {
      /* let user copy .nl file */
      fprintf(stderr,"You can copy .nl file %s for debug purposes or attach debugger\n",saveArgv[1]);
      fprintf(stderr,"Type q to quit, anything else to continue\n");
      char getChar = getc(stdin);
      if (getChar=='q'||getChar=='Q')
        exit(1);
    }
  }
  /* add -quit */
  info->arguments=(char **) realloc(info->arguments,(info->numberArguments+1)*sizeof(char *));
  info->arguments[info->numberArguments++]=strdup("-quit");
  return 0;
}
void freeArrays1(ampl_info * info)
{
  free(info->objective);
  info->objective=NULL;
  free(info->rowLower);
  info->rowLower=NULL;
  free(info->rowUpper);
  info->rowUpper=NULL;
  free(info->columnLower);
  info->columnLower=NULL;
  free(info->columnUpper);
  info->columnUpper=NULL;
  /* this one not freed by ASL_free */
  free(info->elements);
  info->elements=NULL;
  free(info->primalSolution);
  info->primalSolution=NULL;
  free(info->dualSolution);
  info->dualSolution=NULL;
  /*free(info->rowStatus);
  info->rowStatus=NULL;
  free(info->columnStatus);
  info->columnStatus=NULL;*/
}
void freeArrays2(ampl_info * info)
{
  free(info->primalSolution);
  info->primalSolution=NULL;
  free(info->dualSolution);
  info->dualSolution=NULL;
  free(info->rowStatus);
  info->rowStatus=NULL;
  free(info->columnStatus);
  info->columnStatus=NULL;
  ASL_free(&asl);
}
void freeArgs(ampl_info * info)
{
  int i;
  for ( i=0; i<info->numberArguments;i++)
    free(info->arguments[i]);
  free(info->arguments);
}
int ampl_obj_prec()
{
  return obj_prec();
}
void writeAmpl(ampl_info * info)
{
  char buf[1000];
  typedef struct { const char *msg; int code; int wantObj; } Sol_info;
  static Sol_info solinfo[] = {
    { "optimal solution",			000, 1 },
    { "infeasible",     			200, 1 },
    { "unbounded",	        		300, 0 },
    { "iteration limit etc",			400, 1 },
    { "solution limit",				401, 1 },
    { "ran out of space",			500, 0 },
    { "status unknown",				501, 1 },
    { "bug!",					502, 0 },
    { "best MIP solution so far restored",	101, 1 },
    { "failed to restore best MIP solution",	503, 1 },
    { "optimal (?) solution",			100, 1 }
  };
  /* convert status - need info on map */
  static int map[] = {0, 3, 4, 1};
  sprintf(buf,"%s %s",Oinfo.bsname,info->buffer);
  solve_result_num = solinfo[info->problemStatus].code;
  if (info->columnStatus) {
    stat_map(info->columnStatus, n_var, map, 4, "outgoing columnStatus");
    stat_map(info->rowStatus, n_con, map, 4, "outgoing rowStatus");
    suf_iput("sstatus", ASL_Sufkind_var, info->columnStatus);
    suf_iput("sstatus", ASL_Sufkind_con, info->rowStatus);
  }
  write_sol(buf,info->primalSolution,info->dualSolution,&Oinfo);
}

#endif
