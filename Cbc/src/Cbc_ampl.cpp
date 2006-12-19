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

#include "CbcConfig.h"
#ifdef HAVE_UNISTD_H
# include "unistd.h"
#endif
#include "CoinUtilsConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinModel.hpp"
#include "CoinSort.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinMpsIO.hpp"
#include "CoinFloatEqual.hpp"

extern "C" {
# include "getstub.h"
}

#include "Cbc_ampl.h"
#include <string>
#include <cassert>
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
    if (strlen(phrase))
      saveInfo->arguments[saveInfo->numberArguments++]=strdup(phrase);
  } else if (strlen(phrase)) {
    saveInfo->arguments=(char **) realloc(saveInfo->arguments,(saveInfo->numberArguments+1)*sizeof(char *));
    saveInfo->arguments[saveInfo->numberArguments++]=strdup(phrase);
  }
  return 0;
}
static void
sos_kludge(int nsos, int *sosbeg, double *sosref,int * sosind)
{
  // Adjust sosref if necessary to make monotonic increasing
  int i, j, k;
  // first sort
  for (i=0;i<nsos;i++) {
    k = sosbeg[i];
    int end=sosbeg[i+1];
    CoinSort_2(sosref+k,sosref+end,sosind+k);
  }
  double t, t1;
  for(i = j = 0; i++ < nsos; ) {
    k = sosbeg[i];
    t = sosref[j];
    while(++j < k) {
      t1 = sosref[j];
      t += 1e-10;
      if (t1 <= t)
        sosref[j] = t1 = t + 1e-10;
      t = t1;
    }
  }
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
	{ "quit",	checkPhrase2,		(char *) xxxxxx , "-quit"},
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
#endif
	{ "direction", 0, ASL_Sufkind_var },
	{ "downPseudocost", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "priority", 0, ASL_Sufkind_var },
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sos", 0, ASL_Sufkind_var },
	{ "sos", 0, ASL_Sufkind_con },
	{ "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ strdup("sstatus"), 0, ASL_Sufkind_var, 0 },
	{ strdup("sstatus"), 0, ASL_Sufkind_con, 0 },
	{ "upPseudocost", 0, ASL_Sufkind_var | ASL_Sufkind_real }
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
mip_stuff(void)
{
  int i;
  double *pseudoUp, *pseudoDown;
  int *priority, *direction;
  SufDesc *dpup, *dpdown, *dpri, *ddir;
  
  ddir = suf_get("direction", ASL_Sufkind_var);
  direction = ddir->u.i;
  dpri = suf_get("priority", ASL_Sufkind_var);
  priority = dpri->u.i;
  dpdown = suf_get("downPseudocost", ASL_Sufkind_var);
  pseudoDown = dpdown->u.r;
  dpup = suf_get("upPseudocost", ASL_Sufkind_var);
  pseudoUp = dpup->u.r;
  assert(saveInfo);
  int numberColumns = saveInfo->numberColumns;
  if (direction) {
    int baddir=0;
    saveInfo->branchDirection = (int *) malloc(numberColumns*sizeof(int));
    for (i=0;i<numberColumns;i++) {
      int value = direction[i];
      if (value<-1||value>1) {
        baddir++;
        value=0;
      }
      saveInfo->branchDirection[i]=value;
    }
    if (baddir)
      fprintf(Stderr,
              "Treating %d .direction values outside [-1, 1] as 0.\n",
              baddir);
  }
  if (priority) {
    int badpri=0;
    saveInfo->priorities= (int *) malloc(numberColumns*sizeof(int));
    for (i=0;i<numberColumns;i++) {
      int value = priority[i];
      if (value<0) {
        badpri++;
        value=0;
      }
      saveInfo->priorities[i]=value;
    }
    if (badpri)
      fprintf(Stderr,
              "Treating %d negative .priority values as 0\n",
              badpri);
  }
  if (pseudoDown||pseudoUp) {
    int badpseudo=0;
    if (!pseudoDown||!pseudoUp)
      fprintf(Stderr,
              "Only one set of pseudocosts - assumed same\n");
    saveInfo->pseudoDown= (double *) malloc(numberColumns*sizeof(double));
    saveInfo->pseudoUp = (double *) malloc(numberColumns*sizeof(double));
    for (i=0;i<numberColumns;i++) {
      double valueD=0.0, valueU=0.0;
      if (pseudoDown) {
        valueD = pseudoDown[i];
        if (valueD<0) {
          badpseudo++;
          valueD=0.0;
        }
      }
      if (pseudoUp) {
        valueU = pseudoUp[i];
        if (valueU<0) {
          badpseudo++;
          valueU=0.0;
        }
      }
      if (!valueD)
        valueD=valueU;
      if (!valueU)
        valueU=valueD;
      saveInfo->pseudoDown[i]=valueD;
      saveInfo->pseudoUp[i]=valueU;
    }
    if (badpseudo)
      fprintf(Stderr,
              "Treating %d negative pseudoCosts as 0.0\n",badpseudo);
  }
}
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
readAmpl(ampl_info * info, int argc, char **argv, void ** coinModel)
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
  char fileName[1000];
  if (argc>1)
    strcpy(fileName,argv[1]);
  else
    fileName[0]='\0';
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
  // testosi parameter - if >= 10 then go in through coinModel
  int testOsi=-1;
  for (i=0;i<saveArgc;i++) {
    if (!strncmp(saveArgv[i],"testosi",7)) {
      testOsi = atoi(saveArgv[i+1]);
      break;
    }
  }
  if (!(nlvc+nlvo)||testOsi>=10) {
    /* read linear model*/
    f_read(nl,0);
    // see if any sos
    if (true) {
      char *sostype;
      int nsosnz, *sosbeg, *sosind, * sospri;
      double *sosref;
      int nsos;
      int i = ASL_suf_sos_explict_free;
      int copri[2], **p_sospri;
      copri[0] = 0;
      copri[1] = 0;
      p_sospri = &sospri;
      nsos = suf_sos(i, &nsosnz, &sostype, p_sospri, copri,
		     &sosbeg, &sosind, &sosref);
      if (nsos) {
	info->numberSos=nsos;
	info->sosType = (char *) malloc(nsos);
	info->sosPriority = (int *) malloc(nsos*sizeof(int));
	info->sosStart = (int *) malloc((nsos+1)*sizeof(int));
	info->sosIndices = (int *) malloc(nsosnz*sizeof(int));
	info->sosReference = (double *) malloc(nsosnz*sizeof(double));
	sos_kludge(nsos, sosbeg, sosref,sosind);
	for (int i=0;i<nsos;i++) {
	  int ichar = sostype[i];
	  assert (ichar=='1'||ichar=='2');
	  info->sosType[i]=ichar-'0';
	}	
	memcpy(info->sosPriority,sospri,nsos*sizeof(int));
	memcpy(info->sosStart,sosbeg,(nsos+1)*sizeof(int));
	memcpy(info->sosIndices,sosind,nsosnz*sizeof(int));
	memcpy(info->sosReference,sosref,nsosnz*sizeof(double));
      }
    }
    
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
    if (niv+nbv>0)
      mip_stuff(); // get any extra info
    if ((!(niv+nbv)&&(csd->kind & ASL_Sufkind_input))
	||(rsd->kind & ASL_Sufkind_input)) {
      /* convert status - need info on map */
      static int map[] = {1, 3, 1, 1, 2, 1, 1};
      stat_map(info->columnStatus, n_var, map, 6, "incoming columnStatus");
      stat_map(info->rowStatus, n_con, map, 6, "incoming rowStatus");
    } else {
      /* all slack basis */
      // leave status for output */
#if 0
      free(info->rowStatus);
      info->rowStatus=NULL;
      free(info->columnStatus);
      info->columnStatus=NULL;
#endif
    }
  } else {
    // QP
    CoinModel * model = new CoinModel(1,fileName);
    if (model->numberRows()>0)
      *coinModel=(void *) model;
    Oinfo.uinfo = tempBuffer;
    if (getopts(argv, &Oinfo))
      return 1;
    if (objtype[0])
      info->direction=-1.0;
    else
      info->direction=1.0;
    info->offset=objconst(0);
    info->numberRows=n_con;
    info->numberColumns=n_var;
    info->numberElements=nzc;;
    info->numberBinary=nbv;
    info->numberIntegers=niv;
  }
  /* add -solve - unless something there already
   - also check for sleep=yes */
  {
    int found=0;
    int foundLog=0;
    int foundSleep=0;
    const char * something[]={"solve","branch","duals","primals","user"};
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
  free(info->priorities);
  info->priorities=NULL;
  free(info->branchDirection);
  info->branchDirection=NULL;
  free(info->pseudoDown);
  info->pseudoDown=NULL;
  free(info->pseudoUp);
  info->pseudoUp=NULL;
  free(info->sosType);
  info->sosType=NULL;
  free(info->sosPriority);
  info->sosPriority=NULL;
  free(info->sosStart);
  info->sosStart=NULL;
  free(info->sosIndices);
  info->sosIndices=NULL;
  free(info->sosReference);
  info->sosReference=NULL;
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
/* Read a problem from AMPL nl file
 */
CoinModel::CoinModel( int nonLinear, const char * fileName)
 :  numberRows_(0),
    maximumRows_(0),
    numberColumns_(0),
    maximumColumns_(0),
    numberElements_(0),
    maximumElements_(0),
    numberQuadraticElements_(0),
    maximumQuadraticElements_(0),
    optimizationDirection_(1.0),
    objectiveOffset_(0.0),
    rowLower_(NULL),
    rowUpper_(NULL),
    rowType_(NULL),
    objective_(NULL),
    columnLower_(NULL),
    columnUpper_(NULL),
    integerType_(NULL),
    columnType_(NULL),
    start_(NULL),
    elements_(NULL),
    quadraticElements_(NULL),
    sortIndices_(NULL),
    sortElements_(NULL),
    sortSize_(0),
    sizeAssociated_(0),
    associated_(NULL),
    numberSOS_(0),
    startSOS_(NULL),
    memberSOS_(NULL),
    typeSOS_(NULL),
    prioritySOS_(NULL),
    referenceSOS_(NULL),
    priority_(NULL),
    logLevel_(0),
    type_(-1),
    links_(0)
{
  problemName_ = strdup("");
  int status=0;
  if (!strcmp(fileName,"-")||!strcmp(fileName,"stdin")) {
    // stdin
  } else {
    std::string name=fileName;
    bool readable = fileCoinReadable(name);
    if (!readable) {
      std::cerr<<"Unable to open file "
	       <<fileName<<std::endl;
      status = -1;
    }
  }
  if (!status) {
    gdb(nonLinear,fileName);
  }
}
 static real
qterm(ASL *asl, fint *colq, fint *rowq, real *delsq)
{
  double t, t1, *x, *x0, *xe;
  fint *rq0, *rqe;
  
  t = 0.;
  x0 = x = X0;
  xe = x + n_var;
  rq0 = rowq;
  while(x < xe) {
    t1 = *x++;
    rqe = rq0 + *++colq;
    while(rowq < rqe)
      t += t1*x0[*rowq++]**delsq++;
  }
  return 0.5 * t;
}

void
CoinModel::gdb( int nonLinear, const char * fileName)
{
  ograd *og=NULL;
  int i;
  SufDesc *csd=NULL;
  SufDesc *rsd=NULL;
  /*bool *basis, *lower;*/
  /*double *LU, *c, lb, objadj, *rshift, *shift, t, ub, *x, *x0, *x1;*/
  //char tempBuffer[20];
  double * objective=NULL;
  double * columnLower=NULL;
  double * columnUpper=NULL;
  double * rowLower=NULL;
  double * rowUpper=NULL;
  int * columnStatus=NULL;
  int * rowStatus=NULL;
  int numberRows=-1;
  int numberColumns=-1;
  int numberElements=-1;
  int numberBinary=-1;
  int numberIntegers=-1;
  double * primalSolution=NULL;
  double direction=1.0;
  char * stub = strdup(fileName);
  CoinPackedMatrix matrixByRow;
  fint ** colqp=NULL;
  int *z = NULL;
  if (nonLinear==0) {
    // linear
    asl = ASL_alloc(ASL_read_f);
    nl = jac0dim(stub, 0);
    free(stub);
    suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
    
    /* set A_vals to get the constraints column-wise (malloc so can be freed) */
    A_vals = (double *) malloc(nzc*sizeof(double));
    if (!A_vals) {
      printf("no memory\n");
      return ;
    }
    /* say we want primal solution */
    want_xpi0=1;
    /* for basis info */
    columnStatus = (int *) malloc(n_var*sizeof(int));
    rowStatus = (int *) malloc(n_con*sizeof(int));
    csd = suf_iput("sstatus", ASL_Sufkind_var, columnStatus);
    rsd = suf_iput("sstatus", ASL_Sufkind_con, rowStatus);
    /* read linear model*/
    f_read(nl,0);
    // see if any sos
    if (true) {
      char *sostype;
      int nsosnz, *sosbeg, *sosind, * sospri;
      double *sosref;
      int nsos;
      int i = ASL_suf_sos_explict_free;
      int copri[2], **p_sospri;
      copri[0] = 0;
      copri[1] = 0;
      p_sospri = &sospri;
      nsos = suf_sos(i, &nsosnz, &sostype, p_sospri, copri,
		     &sosbeg, &sosind, &sosref);
      if (nsos) {
	abort();
#if 0
	info->numberSos=nsos;
	info->sosType = (char *) malloc(nsos);
	info->sosPriority = (int *) malloc(nsos*sizeof(int));
	info->sosStart = (int *) malloc((nsos+1)*sizeof(int));
	info->sosIndices = (int *) malloc(nsosnz*sizeof(int));
	info->sosReference = (double *) malloc(nsosnz*sizeof(double));
	sos_kludge(nsos, sosbeg, sosref,sosind);
	for (int i=0;i<nsos;i++) {
	  int ichar = sostype[i];
	  assert (ichar=='1'||ichar=='2');
	  info->sosType[i]=ichar-'0';
	}	
	memcpy(info->sosPriority,sospri,nsos*sizeof(int));
	memcpy(info->sosStart,sosbeg,(nsos+1)*sizeof(int));
	memcpy(info->sosIndices,sosind,nsosnz*sizeof(int));
	memcpy(info->sosReference,sosref,nsosnz*sizeof(double));
#endif
      }
    }
    
    /*sos_finish(&specialOrderedInfo, 0, &j, 0, 0, 0, 0, 0);*/
    //Oinfo.uinfo = tempBuffer;
    //if (getopts(argv, &Oinfo))
    //return 1;
    /* objective*/
    objective = (double *) malloc(n_var*sizeof(double));
    for (i=0;i<n_var;i++)
      objective[i]=0.0;;
    if (n_obj) {
      for (og = Ograd[0];og;og = og->next)
	objective[og->varno] = og->coef;
    }
    if (objtype[0])
      direction=-1.0;
    else
      direction=1.0;
    objectiveOffset_=objconst(0);
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
    numberRows=n_con;
    numberColumns=n_var;
    numberElements=nzc;;
    numberBinary=nbv;
    numberIntegers=niv;
    /* put in primalSolution if exists */
    if (X0) {
      primalSolution=(double *) malloc(n_var*sizeof(double));
      memcpy( primalSolution,X0,n_var*sizeof(double));
    }
    //double * dualSolution=NULL;
    if (niv+nbv>0)
      mip_stuff(); // get any extra info
    if ((!(niv+nbv)&&(csd->kind & ASL_Sufkind_input))
	||(rsd->kind & ASL_Sufkind_input)) {
      /* convert status - need info on map */
      static int map[] = {1, 3, 1, 1, 2, 1, 1};
      stat_map(columnStatus, n_var, map, 6, "incoming columnStatus");
      stat_map(rowStatus, n_con, map, 6, "incoming rowStatus");
    } else {
      /* all slack basis */
      // leave status for output */
#if 0
      free(rowStatus);
      rowStatus=NULL;
      free(columnStatus);
      columnStatus=NULL;
#endif
    }
    CoinPackedMatrix columnCopy(true,numberRows,numberColumns,numberElements,
				A_vals,A_rownos,A_colstarts,NULL);
    matrixByRow.reverseOrderedCopyOf(columnCopy);
  } else if (nonLinear==1) {
    // quadratic
    asl = ASL_alloc(ASL_read_fg);
    nl = jac0dim(stub, (long) strlen(stub));
    free(stub);
    suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
    /* read  model*/
    X0 = (double*) malloc(n_var*sizeof(double));
    CoinZeroN(X0,n_var);
    qp_read(nl,0);
    assert (n_obj==1);
    int nz = 1 + n_con;
    colqp = (fint**) malloc(nz*(2*sizeof(int*)
				      + sizeof(double*)));
    fint ** rowqp = colqp + nz;
    double ** delsqp = (double **)(rowqp + nz);
    z = (int*) malloc(nz*sizeof(int));
    for (i=0;i<=n_con;i++) {
      z[i] = nqpcheck(-i, rowqp+i, colqp+i, delsqp+i);
    }
    qp_opify();
    /* objective*/
    objective = (double *) malloc(n_var*sizeof(double));
    for (i=0;i<n_var;i++)
      objective[i]=0.0;;
    if (n_obj) {
      for (og = Ograd[0];og;og = og->next)
	objective[og->varno] = og->coef;
    }
    if (objtype[0])
      direction=-1.0;
    else
      direction=1.0;
    objectiveOffset_=objconst(0);
    /* Column bounds*/
    columnLower = (double *) malloc(n_var*sizeof(double));
    columnUpper = (double *) malloc(n_var*sizeof(double));
    for (i=0;i<n_var;i++) {
      columnLower[i]=LUv[2*i];
      if (columnLower[i]<= negInfinity)
	columnLower[i]=-COIN_DBL_MAX;
      columnUpper[i]=LUv[2*i+1];
      if (columnUpper[i]>= Infinity)
	columnUpper[i]=COIN_DBL_MAX;
    }
    // Build by row from scratch
    //matrixByRow.reserve(n_var,nzc,true);
    // say row orderded
    matrixByRow.transpose();
    /* Row bounds*/
    rowLower = (double *) malloc(n_con*sizeof(double));
    rowUpper = (double *) malloc(n_con*sizeof(double));
    int * column = new int [n_var];
    double * element = new double [n_var];
    for (i=0;i<n_con;i++) {
      rowLower[i]=LUrhs[2*i];
      if (rowLower[i]<= negInfinity)
	rowLower[i]=-COIN_DBL_MAX;
      rowUpper[i]=LUrhs[2*i+1];
      if (rowUpper[i]>= Infinity)
	rowUpper[i]=COIN_DBL_MAX;
      int k=0;
      for(cgrad * cg = Cgrad[i]; cg; cg = cg->next) {
	column[k]=cg->varno;
	element[k++]=cg->coef;
      }
      matrixByRow.appendRow(k,column,element);
    }
    numberRows=n_con;
    numberColumns=n_var;
    numberElements=nzc;;
    numberBinary=nbv;
    numberIntegers=niv+nlvci;
    //double * dualSolution=NULL;
  } else {
    abort();
  }
  // set problem name
  free(problemName_);
  problemName_=strdup("???");

  // Build by row from scratch
  const double * element = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();
  for (i=0;i<numberRows;i++) {
    addRow(rowLength[i],column+rowStart[i],
	   element+rowStart[i],rowLower[i],rowUpper[i]);
  }
  // Now do column part
  for (i=0;i<numberColumns;i++) {
    setColumnBounds(i,columnLower[i],columnUpper[i]);
    setColumnObjective(i,objective[i]);
  }
  for ( i=numberColumns-numberBinary-numberIntegers;
	i<numberColumns;i++) {
      setColumnIsInteger(i,true);;
  }
  // do names
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    char name[9];
    sprintf(name,"r%7.7d",iRow);
    setRowName(iRow,name);
  }
    
  int iColumn;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    char name[9];
    sprintf(name,"c%7.7d",iColumn);
    setColumnName(iColumn,name);
  }
  if (colqp) {
    // add in quadratic
    int nz = 1 + n_con;
    int nOdd=0;
    fint ** rowqp = colqp + nz;
    double ** delsqp = (double **)(rowqp + nz);
    for (i=0;i<=n_con;i++) {
      int nels = z[i];
      if (nels) {
	double * element = delsqp[i];
	int * start = (int *) colqp[i];
	int * row = (int *) rowqp[i];
	if (!element) {
	  // odd row - probably not quadratic
	  nOdd++;
	  continue;
	}
#if 0
	printf("%d quadratic els\n",nels);
	for (int j=0;j<n_var;j++) {
	  for (int k=start[j];k<start[j+1];k++)
	    printf("%d %d %g\n",j,row[k],element[k]);
	}
#endif
	if (i) {
	  int iRow = i-1;
	  for (int j=0;j<n_var;j++) {
	    for (int k=start[j];k<start[j+1];k++) {
	      int kColumn = row[k];
	      double value = element[k];
	      // ampl gives twice with assumed 0.5
	      if (kColumn>j)
		continue;
	      else if (kColumn==j)
		value *= 0.5;
	      const char * expr = getElementAsString(iRow,j);
	      double constant=0.0;
	      bool linear;
	      if (expr&&strcmp(expr,"Numeric")) {
		linear=false;
	      } else {
		constant = getElement(iRow,j);
		linear=true;
	      }
	      char temp[1000];
	      char temp2[30];
	      if (value==1.0) 
		sprintf(temp2,"c%7.7d",kColumn);
	      else
		sprintf(temp2,"%g*c%7.7d",value,kColumn);
	      if (linear) {
		if (!constant)
		  strcpy(temp,temp2);
		else if (value>0.0) 
		  sprintf(temp,"%g+%s",constant,temp2);
		else
		  sprintf(temp,"%g%s",constant,temp2);
	      } else {
		if (value>0.0) 
		  sprintf(temp,"%s+%s",expr,temp2);
		else
		  sprintf(temp,"%s%s",expr,temp2);
	      }
	      assert (strlen(temp)<1000);
	      setElement(iRow,j,temp);
	      printf("el for row %d column c%7.7d is %s\n",iRow,j,temp);
	    }
	  }
	} else {
	  // objective
	  for (int j=0;j<n_var;j++) {
	    for (int k=start[j];k<start[j+1];k++) {
	      int kColumn = row[k];
	      double value = element[k];
	      // ampl gives twice with assumed 0.5
	      if (kColumn>j)
		continue;
	      else if (kColumn==j)
		value *= 0.5;
	      const char * expr = getColumnObjectiveAsString(j);
	      double constant=0.0;
	      bool linear;
	      if (expr&&strcmp(expr,"Numeric")) {
		linear=false;
	      } else {
		constant = getColumnObjective(j);
		linear=true;
	      }
	      char temp[1000];
	      char temp2[30];
	      if (value==1.0) 
		sprintf(temp2,"c%7.7d",kColumn);
	      else
		sprintf(temp2,"%g*c%7.7d",value,kColumn);
	      if (linear) {
		if (!constant)
		  strcpy(temp,temp2);
		else if (value>0.0) 
		  sprintf(temp,"%g+%s",constant,temp2);
		else
		  sprintf(temp,"%g%s",constant,temp2);
	      } else {
		if (value>0.0) 
		  sprintf(temp,"%s+%s",expr,temp2);
		else
		  sprintf(temp,"%s%s",expr,temp2);
	      }
	      assert (strlen(temp)<1000);
	      setObjective(j,temp);
	      printf("el for objective column c%7.7d is %s\n",j,temp);
	    }
	  }
	}
      }
    }
    if (nOdd) {
      printf("%d non-linear constraints could not be converted to quadratic\n",nOdd);
      exit(77);
    }
  }
  // see if any sos
  {
    char *sostype;
    int nsosnz, *sosbeg, *sosind, * sospri;
    double *sosref;
    int nsos;
    int i = ASL_suf_sos_explict_free;
    int copri[2], **p_sospri;
    copri[0] = 0;
    copri[1] = 0;
    p_sospri = &sospri;
    nsos = suf_sos(i, &nsosnz, &sostype, p_sospri, copri,
		   &sosbeg, &sosind, &sosref);
    if (nsos) {
      numberSOS_=nsos;
      typeSOS_ = new int [numberSOS_];
      prioritySOS_ = new int [numberSOS_];
      startSOS_ = new int [numberSOS_+1];
      memberSOS_ = new int[nsosnz];
      referenceSOS_ = new double [nsosnz];
      sos_kludge(nsos, sosbeg, sosref,sosind);
      for (int i=0;i<nsos;i++) {
	int ichar = sostype[i];
	assert (ichar=='1'||ichar=='2');
	typeSOS_[i]=ichar-'0';
      }	
      memcpy(prioritySOS_,sospri,nsos*sizeof(int));
      memcpy(startSOS_,sosbeg,(nsos+1)*sizeof(int));
      memcpy(memberSOS_,sosind,nsosnz*sizeof(int));
      memcpy(referenceSOS_,sosref,nsosnz*sizeof(double));
    }
  }
}
