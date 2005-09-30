// copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <typeinfo>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>

#define CBCVERSION "0.97"

#include "CoinMpsIO.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

#include "CglCutGenerator.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"

#include "CbcModel.hpp"
#include "CbcTree.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristic.hpp"
#include "CbcCompareActual.hpp"
#include  "CbcParam.hpp"

#ifdef COIN_USE_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#ifdef COIN_USE_DYLP
#include "OsiDylpSolverInterface.hpp"
#endif




/* Before first solution do depth first,
   then it is computed to hit first solution less 2%
*/
class CbcCompareUser  : public CbcCompareBase {
public:
  // Weight for each infeasibility
  double weight_;
  // Weight for each infeasibility - computed from solution
  double saveWeight_;
  // Number of solutions
  int numberSolutions_;
  // Tree size (at last check)
  int treeSize_;
  // Default Constructor 
  CbcCompareUser () : CbcCompareBase(),
                      weight_(-1.0),saveWeight_(0.0),numberSolutions_(0),
                      treeSize_(0)
  { test_=this;};

  // Copy constructor 
  CbcCompareUser ( const CbcCompareUser &rhs)
    : CbcCompareBase(rhs)
  {
    weight_=rhs.weight_;
    saveWeight_ = rhs.saveWeight_;
    numberSolutions_=rhs.numberSolutions_;
    treeSize_ = rhs.treeSize_;
  };
   
  // Assignment operator 
  CbcCompareUser & operator=( const CbcCompareUser& rhs)
  {
    if (this!=&rhs) { 
      CbcCompareBase::operator=(rhs);
      weight_=rhs.weight_;
      saveWeight_ = rhs.saveWeight_;
      numberSolutions_=rhs.numberSolutions_;
      treeSize_ = rhs.treeSize_;
    }
    return *this;
  };

  /// Clone
  virtual CbcCompareBase * clone() const
  { 
    return new CbcCompareUser (*this);
  };

  ~CbcCompareUser() {};

  /* 
     Return true if y better than x
     Node y is better than node x if y has fewer unsatisfied (greater depth on tie) or
     after solution weighted value of y is less than weighted value of x
  */
  virtual bool test (CbcNode * x, CbcNode * y) {
    if (weight_==-1.0) {
      // before solution
      if (x->numberUnsatisfied() > y->numberUnsatisfied())
        return true;
      else if (x->numberUnsatisfied() < y->numberUnsatisfied())
        return false;
      else
        return x->depth() < y->depth();
    } else {
      // after solution
      double weight = CoinMax(weight_,0.0);
      return x->objectiveValue()+ weight*x->numberUnsatisfied() > 
        y->objectiveValue() + weight*y->numberUnsatisfied();
    }
  }
  // This allows method to change behavior as it is called
  // after each solution
  virtual void newSolution(CbcModel * model,
			   double objectiveAtContinuous,
			   int numberInfeasibilitiesAtContinuous) 
  {
    if (model->getSolutionCount()==model->getNumberHeuristicSolutions())
      return; // solution was got by rounding
    // set to get close to this solution
    double costPerInteger = 
      (model->getObjValue()-objectiveAtContinuous)/
      ((double) numberInfeasibilitiesAtContinuous);
    weight_ = 0.98*costPerInteger;
    saveWeight_=weight_;
    numberSolutions_++;
    if (numberSolutions_>5)
      weight_ =0.0; // this searches on objective
  }
  // This allows method to change behavior 
  virtual bool every1000Nodes(CbcModel * model, int numberNodes)
  {
    if (numberNodes>10000)
      weight_ =0.0; // this searches on objective
    else if (numberNodes==1000&&weight_==-2.0)
      weight_=-1.0; // Go to depth first
    // get size of tree
    treeSize_ = model->tree()->size();
    if (treeSize_>10000) {
      // set weight to reduce size most of time
      if (treeSize_>20000)
        weight_=-1.0;
      else if ((numberNodes%4000)!=0)
        weight_=-1.0;
      else
        weight_=saveWeight_;
    }
    return numberNodes==11000; // resort if first time
  }
};


#define MAXPARAMETERS 100

namespace {

void establishParams (int &numberParameters, CbcParam *const parameters)
/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcParam for details. Pulled from main() for clarity. 
*/
{ numberParameters = 0 ;

  parameters[numberParameters++]=
      CbcParam("?","For help",GENERALQUERY);
  parameters[numberParameters++]=
      CbcParam("dualT!olerance",
	       "For an optimal solution no dual infeasibility may "
	       "exceed this value",
	       1.0e-20,1.0e12,DUALTOLERANCE);
  parameters[numberParameters++]=
      CbcParam("primalT!olerance",
	       "For an optimal solution no primal infeasibility may "
	       " exceed this value",
	      1.0e-20,1.0e12,PRIMALTOLERANCE);
    parameters[numberParameters++]=
      CbcParam("inf!easibilityWeight","Each integer infeasibility is expected \
to cost this much",
	      0.0,1.0e20,INFEASIBILITYWEIGHT);
    parameters[numberParameters++]=
      CbcParam("integerT!olerance","For an optimal solution \
no integer variable may be this away from an integer value",
	      1.0e-20,0.5,INTEGERTOLERANCE);
    parameters[numberParameters++]=
      CbcParam("inc!rement","A valid solution must be at least this \
much better than last integer solution",
	      -1.0e20,1.0e20,INCREMENT);
    parameters[numberParameters++]=
      CbcParam("allow!ableGap","Stop when gap between best possible and \
best less than this",
	      0.0,1.0e20,ALLOWABLEGAP);
    parameters[numberParameters++]=
      CbcParam("ratio!Gap","Stop when gap between best possible and \
best less than this fraction of larger of two",
	      0.0,1.0e20,GAPRATIO);
    parameters[numberParameters++]=
      CbcParam("fix!OnDj","Try heuristic based on fixing variables with \
reduced costs greater than this",
	      -1.0e20,1.0e20,DJFIX);
    parameters[numberParameters++]=
      CbcParam("tighten!Factor","Tighten bounds using this times largest \
activity at continuous solution",
	      1.0,1.0e20,TIGHTENFACTOR);
    parameters[numberParameters++]=
      CbcParam("log!Level","Level of detail in BAB output",
	      0,63,LOGLEVEL);
    parameters[numberParameters++]=
      CbcParam("slog!Level","Level of detail in Solver output",
	      0,63,SOLVERLOGLEVEL);
    parameters[numberParameters++]=
      CbcParam("maxN!odes","Maximum number of nodes to do",
	      1,999999,MAXNODES);
    parameters[numberParameters++]=
      CbcParam("strong!Branching","Number of variables to look at in strong branching",
	      0,999999,STRONGBRANCHING);
    parameters[numberParameters++]=
      CbcParam("direction","Minimize or Maximize",
	      "min!imize",DIRECTION);
    parameters[numberParameters-1].append("max!imize");
    parameters[numberParameters++]=
      CbcParam("error!sAllowed","Whether to allow import errors",
	      "off",ERRORSALLOWED);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters++]=
      CbcParam("gomory!Cuts","Whether to use Gomory cuts",
	      "off",GOMORYCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("probing!Cuts","Whether to use Probing cuts",
	      "off",PROBINGCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("knapsack!Cuts","Whether to use Knapsack cuts",
	      "off",KNAPSACKCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("oddhole!Cuts","Whether to use Oddhole cuts",
	      "off",ODDHOLECUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("clique!Cuts","Whether to use Clique cuts",
	      "off",CLIQUECUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("mixed!IntegerRoundingCuts","Whether to use Mixed Integer Rounding cuts",
	      "off",MIXEDCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("flow!CoverCuts","Whether to use Flow Cover cuts",
	      "off",FLOWCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("two!MirCuts","Whether to use Two phase Mixed Integer Rounding cuts",
	      "off",TWOMIRCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters++]=
      CbcParam("round!ingHeuristic","Whether to use Rounding heuristic",
	      "off",ROUNDING);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters++]=
      CbcParam("cost!Strategy","How to use costs",
	      "off",COSTSTRATEGY);
    parameters[numberParameters-1].append("pri!orities");
    parameters[numberParameters-1].append("pseudo!costs(not implemented yet)");
    parameters[numberParameters++]=
      CbcParam("keepN!ames","Whether to keep names from import",
	      "on",KEEPNAMES);
    parameters[numberParameters-1].append("off");
    parameters[numberParameters++]=
      CbcParam("scaling","Whether to do scaling",
	      "on",SCALING);
    parameters[numberParameters-1].append("off");
    parameters[numberParameters++]=
      CbcParam("directory","Set Default import directory",
	      DIRECTORY);
    parameters[numberParameters++]=
      CbcParam("solver!","Set the solver used by cbc",
	       SOLVER) ;
    parameters[numberParameters++]=
      CbcParam("import","Import model from mps file",
	      IMPORT);
    parameters[numberParameters++]=
      CbcParam("export","Export model as mps file",
	      EXPORT);
    parameters[numberParameters++]=
      CbcParam("save!Model","Save model to binary file",
	      SAVE);
    parameters[numberParameters++]=
      CbcParam("restore!Model","Restore model from binary file",
	      RESTORE);
    parameters[numberParameters++]=
      CbcParam("presolve","Whether to use integer presolve - be careful",
	      "off",PRESOLVE);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters++]=
      CbcParam("preprocess","Whether to use integer preprocessing",
	      "off",PREPROCESS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters++]=
      CbcParam("initialS!olve","Solve to continuous",
	      SOLVECONTINUOUS);
    parameters[numberParameters++]=
      CbcParam("branch!AndBound","Do Branch and Bound",
	      BAB);
    parameters[numberParameters++]=
      CbcParam("sol!ution","Prints solution to file",
	      SOLUTION);
    parameters[numberParameters++]=
      CbcParam("max!imize","Set optimization direction to maximize",
	      MAXIMIZE);
    parameters[numberParameters++]=
      CbcParam("min!imize","Set optimization direction to minimize",
	      MINIMIZE);
    parameters[numberParameters++] =
      CbcParam("time!Limit","Set a time limit for solving this problem",
	      1.0,(double)(60*60*24*365*10),TIMELIMIT) ;
    parameters[numberParameters++]=
      CbcParam("exit","Stops cbc execution",
	      EXIT);
    parameters[numberParameters++]=
      CbcParam("stop","Stops cbc execution",
	      EXIT);
    parameters[numberParameters++]=
      CbcParam("quit","Stops cbc execution",
	      EXIT);
    parameters[numberParameters++]=
      CbcParam("-","From stdin",
	      STDIN);
    parameters[numberParameters++]=
      CbcParam("stdin","From stdin",
	      STDIN);
    parameters[numberParameters++]=
      CbcParam("unitTest","Do unit test",
	      UNITTEST);
    parameters[numberParameters++]=
      CbcParam("miplib","Do some of miplib test set",
	      MIPLIB);
    parameters[numberParameters++]=
      CbcParam("ver!sion","Print out version",
	      VERSION);
    parameters[numberParameters++]=
      CbcParam("alg!orithm","Whether to use dual or primal",
	      "dual",ALGORITHM);
    parameters[numberParameters-1].append("primal");

    assert(numberParameters<MAXPARAMETERS);

  return ; }

#ifdef COIN_USE_READLINE     
#include <readline/readline.h>
#include <readline/history.h>
#endif

// Returns next valid field

int read_mode=1;
char line[1000];
char * where=NULL;

std::string
nextField()
{
  std::string field;
  if (!where) {
    // need new line
#ifdef COIN_USE_READLINE     
    // Get a line from the user. 
    where = readline ("Cbc:");
     
    // If the line has any text in it, save it on the history.
    if (where) {
      if ( *where)
	add_history (where);
      strcpy(line,where);
    }
#else
    fprintf(stdout,"Cbc:");
    fflush(stdout);
    where = fgets(line,1000,stdin);
#endif
    if (!where)
      return field; // EOF
    where = line;
    // clean image
    char * lastNonBlank = line-1;
    while ( *where != '\0' ) {
      if ( *where != '\t' && *where < ' ' ) {
	break;
      } else if ( *where != '\t' && *where != ' ') {
	lastNonBlank = where;
      }
      where++;
    }
    where=line;
    *(lastNonBlank+1)='\0';
  }
  // munch white space
  while(*where==' '||*where=='\t')
    where++;
  char * saveWhere = where;
  while (*where!=' '&&*where!='\t'&&*where!='\0')
    where++;
  if (where!=saveWhere) {
    char save = *where;
    *where='\0';
    //convert to string
    field=saveWhere;
    *where=save;
  } else {
    where=NULL;
    field="EOL";
  }
  return field;
}

std::string
getCommand(int argc, const char *argv[])
{
  std::string field="EOL";
  while (field=="EOL") {
    if (read_mode>0) {
      if (read_mode<argc) {
	field = argv[read_mode++];
	if (field=="-") {
	  std::cout<<"Switching to line mode"<<std::endl;
	  read_mode=-1;
	  field=nextField();
	} else if (field[0]!='-') {
	  if (read_mode!=2) {
	    std::cout<<"skipping non-command "<<field<<std::endl;
	    field="EOL"; // skip
	  } else {
	    // special dispensation - taken as -import name
	    read_mode--;
	    field="import";
	  }
	} else {
	  if (field!="--") {
	    // take off -
	    field = field.substr(1);
	  } else {
	    // special dispensation - taken as -import --
	    read_mode--;
	    field="import";
	  }
	}
      } else {
	field="";
      }
    } else {
      field=nextField();
    }
  }
  //std::cout<<field<<std::endl;
  return field;
}
std::string
getString(int argc, const char *argv[])
{
  std::string field="EOL";
  if (read_mode>0) {
    if (read_mode<argc) {
      if (argv[read_mode][0]!='-') { 
	field = argv[read_mode++];
      } else if (!strcmp(argv[read_mode],"--")) {
	field = argv[read_mode++];
	// -- means import from stdin
	field = "-";
      }
    }
  } else {
    field=nextField();
  }
  //std::cout<<field<<std::endl;
  return field;
}

// valid = 0 okay, 1 bad, 2 not there
int
getIntField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (read_mode>0) {
    if (read_mode<argc) {
      // may be negative value so do not check for -
      field = argv[read_mode++];
    }
  } else {
    field=nextField();
  }
  int value=0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    // how do I check valid
    value =  atoi(field.c_str());
    *valid=0;
  } else {
    *valid=2;
  }
  return value;
}


// valid = 0 okay, 1 bad, 2 not there
double
getDoubleField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (read_mode>0) {
    if (read_mode<argc) {
      // may be negative value so do not check for -
      field = argv[read_mode++];
    }
  } else {
    field=nextField();
  }
  double value=0.0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    // how do I check valid
    value = atof(field.c_str());
    *valid=0;
  } else {
    *valid=2;
  }
  return value;
}

/// For run timing

double totalTime=0.0;

}	/* end unnamed namespace */


int main (int argc, const char *argv[])
{
  // next {} is just to make sure all memory should be freed - for debug
  {
    std::ios::sync_with_stdio() ;
/*
  Create a vector of solver prototypes and establish a default solver. After
  this the code is solver independent.

  Creating multiple solvers is moderately expensive. If you don't want to
  make use of this feature, best to define just one. The businesss with
  CBC_DEFAULT_SOLVER will select the first available solver as the default,
  unless overridden at compile time.

  NOTE that processing of string parameters is case-independent, but maps are
       case-sensitive. The solver name given here must contain only lower case
       letters.
*/
    typedef std::map<std::string,OsiSolverInterface*> solverMap_t ;
    typedef solverMap_t::const_iterator solverMapIter_t ;

    solverMap_t solvers ;

#   ifdef COIN_USE_CLP
#     ifndef CBC_DEFAULT_SOLVER
#       define CBC_DEFAULT_SOLVER "clp"
#     endif
      solvers["clp"] = new OsiClpSolverInterface ;
#   else
      solvers["clp"] = 0 ;
#   endif
#   ifdef COIN_USE_DYLP
#     ifndef CBC_DEFAULT_SOLVER
#       define CBC_DEFAULT_SOLVER "dylp"
#     endif
      solvers["dylp"] = new OsiDylpSolverInterface  ;
#   else
      solvers["dylp"] = 0 ;
#   endif
/*
  If we don't have a default solver, we're deeply confused.
*/
    OsiSolverInterface *dflt_solver = solvers[CBC_DEFAULT_SOLVER] ;
    if (dflt_solver)
    { std::cout << "Default solver is " << CBC_DEFAULT_SOLVER << std::endl ; }
    else
    { std::cerr << "No solvers! Aborting." << std::endl ;
      return (1) ; }
/*
  For calculating run time.
*/
    double time1 = CoinCpuTime() ;
    double time2;
/*
  Establish the command line interface: parameters, with associated info
  messages, ranges, defaults. See CbcParam for details. Scan the vector of
  solvers and add the names to the parameter object.
*/
    CbcParam parameters[MAXPARAMETERS];
    int numberParameters ;
    establishParams(numberParameters,parameters) ;

    { int iSolver = 0 ;
      for ( ; iSolver < numberParameters ; iSolver++)
      { int match = parameters[iSolver].matches("solver") ;
	if (match==1) break ; }
      for (solverMapIter_t solverIter = solvers.begin() ;
	   solverIter != solvers.end() ;
	   solverIter++)
      { if (solverIter->second)
	  parameters[iSolver].append(solverIter->first) ; }
      int iKwd = parameters[iSolver].parameterOption(CBC_DEFAULT_SOLVER) ;
      parameters[iSolver].setCurrentOption(iKwd) ; }
/*
  The rest of the default setup: establish a model, instantiate cut generators
  and heuristics, set various default values.
*/
    dflt_solver->messageHandler()->setLogLevel(0) ;
    CbcModel *model = new CbcModel(*dflt_solver) ;
    model->messageHandler()->setLogLevel(1);
    bool goodModel=false;
    
// Set up likely cut generators and defaults

    CglGomory gomoryGen;
    // try larger limit
    gomoryGen.setLimit(300);
    // set default action (0=off,1=on,2=root)
    int gomoryAction=1;

    CglProbing probingGen;
    probingGen.setUsingObjective(true);
    probingGen.setMaxPass(3);
    probingGen.setMaxPassRoot(3);
    // Number of unsatisfied variables to look at
    probingGen.setMaxProbe(10);
    probingGen.setMaxProbeRoot(50);
    // How far to follow the consequences
    probingGen.setMaxLook(10);
    probingGen.setMaxLookRoot(50);
    // Only look at rows with fewer than this number of elements
    probingGen.setMaxElements(200);
    probingGen.setRowCuts(3);
    // set default action (0=off,1=on,2=root)
    int probingAction=1;

    CglKnapsackCover knapsackGen;
    // set default action (0=off,1=on,2=root)
    int knapsackAction=1;

    CglOddHole oddholeGen;
    oddholeGen.setMinimumViolation(0.005);
    oddholeGen.setMinimumViolationPer(0.0002);
    oddholeGen.setMaximumEntries(100);
    // set default action (0=off,1=on,2=root)
    int oddholeAction=0;

    CglClique cliqueGen;
    cliqueGen.setStarCliqueReport(false);
    cliqueGen.setRowCliqueReport(false);
    // set default action (0=off,1=on,2=root)
    int cliqueAction=1;

    CglMixedIntegerRounding mixedGen;
    // set default action (0=off,1=on,2=root)
    int mixedAction=1;

    CglFlowCover flowGen;
    // set default action (0=off,1=on,2=root)
    int flowAction=1;

    CglTwomir twomirGen;
    // set default action (0=off,1=on,2=root)
    int twomirAction=0;

    bool useRounding=true;
   
    int allowImportErrors=0;
    int algorithm=0; // dual
    int keepImportNames=1;	// not implemented
    int doScaling=1;
    int preSolve=0;
    int preProcess=1;
    double djFix=1.0e100;
    double gapRatio=1.0e100;
    double tightenFactor=0.0;

    // Set false if user does anything advanced
    bool defaultSettings=true;

    std::string directory ="./";
    std::string field;
/*
  The main command parsing loop.
*/
    // total number of commands read
    int numberGoodCommands=0;
    
    while (1) {
      // next command
      field=getCommand(argc,argv);
      
      // exit if null or similar
      if (!field.length()) {
	if (numberGoodCommands==1&&goodModel) {
	  // we just had file name
          field="branchAndBound";
	} else if (!numberGoodCommands) {
	  // let's give the sucker a hint
	  std::cout
	    <<"Cbc takes input from arguments ( - switches to stdin)"
	    <<std::endl
	    <<"Enter ? for list of commands or (-)unitTest or -miplib"
	    <<" for tests"<<std::endl;
          break;
	} else {
          break;
        }
      }
      
      // see if ? at end
      int numberQuery=0;
      if (field!="?") {
	int length = field.length();
	int i;
	for (i=length-1;i>0;i--) {
	  if (field[i]=='?') 
	    numberQuery++;
	  else
	    break;
	}
	field=field.substr(0,length-numberQuery);
      }
      // find out if valid command
      int iParam;
      int numberMatches=0;
      for ( iParam=0; iParam<numberParameters; iParam++ ) {
	int match = parameters[iParam].matches(field);
	if (match==1) {
	  numberMatches = 1;
	  break;
	} else {
	  numberMatches += match>>1;
	}
      }
      if (iParam<numberParameters&&!numberQuery) {
	// found
	CbcParam found = parameters[iParam];
	CbcParameterType type = found.type();
	int valid;
	numberGoodCommands++;
	if (type==GENERALQUERY) {
	  std::cout<<"In argument list keywords have leading - "
	    ", -stdin or just - switches to stdin"<<std::endl;
	  std::cout<<"One command per line (and no -)"<<std::endl;
	  std::cout<<"abcd? gives list of possibilities, if only one + explanation"<<std::endl;
	  std::cout<<"abcd?? adds explanation, if only one fuller help(LATER)"<<std::endl;
	  std::cout<<"abcd without value (where expected) gives current value"<<std::endl;
	  std::cout<<"abcd value or abcd = value sets value"<<std::endl;
	  std::cout<<"Commands are:"<<std::endl;
	  for ( iParam=0; iParam<numberParameters; iParam+=4 ) {
	    int i;
	    for (i=iParam;i<min(numberParameters,iParam+4);i++) 
	      std::cout<<parameters[i].matchName()<<"  ";
	    std::cout<<std::endl;
	  }
	} else if (type<81) {
	  // get next field as double
	  double value = getDoubleField(argc,argv,&valid);
	  if (!valid) {
	    parameters[iParam].setDoubleParameter(*model,value);
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].doubleParameter(*model)<<std::endl;
	  }
	} else if (type<101) {
	  // get next field as double for local use
	  double value = getDoubleField(argc,argv,&valid);
	  if (!valid) {
	    if (!parameters[iParam].checkDoubleParameter(value)) {
	      switch(type) {
	      case DJFIX:
		djFix=value;
		preSolve=5;
                defaultSettings=false; // user knows what she is doing
		break;
	      case GAPRATIO:
		gapRatio=value;
		break;
	      case TIGHTENFACTOR:
		tightenFactor=value;
                defaultSettings=false; // user knows what she is doing
		break;
	      default:
		abort();
	      }
	    }
	  } else if (valid==1) {
	    abort();
	  } else {
	    switch(type) {
	    case DJFIX:
	      value = djFix ;
	      break;
	    case GAPRATIO:
	      value = gapRatio ;
	      break;
	    case TIGHTENFACTOR:
	      value = tightenFactor ;
	      break;
	    default:
	      abort();
	    }
	    std::cout << parameters[iParam].name() << " has value " <<
			 value << std::endl ;
	  }
	} else if (type<201) {
	  // get next field as int
	  int value = getIntField(argc,argv,&valid);
	  if (!valid) {
	    parameters[iParam].setIntParameter(*model,value);
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].intParameter(*model)<<std::endl;
	  }
	} else if (type<301) {
	  // one of several strings
	  std::string value = getString(argc,argv);
	  int action = parameters[iParam].parameterOption(value);
	  if (action<0) {
	    if (value!="EOL") {
	      // no match
	      parameters[iParam].printOptions();
	    } else {
	      // print current value
	      std::cout<<parameters[iParam].name()<<" has value "<<
		parameters[iParam].currentOption()<<std::endl;
	    }
	  } else {
	    parameters[iParam].setCurrentOption(action);
	    // for now hard wired
	    switch (type) {
	    case DIRECTION:
	      if (action==0)
		model->solver()->setObjSense(1);
	      else
		model->solver()->setObjSense(-1);
	      break;
	    case ERRORSALLOWED:
	      allowImportErrors = action;
	      break;
	    case ALGORITHM:
	      algorithm  = action;
              defaultSettings=false; // user knows what she is doing
	      break;
	    case KEEPNAMES:
	      keepImportNames = 1-action;
	      break;
	    case SCALING:
	      doScaling = 1-action;
	      break;
	    case GOMORYCUTS:
              defaultSettings=false; // user knows what she is doing
	      gomoryAction = action;
	      break;
	    case PROBINGCUTS:
              defaultSettings=false; // user knows what she is doing
	      probingAction = action;
	      break;
	    case KNAPSACKCUTS:
              defaultSettings=false; // user knows what she is doing
	      knapsackAction = action;
	      break;
	    case ODDHOLECUTS:
              defaultSettings=false; // user knows what she is doing
	      oddholeAction = action;
	      break;
	    case CLIQUECUTS:
              defaultSettings=false; // user knows what she is doing
	      cliqueAction = action;
	      break;
	    case FLOWCUTS:
              defaultSettings=false; // user knows what she is doing
	      flowAction = action;
	      break;
	    case MIXEDCUTS:
              defaultSettings=false; // user knows what she is doing
	      mixedAction = action;
	      break;
	    case TWOMIRCUTS:
              defaultSettings=false; // user knows what she is doing
	      twomirAction = action;
	      break;
	    case ROUNDING:
              defaultSettings=false; // user knows what she is doing
	      useRounding = action;
	      break;
	    case COSTSTRATEGY:
	      if (action!=1) {
		printf("Pseudo costs not implemented yet\n");
	      } else {
		int numberColumns = model->getNumCols();
		int * sort = new int[numberColumns];
		double * dsort = new double[numberColumns];
		int * priority = new int [numberColumns];
		const double * objective = model->getObjCoefficients();
		int iColumn;
		int n=0;
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  if (model->isInteger(iColumn)) {
		    sort[n]=n;
		    dsort[n++]=-objective[iColumn];
		  }
		}
		CoinSort_2(dsort,dsort+n,sort);
		int level=0;
		double last = -1.0e100;
		for (int i=0;i<n;i++) {
		  int iPut=sort[i];
		  if (dsort[i]!=last) {
		    level++;
		    last=dsort[i];
		  }
		  priority[iPut]=level;
		}
		model->passInPriorities( priority,false);
		delete [] priority;
		delete [] sort;
		delete [] dsort;
	      }
	      break;
	    case PRESOLVE:
                defaultSettings=false; // user knows what she is doing
	      preSolve = action*5;
	      break;
	    case PREPROCESS:
	      preProcess = action;
	      break;
	    case SOLVER:
	    { for (int i = 0 ; i < (int) value.length() ; i++)
		value[i] = tolower(value[i]) ;
	      OsiSolverInterface *newSolver = solvers[value]->clone() ;
	      model->assignSolver(newSolver) ;
	      std::cout << "Solver set to " << value << "." << std::endl ;
	      break ; }
	    default:
	    { std::cerr << "Unrecognized action. Aborting." << std::endl ;
	      abort(); }
	    }
	  }
	} else {
	  // action
	  if (type==EXIT)
	    break; // stop all
	  switch (type) {
	  case SOLVECONTINUOUS:
	    if (goodModel) {
	      model->initialSolve();
	      time2 = CoinCpuTime();
	      totalTime += time2-time1;
	      std::cout<<"Result "<<model->solver()->getObjValue()<<
		" iterations "<<model->solver()->getIterationCount()<<
		" took "<<time2-time1<<" seconds - total "<<totalTime<<std::endl;
	      time1=time2;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
/*
  Run branch-and-cut. First set a few options -- node comparison, scaling. If
  the solver is Clp, consider running some presolve code (not yet converted
  this to generic OSI) with branch-and-cut. If presolve is disabled, or the
  solver is not Clp, simply run branch-and-cut. Print elapsed time at the end.
*/
	  case BAB: // branchAndBound
	  { if (goodModel)
	    { CbcCompareUser compare; // Definition of node choice
            // If user made settings then use them
            if (!defaultSettings) {
	      model->setNodeComparison(compare);
	      OsiSolverInterface * solver = model->solver();
	      if (!doScaling)
		solver->setHintParam(OsiDoScale,false,OsiHintTry);
#ifdef COIN_USE_CLP
	      OsiClpSolverInterface * si =
		dynamic_cast<OsiClpSolverInterface *>(solver) ;
	      if (preSolve&&si != NULL) {
		// get clp itself
		ClpSimplex * modelC = si->getModelPtr();
		if (si->messageHandler()->logLevel())
		  si->messageHandler()->setLogLevel(1);
		if (modelC->tightenPrimalBounds()!=0) {
		  std::cout<<"Problem is infeasible!"<<std::endl;
		  break;
		}
		model->initialSolve();
		// bounds based on continuous
		if (tightenFactor) {
		  if (modelC->tightenPrimalBounds(tightenFactor)!=0) {
		    std::cout<<"Problem is infeasible!"<<std::endl;
		    break;
		  }
		}
		if (gapRatio<1.0e100) {
		  double value = si->getObjValue();
		  double value2 = gapRatio*(1.0e-5+fabs(value));
		  model->setAllowableGap(value2);
		  model->setAllowableFractionGap(gapRatio);
		  std::cout<<"Continuous "<<value
			   <<", so allowable gap set to "<<value2<<std::endl;
		}
		if (djFix<1.0e20) {
		  // do some fixing
		  int numberColumns = modelC->numberColumns();
		  int i;
		  const char * type = modelC->integerInformation();
		  double * lower = modelC->columnLower();
		  double * upper = modelC->columnUpper();
		  double * solution = modelC->primalColumnSolution();
		  double * dj = modelC->dualColumnSolution();
                  int numberFixed=0;
		  for (i=0;i<numberColumns;i++) {
		    if (type[i]) {
		      double value = solution[i];
		      if (value<lower[i]+1.0e-5&&dj[i]>djFix) {
			solution[i]=lower[i];
			upper[i]=lower[i];
                        numberFixed++;
		      } else if (value>upper[i]-1.0e-5&&dj[i]<-djFix) {
			solution[i]=upper[i];
			lower[i]=upper[i];
                        numberFixed++;
		      }
		    }
		  }
                  printf("%d columns fixed\n",numberFixed);
		}
		{
		  // integer presolve
		  CbcModel * model2 = model->integerPresolve();
		  if (model2) {
		    // Do complete search
		    
		    CbcRounding heuristic1(*model2);
		    if (useRounding)
		      model2->addHeuristic(&heuristic1);
                    if (algorithm) {
                      // user wants primal
                      model2->solver()->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
                      model2->solver()->setHintParam(OsiDoDualInResolve,false,OsiHintTry);
                    }
		    model2->branchAndBound();
		    // get back solution
		    model->originalModel(model2,false);
		  } else {
		    // infeasible
		    exit(1);
		  }
		}
	      } else
#endif
	      { if (model->solver()->messageHandler()->logLevel())
		  model->solver()->messageHandler()->setLogLevel(1) ;
		model->initialSolve() ;
		if (gapRatio < 1.0e100)
		{ double value = model->solver()->getObjValue() ;
		  double value2 = gapRatio*(1.0e-5+fabs(value)) ;
		  model->setAllowableGap(value2) ;
		  std::cout << "Continuous " << value
			    << ", so allowable gap set to "
			    << value2 << std::endl ; }
		CbcRounding heuristic1(*model) ;
		if (useRounding)
		  model->addHeuristic(&heuristic1) ;
		// add cut generators if wanted
		if (probingAction==1)
		  model->addCutGenerator(&probingGen,-1,"Probing");
		else if (probingAction==2)
		  model->addCutGenerator(&probingGen,-99,"Probing");
		if (gomoryAction==1)
		  model->addCutGenerator(&gomoryGen,-1,"Gomory");
		else if (gomoryAction==2)
		  model->addCutGenerator(&gomoryGen,-99,"Gomory");
		if (knapsackAction==1)
		  model->addCutGenerator(&knapsackGen,-1,"Knapsack");
		else if (knapsackAction==2)
		  model->addCutGenerator(&knapsackGen,-99,"Knapsack");
		if (oddholeAction==1)
		  model->addCutGenerator(&oddholeGen,-1,"OddHole");
		else if (oddholeAction==2)
		  model->addCutGenerator(&oddholeGen,-99,"OddHole");
		if (cliqueAction==1)
		  model->addCutGenerator(&cliqueGen,-1,"Clique");
		else if (cliqueAction==2)
		  model->addCutGenerator(&cliqueGen,-99,"Clique");
		if (mixedAction==1)
		  model->addCutGenerator(&mixedGen,-1,"MixedintegerRounding");
		else if (mixedAction==2)
		  model->addCutGenerator(&mixedGen,-99,"MixedintegerRounding");
		if (flowAction==1)
		  model->addCutGenerator(&flowGen,-1,"FlowCover");
		else if (flowAction==2)
		  model->addCutGenerator(&flowGen,-99,"FlowCover");
		if (twomirAction==1)
		  model->addCutGenerator(&twomirGen,-1,"TwoMirCuts");
		else if (twomirAction==2)
		  model->addCutGenerator(&twomirGen,-99,"TwoMirCuts");
		model->branchAndBound() ; }
	      if (model->bestSolution())
	      { std::cout << "Optimal solution "
			  << model->solver()->getObjValue() << std::endl ; }
	      else
	      { std::cout << "No integer solution found." << std::endl ; }
			
	      time2 = CoinCpuTime() ;
	      totalTime += time2-time1 ;
	      std::cout << "Result " << model->solver()->getObjValue()
			<< " took " << time2-time1 << " seconds - total "
			<< totalTime << std::endl ;
	      time1 = time2 ;
            } else {
              // User is going to get what I think best
	      if (!doScaling)
		model->solver()->setHintParam(OsiDoScale,false,OsiHintTry);
              model->initialSolve();
              // See if we want preprocessing
              OsiSolverInterface * saveSolver=NULL;
              CglPreProcess process;
              if (preProcess) {
                saveSolver=model->solver()->clone();
                /* Do not try and produce equality cliques and
                   do up to 10 passes */
                OsiSolverInterface * solver2 = process.preProcess(*saveSolver,false,10);
                if (!solver2) {
                  printf("Pre-processing says infeasible\n");
                  break;
                } else {
                  printf("processed model has %d rows and %d columns\n",
                         solver2->getNumRows(),solver2->getNumCols());
                }
                //solver2->resolve();
                // we have to keep solver2 so pass clone
                solver2 = solver2->clone();
                model->assignSolver(solver2);
                model->initialSolve();
              }
	      model->setNodeComparison(compare);
              CbcRounding heuristic1(*model);
              if (useRounding)
                model->addHeuristic(&heuristic1) ;
              // add cut generators if wanted
              if (probingAction==1)
                model->addCutGenerator(&probingGen,-1,"Probing");
              else if (probingAction==2)
                model->addCutGenerator(&probingGen,-99,"Probing");
              if (gomoryAction==1)
                model->addCutGenerator(&gomoryGen,-1,"Gomory");
              else if (gomoryAction==2)
                model->addCutGenerator(&gomoryGen,-99,"Gomory");
              if (knapsackAction==1)
                model->addCutGenerator(&knapsackGen,-1,"Knapsack");
              else if (knapsackAction==2)
                model->addCutGenerator(&knapsackGen,-99,"Knapsack");
              if (oddholeAction==1)
                model->addCutGenerator(&oddholeGen,-1,"OddHole");
              else if (oddholeAction==2)
                model->addCutGenerator(&oddholeGen,-99,"OddHole");
              if (cliqueAction==1)
                model->addCutGenerator(&cliqueGen,-1,"Clique");
              else if (cliqueAction==2)
                model->addCutGenerator(&cliqueGen,-99,"Clique");
              if (mixedAction==1)
                model->addCutGenerator(&mixedGen,-1,"MixedintegerRounding");
              else if (mixedAction==2)
                model->addCutGenerator(&mixedGen,-99,"MixedintegerRounding");
              if (flowAction==1)
                model->addCutGenerator(&flowGen,-1,"FlowCover");
              else if (flowAction==2)
                model->addCutGenerator(&flowGen,-99,"FlowCover");
              if (twomirAction==1)
                model->addCutGenerator(&twomirGen,-1,"TwoMirCuts");
              else if (twomirAction==2)
                model->addCutGenerator(&twomirGen,-99,"TwoMirCuts");
              // Say we want timings
              int numberGenerators = model->numberCutGenerators();
              int iGenerator;
              for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
                CbcCutGenerator * generator = model->cutGenerator(iGenerator);
                generator->setTiming(true);
              }
              // Could tune more
              model->setMinimumDrop(min(1.0,
                                        fabs(model->getMinimizationObjValue())*1.0e-3+1.0e-4));
              
              if (model->getNumCols()<500)
                model->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
              else if (model->getNumCols()<5000)
                model->setMaximumCutPassesAtRoot(100); // use minimum drop
              else
                model->setMaximumCutPassesAtRoot(20);
              model->setMaximumCutPasses(2);
              
              // Do more strong branching if small
              if (model->getNumCols()<5000)
                model->setNumberStrong(20);
              // Switch off strong branching if wanted
              //if (model->getNumCols()>10*model->getNumRows())
              //model->setNumberStrong(0);
              if (model->getNumCols()>2000||model->getNumRows()>1500||
                  model->messageHandler()->logLevel()>1)
                model->setPrintFrequency(100);
              
              model->solver()->setIntParam(OsiMaxNumIterationHotStart,100);
#ifdef COIN_USE_CLP
              OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (model->solver());
              // go faster stripes
              if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
                osiclp->setupForRepeatedUse(2,0);
              }
#endif
              if (gapRatio < 1.0e100)
		{ double value = model->solver()->getObjValue() ;
                double value2 = gapRatio*(1.0e-5+fabs(value)) ;
                model->setAllowableGap(value2) ;
                std::cout << "Continuous " << value
                          << ", so allowable gap set to "
                          << value2 << std::endl ; }
              model->branchAndBound();
              time2 = CoinCpuTime();
              totalTime += time2-time1;
              if (model->getMinimizationObjValue()<1.0e50) {
                // post process
                if (preProcess) {
                  process.postProcess(*model->solver());
                  // Solution now back in saveSolver
                  model->assignSolver(saveSolver);
                }
              }
              std::cout<<"Result "<<model->getObjValue()<<
                " iterations "<<model->getIterationCount()<<
                " nodes "<<model->getNodeCount()<<
                " took "<<time2-time1<<" seconds - total "<<totalTime<<std::endl;
              time1 = time2;
            }
            } else
	    { std::cout << "** Current model not valid" << std::endl ; }
	    break ; }
	  case IMPORT:
	    {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      bool canOpen=false;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
		if (field[0]=='/'||field[0]=='~')
		  fileName = field;
		else
		  fileName = directory+field;
		FILE *fp=fopen(fileName.c_str(),"r");
		if (fp) {
		  // can open - lets go for it
		  fclose(fp);
		  canOpen=true;
		} else {
		  std::cout<<"Unable to open file "<<fileName<<std::endl;
		}
	      }
	      if (canOpen) {
		model->gutsOfDestructor();
		int status =model->solver()->readMps(fileName.c_str(),"");
		if (!status||(status>0&&allowImportErrors)) {
		  // I don't think there is any need for this but ..
		  //OsiWarmStartBasis allSlack;
		  goodModel=true;
		  //model->setBasis(allSlack);
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		} else {
		  // errors
		  std::cout<<"There were "<<status<<
		    " errors on input"<<std::endl;
		}
	      }
	    }
	    break;
	  case EXPORT:
	    {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='~')
		fileName = field;
	      else
		fileName = directory+field;
	      FILE *fp=fopen(fileName.c_str(),"w");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		model->solver()->writeMps(fileName.c_str(),"");
		time2 = CoinCpuTime();
		totalTime += time2-time1;
		time1=time2;
	      }
	    }
	    break;
	  case MAXIMIZE:
	    model->solver()->setObjSense(-1);
	    break;
	  case MINIMIZE:
	    model->solver()->setObjSense(1);
	    break;
	  case DIRECTORY:
	  { directory = getString(argc,argv);
	    if (directory[directory.length()-1] != '/')
	      directory += '/' ;
	    break ; }
	  case STDIN:
	    read_mode=-1;
	    break;
	  case VERSION:
	    std::cout<<"Coin LP version "<<CBCVERSION
		     <<", build "<<__DATE__<<std::endl;
	    break;
	  case UNITTEST:
	    {
	      // okay so there is not a real unit test

	      int status =model->solver()->readMps("../Mps/Sample/p0033.mps",
						   "");
	      assert(!status);
	      model->branchAndBound();
	      model->solver()->resolve();
	      std::cout<<"Optimal solution "<<model->solver()->getObjValue()<<std::endl;
	      assert(fabs(model->solver()->getObjValue()-3089.0)<1.0e-5);
	      fprintf(stderr,"Test was okay\n");
	      status =model->solver()->readMps("../Mps/Sample/p0033.mps",
						   "");
	      assert(!status);
	      model->setCutoff(1.0e20);
	      model->setBestObjectiveValue(1.0e20);
	      model->setMaximumSolutions(1);
	      model->setSolutionCount(0);
	      // Switch off strong branching to give better chance of NOT finding best
	      model->setNumberStrong(0);
	      // Definition of node choice
	      CbcCompareDefault compare(100.0);
	      model->setNodeComparison(compare);
	      model->solver()->resolve();
	      model->branchAndBound();
	      model->solver()->resolve();
	      std::cout<<"partial solution "<<model->solver()->getObjValue()<<std::endl;
	      if (model->solver()->getObjValue()<3090.0) {
		std::cout<<"Got optimal solution by mistake!"<<std::endl;
	      }
	    }
	    break;
	  case MIPLIB:
	    {
	      int mainTest (int argc, const char *argv[]);
	      // create fields for test
	      const char * fields[3];
	      int nFields=1;
	      fields[0]="fake main for miplib";
	      if (directory!="./") {
		fields[1]=("-miplibDir="+directory).c_str();
		nFields=2;
	      }
	      mainTest(nFields,fields);
	    }
	    break;
	  case SOLUTION:
	    if (goodModel) {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      FILE *fp=NULL;
	      if (field=="-"||field=="EOL") {
		// stdout
		fp=stdout;
	      } else {
		if (field[0]=='/'||field[0]=='~')
		  fileName = field;
		else
		  fileName = directory+field;
		fp=fopen(fileName.c_str(),"w");
	      }
	      if (fp) {
		// make fancy later on
		int iRow;
		int numberRows=model->solver()->getNumRows();
		const double * dualRowSolution = model->getRowPrice();
		const double * primalRowSolution =  model->getRowActivity();
		for (iRow=0;iRow<numberRows;iRow++) {
		  fprintf(fp,"%7d ",iRow);
		  fprintf(fp,"%15.8g        %15.8g\n",primalRowSolution[iRow],
			  dualRowSolution[iRow]);
		}
		int iColumn;
		int numberColumns=model->solver()->getNumCols();
		const double * dualColumnSolution = 
		  model->getReducedCost();
		const double * primalColumnSolution = 
		  model->getColSolution();
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  fprintf(fp,"%7d ",iColumn);
		  fprintf(fp,"%15.8g        %15.8g\n",
			  primalColumnSolution[iColumn],
			  dualColumnSolution[iColumn]);
		}
		if (fp!=stdout)
		  fclose(fp);
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	      
	    }
	    break;
	  default:
	    abort();
	  }
	} 
      } else if (!numberMatches) {
	std::cout<<"No match for "<<field<<" - ? for list of commands"
		 <<std::endl;
      } else if (numberMatches==1) {
	if (!numberQuery) {
	  std::cout<<"Short match for "<<field<<" possible completion:"
		   <<std::endl;
	  for ( iParam=0; iParam<numberParameters; iParam++ ) {
	    int match = parameters[iParam].matches(field);
	    if (match) 
	      std::cout<<parameters[iParam].matchName()<<std::endl;
	  }
	} else if (numberQuery) {
	  std::cout<<"Short match for "<<field<<" completion:"
		   <<std::endl;
	  for ( iParam=0; iParam<numberParameters; iParam++ ) {
	    int match = parameters[iParam].matches(field);
	    if (match) {
	      std::cout<<parameters[iParam].matchName()<<" : ";
	      std::cout<<parameters[iParam].shortHelp()<<std::endl;
	    }
	  }
	}
      } else {
	if (!numberQuery) 
	  std::cout<<"Multiple matches for "<<field<<" - possible completions:"
		   <<std::endl;
	else
	  std::cout<<"Completions of "<<field<<":"<<std::endl;
	for ( iParam=0; iParam<numberParameters; iParam++ ) {
	  int match = parameters[iParam].matches(field);
	  if (match) {
	    std::cout<<parameters[iParam].matchName();
	    if (numberQuery>=2) 
	      std::cout<<" : "<<parameters[iParam].shortHelp();
	    std::cout<<std::endl;
	  }
	}
      }
    }
/*
  Final cleanup. Delete the model and the vector of available solvers.
*/
    delete model;

    for (solverMapIter_t solverIter = solvers.begin() ;
	 solverIter != solvers.end() ;
	 solverIter++)
    { if (solverIter->second) delete solverIter->second ; }
  }
  return 0;
}    
