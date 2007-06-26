// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <math.h>
#include <cfloat>

#include "CoinPragma.hpp"
//#include "CoinHelperFunctions.hpp"
//#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"

#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"

#include "CoinMessageHandler.hpp"
#include "OsiClpSolverInterface.hpp"

//  bobe including extras.h to get strdup()
#if defined(__MWERKS__) 
// #include <extras.h>  // bobe 06-02-14
#endif

// Get C stuff but with extern C
#define CBC_EXTERN_C
#include "Coin_C_defines.h"

const int  VERBOSE = 0;

// To allow call backs
class Cbc_MessageHandler 
 : public CoinMessageHandler 
  {
  
public:
  /**@name Overrides */
  //@{
  virtual int print();
  //@}
  /**@name set and get */
  //@{
  /// Model
  const Cbc_Model * model() const;
  void setModel(Cbc_Model * model);
  /// Call back
  void setCallBack(cbc_callback callback);
  //@}

  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  Cbc_MessageHandler();
  /// Constructor with pointer to model
  Cbc_MessageHandler(Cbc_Model * model,
         FILE * userPointer=NULL);
  /** Destructor */
  virtual ~Cbc_MessageHandler();
  //@}

  /**@name Copy method */
  //@{
  /** The copy constructor. */
  Cbc_MessageHandler(const Cbc_MessageHandler&);
  /** The copy constructor from an CoinSimplexMessageHandler. */
  Cbc_MessageHandler(const CoinMessageHandler&);
  
  Cbc_MessageHandler& operator=(const Cbc_MessageHandler&);
  /// Clone
  virtual CoinMessageHandler * clone() const ;
  //@}
   
    
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer back to model
  Cbc_Model * model_;
  /// call back
  cbc_callback callback_;
  //@}
};


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
Cbc_MessageHandler::Cbc_MessageHandler () 
  : CoinMessageHandler(),
    model_(NULL),
    callback_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
Cbc_MessageHandler::Cbc_MessageHandler (const Cbc_MessageHandler & rhs) 
 : CoinMessageHandler(rhs),
    model_(rhs.model_),
    callback_(rhs.callback_)
{  
}

Cbc_MessageHandler::Cbc_MessageHandler (const CoinMessageHandler & rhs) 
  : CoinMessageHandler(),
    model_(NULL),
    callback_(NULL)
{  
}

// Constructor with pointer to model
Cbc_MessageHandler::Cbc_MessageHandler(Cbc_Model * model,
               FILE * userPointer)
  : CoinMessageHandler(),
    model_(model),
    callback_(NULL)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
Cbc_MessageHandler::~Cbc_MessageHandler ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
Cbc_MessageHandler &
Cbc_MessageHandler::operator=(const Cbc_MessageHandler& rhs)
{
  if (this != &rhs) {
    CoinMessageHandler::operator=(rhs);
    model_ = rhs.model_;
    callback_ = rhs.callback_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CoinMessageHandler * Cbc_MessageHandler::clone() const
{
  return new Cbc_MessageHandler(*this);
}
int 
Cbc_MessageHandler::print()
{
  if (callback_) {
    int messageNumber = currentMessage().externalNumber();
    if (currentSource()!="Cbc")
      messageNumber += 1000000;
    int i;
    int nDouble=numberDoubleFields();
    assert (nDouble<=10);
    double vDouble[10];
    for (i=0;i<nDouble;i++)
      vDouble[i]=doubleValue(i);
    int nInt=numberIntFields();
    assert (nInt<=10);
    int vInt[10];
    for (i=0;i<nInt;i++)
      vInt[i]=intValue(i);
    int nString=numberStringFields();
    assert (nString<=10);
    char * vString[10];
    for (i=0;i<nString;i++) {
      std::string value = stringValue(i);
      vString[i]=strdup(value.c_str());
    }
    callback_(model_,messageNumber,
        nDouble,vDouble,
        nInt,vInt,
        nString,vString);
    for (i=0;i<nString;i++) 
      free(vString[i]);
    
  }
  return CoinMessageHandler::print();
  return 0;
}
const Cbc_Model *
Cbc_MessageHandler::model() const
{
  return model_;
}
void 
Cbc_MessageHandler::setModel(Cbc_Model * model)
{
  model_ = model;
}
// Call back
void 
Cbc_MessageHandler::setCallBack(cbc_callback callback)
{
  callback_ = callback;
}
/**
  *
  *  C Interface Routines
  *
  */
#include "Cbc_C_Interface.h"
#include <string>
#include <stdio.h>
#include <iostream>

#if defined(__MWERKS__) 
#pragma export on
#endif

/* Version */
COINLIBAPI double COINLINKAGE Cbc_getVersion()
{
  double v = 1.0;
  return v;
}

/* Default Cbc_Model constructor */
COINLIBAPI Cbc_Model *  COINLINKAGE 
Cbc_newModel()
{
  const char prefix[] = "Cbc_C_Interface::Cbc_newModel(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  Cbc_Model * model = new Cbc_Model;
  OsiClpSolverInterface solver1;
  model->solver_    = &solver1;
  model->solver_->OsiClpSolverInterface::setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model->model_     = new CbcModel(solver1);
  model->handler_   = NULL;
  model->information_ = NULL;
  
  if (VERBOSE>0) printf("%s return\n",prefix);
  return model;
}
/* Cbc_Model Destructor */
COINLIBAPI void COINLINKAGE 
Cbc_deleteModel(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_deleteModel(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix); fflush(stdout);
  
  if (VERBOSE>1) printf("%s delete model->model_\n",prefix); fflush(stdout);
  delete model->model_;
  
  if (VERBOSE>1) printf("%s delete model->handler_\n",prefix); fflush(stdout);
  delete model->handler_;
  
  if (VERBOSE>1) printf("%s free model->information_\n",prefix); fflush(stdout);
  if (model->information_) free(model->information_);
  
  if (VERBOSE>1) printf("%s delete model\n",prefix); fflush(stdout);
  delete model;

  if (VERBOSE>0) printf("%s return\n",prefix); fflush(stdout);
}

/* Loads a problem (the constraints on the
    rows are given by lower and upper bounds). If a pointer is NULL then the
    following values are the default:
    <ul>
    <li> <code>colub</code>: all columns have upper bound infinity
    <li> <code>collb</code>: all columns have lower bound 0 
    <li> <code>rowub</code>: all rows have upper bound infinity
    <li> <code>rowlb</code>: all rows have lower bound -infinity
    <li> <code>obj</code>: all variables have 0 objective coefficient
    </ul>

   Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps). 
*/
COINLIBAPI void COINLINKAGE 
Cbc_loadProblem (Cbc_Model * model,  const int numcols, const int numrows,
     const CoinBigIndex * start, const int* index,
     const double* value,
     const double* collb, const double* colub,   
     const double* obj,
     const double* rowlb, const double* rowub)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_loadProblem(): ";
//  const int  VERBOSE = 2;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  OsiSolverInterface * solver = model->model_->solver();

  if (VERBOSE>1) {
    printf("%s numcols = %i, numrows = %i\n", 
      prefix, numcols, numrows);
    printf("%s model = %p, start = %p, index = %p, value = %p\n",
      prefix, static_cast<void*>(model), static_cast<const void*>(start),
	   static_cast<const void*>(index), static_cast<const void*>(value));
    printf("%s collb = %p, colub = %p, obj = %p, rowlb = %p, rowub = %p\n",
      prefix, static_cast<const void*>(collb),
	   static_cast<const void*>(colub), static_cast<const void*>(obj),
	   static_cast<const void*>(rowlb), static_cast<const void*>(rowub));
  }

  if (VERBOSE>1) printf("%s Calling solver->loadProblem()\n",prefix);
  fflush(stdout);
  
  if (1) {
    solver->loadProblem(numcols,numrows,start,index,value,
             collb,colub,obj,rowlb,rowub);
  } else {
    solver->loadProblem(0,0,NULL,NULL,NULL,
             NULL,NULL,NULL,NULL,NULL);
  }
  if (VERBOSE > 1) printf("%s Finished solver->loadProblem()\n",prefix);
  fflush(stdout);

  if (VERBOSE>0) printf("%s return\n",prefix);
} //  Cbc_loadProblem()

/* Read an mps file from the given filename */
COINLIBAPI int COINLINKAGE 
Cbc_readMps(Cbc_Model * model,const char *filename)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_readMps(): ";
//  const int  VERBOSE = 2;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  if (VERBOSE>1) printf("%s filename = '%s'\n",prefix,filename);

  int result = 1;
  result = model->model_->solver()->readMps(filename);
  assert(result==0);

  if (VERBOSE>0) printf("%s return %i\n",prefix, result);
  return result;
}
/* Write an mps file from the given filename */
COINLIBAPI void COINLINKAGE 
Cbc_writeMps(Cbc_Model * model,const char *filename)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_writeMps(): ";
//  const int  VERBOSE = 2;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  if (VERBOSE>1) printf("%s filename = '%s'\n",prefix,filename);

  model->model_->solver()->writeMps(filename,"mps",Cbc_optimizationDirection(model));

  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
}
/* Integer information */
COINLIBAPI char * COINLINKAGE 
Cbc_integerInformation(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_integerInformation(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int col;
  int numcols = Cbc_getNumCols(model);

  // allocate model->information_ if null
  // this is freed in Cbc_deleteModel() if not null
  if (!model->information_) 
    model->information_ = (char *) malloc(numcols*sizeof(char));

  for (col=0; col<numcols; col++) 
    if (model->model_->solver()->isContinuous(col)) 
      model->information_[col] = 0;
    else
      model->information_[col] = 1;

  char * result = model->information_;

  if (VERBOSE>0) printf("%s return %p\n",prefix,result);
  return result;
}
/* Copy in integer information */
COINLIBAPI void COINLINKAGE 
Cbc_copyInIntegerInformation(Cbc_Model * model, const char * information)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_copyInIntegerInformation(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int col;
  int numcols = Cbc_getNumCols(model);
  for (col=0; col<numcols; col++) 
    if (information[col]) 
      model->model_->solver()->setInteger(col);
    else
      model->model_->solver()->setContinuous(col);
  
  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Drop integer informations */
COINLIBAPI void COINLINKAGE 
Cbc_deleteIntegerInformation(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_deleteIntegerInformation(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

//  available through 
//    OsiClpSolverInterface::setContinuous
//tbd  model->model_->deleteIntegerInformation();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Resizes rim part of model  */
COINLIBAPI void COINLINKAGE 
Cbc_resize (Cbc_Model * model, int newNumberRows, int newNumberColumns)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_resize(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find this in Cbc, Osi, or OsiClp
//tbd  model->model_->solver()->resize(newNumberRows,newNumberColumns);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Deletes rows */
COINLIBAPI void COINLINKAGE 
Cbc_deleteRows(Cbc_Model * model, int number, const int * which)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_deleteRows(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  OsiSolverInterface * solver = model->model_->solver();
  solver->deleteRows(number,which);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Add rows */
COINLIBAPI void COINLINKAGE 
Cbc_addRows(Cbc_Model * model, const int number, const double * rowLower, 
         const double * rowUpper,
         const int * rowStarts, const int * columns,
         const double * elements)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_addRows(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// available through OsiClp
//tbd  model->model_->addRows(number,rowLower,rowUpper,rowStarts,columns,elements);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}

/* Deletes columns */
COINLIBAPI void COINLINKAGE 
Cbc_deleteColumns(Cbc_Model * model, int number, const int * which)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_deleteColumns(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  OsiSolverInterface * solver = model->model_->solver();
  solver->deleteCols(number,which);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Add columns */
COINLIBAPI void COINLINKAGE 
Cbc_addColumns(Cbc_Model * model, int number, const double * columnLower, 
      const double * columnUpper,
      const double * objective,
      const int * columnStarts, const int * rows,
      const double * elements)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_addColumns(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// available through OsiClp
//tbd  model->model_->addColumns(number,columnLower,columnUpper,objective,
//tbd			    columnStarts,rows,elements);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Drops names - makes lengthnames 0 and names empty */
COINLIBAPI void COINLINKAGE 
Cbc_dropNames(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dropNames(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
// cannot find names in Cbc, Osi, or OsiClp
//tbd  model->model_->dropNames();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Copies in names */
COINLIBAPI void COINLINKAGE 
Cbc_copyNames(Cbc_Model * model, const char * const * rowNamesIn,
        const char * const * columnNamesIn)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_copyNames(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find names in Cbc, Osi, or OsiClp
/*clean
  int iRow;
  std::vector<std::string> rowNames;
  int numberRows = model->model_->getNumRows();
  rowNames.reserve(numberRows);
  for (iRow=0;iRow<numberRows;iRow++) {
    rowNames.push_back(rowNamesIn[iRow]);
  }
  
  int iColumn;
  std::vector<std::string> columnNames;
  int numberColumns = model->model_->getNumCols();
  columnNames.reserve(numberColumns);
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    columnNames.push_back(columnNamesIn[iColumn]);
  }
  model->model_->copyNames(rowNames,columnNames);
*/

  if (VERBOSE>0) printf("%s return\n",prefix);
}

/* Number of rows */
COINLIBAPI int COINLINKAGE 
Cbc_numberRows(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_numberRows(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->getNumRows();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Number of columns */
COINLIBAPI int COINLINKAGE 
Cbc_numberColumns(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_numberColumns(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->getNumCols();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Primal tolerance to use */
COINLIBAPI double COINLINKAGE 
Cbc_primalTolerance(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_primalTolerance(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  model->model_->solver()->getDblParam(OsiPrimalTolerance,result) ;

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setPrimalTolerance(Cbc_Model * model,  double value) 
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setPrimalTolerance(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->solver()->setDblParam(OsiPrimalTolerance,value) ;

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Dual tolerance to use */
COINLIBAPI double COINLINKAGE 
Cbc_dualTolerance(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dualTolerance(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  model->model_->solver()->getDblParam(OsiDualTolerance,result) ;

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setDualTolerance(Cbc_Model * model,  double value) 
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setDualTolerance(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->solver()->setDblParam(OsiDualTolerance,value) ;

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Integer tolerance to use */
COINLIBAPI double COINLINKAGE 
Cbc_integerTolerance(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_primalTolerance(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  result = model->model_->getDblParam(CbcModel::CbcIntegerTolerance) ;

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setIntegerTolerance(Cbc_Model * model,  double value) 
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setPrimalTolerance(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->setDblParam(CbcModel::CbcIntegerTolerance,value) ;

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Dual objective limit */
COINLIBAPI double COINLINKAGE 
Cbc_dualObjectiveLimit(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dualObjectiveLimit(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  model->model_->solver()->getDblParam(OsiDualObjectiveLimit,result) ;

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setDualObjectiveLimit(Cbc_Model * model, double value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setDualObjectiveLimit(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->solver()->setDblParam(OsiDualObjectiveLimit,value) ;

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Objective offset */
COINLIBAPI double COINLINKAGE 
Cbc_objectiveOffset(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_objectiveOffset(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  return model->model_->objectiveOffset();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setObjectiveOffset(Cbc_Model * model, double value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setObjectiveOffset(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find names in Cbc, Osi, or OsiClp
//tbd  model->model_->solver()->setObjectiveOffset(value);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Fills in array with problem name  */
COINLIBAPI void COINLINKAGE 
Cbc_problemName(Cbc_Model * model, int maxNumberCharacters, char * array)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_problemName(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  std::string name;
  model->model_->solver()->getStrParam(OsiProbName,name);
  maxNumberCharacters = CoinMin(maxNumberCharacters,(int)strlen(name.c_str()));
  strncpy(array,name.c_str(),maxNumberCharacters-1);
  array[maxNumberCharacters-1]='\0';

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Sets problem name.  Must have \0 at end.  */
COINLIBAPI int COINLINKAGE 
Cbc_setProblemName(Cbc_Model * model, int maxNumberCharacters, char * array)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setProblemName(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  bool result = false;
  result = model->model_->solver()->setStrParam(OsiProbName, array);

  if (VERBOSE>0) printf("%s return\n",prefix);
  return (result) ? 1 : 0;
}
/* Number of iterations */
COINLIBAPI int COINLINKAGE 
Cbc_numberIterations(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_numberIterations(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->getIterationCount();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setNumberIterations(Cbc_Model * model, int numberIterations)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setNumberIterations(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find this in Cbc, Osi, or OsiClp
//tbd  model->model_->setNumberIterations(numberIterations);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
}
/* Maximum number of iterations */
COINLIBAPI int 
Cbc_maximumIterations(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_maximumIterations(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
// cannot find this in Cbc, Osi, or OsiClp
//tbd  result = model->model_->solver()->maximumIterations();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setMaximumIterations(Cbc_Model * model, int value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setMaximumIterations(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find this in Cbc, Osi, or OsiClp
//tbd  model->model_->setMaximumIterations(value);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Maximum number of nodes */
COINLIBAPI int 
Cbc_maxNumNode(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_maxNumNode(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->getIntParam(CbcModel::CbcMaxNumNode);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setMaxNumNode(Cbc_Model * model, int value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setMaxNumNode(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->setIntParam(CbcModel::CbcMaxNumNode, value);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Maximum number of solutions */
COINLIBAPI int 
Cbc_maxNumSol(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::maxNumSol(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->getIntParam(CbcModel::CbcMaxNumSol);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setMaxNumSol(Cbc_Model * model, int value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setMaxNumSol(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->setIntParam(CbcModel::CbcMaxNumSol, value);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Maximum time in seconds (from when set called) */
COINLIBAPI double COINLINKAGE 
Cbc_maximumSeconds(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_maximumSeconds(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  result = model->model_->getDblParam(CbcModel::CbcMaximumSeconds);

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setMaximumSeconds(Cbc_Model * model, double value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setMaximumSeconds(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->setDblParam(CbcModel::CbcMaximumSeconds, value);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Returns true if hit maximum iteratio`ns (or time) */
COINLIBAPI int COINLINKAGE 
Cbc_hitMaximumIterations(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_hitMaximumIterations(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->solver()->hitMaximumIterations() ? 1 : 0;
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Status of problem:
   0 - optimal
   1 - primal infeasible
   2 - dual infeasible
   3 - stopped on iterations etc
   4 - stopped due to errors
*/
COINLIBAPI int COINLINKAGE 
Cbc_status(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_status(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->status();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Set problem status */
COINLIBAPI void COINLINKAGE 
Cbc_setProblemStatus(Cbc_Model * model, int problemStatus)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setProblemStatus(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find this in Cbc, Osi, or OsiClp
//tbd  model->model_->setProblemStatus(problemStatus);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Secondary status of problem - may get extended
   0 - none
   1 - primal infeasible because dual limit reached
   2 - scaled problem optimal - unscaled has primal infeasibilities
   3 - scaled problem optimal - unscaled has dual infeasibilities
   4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
*/
COINLIBAPI int COINLINKAGE 
Cbc_secondaryStatus(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_secondaryStatus(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
// cannot find this in Cbc, Osi, or OsiClp
  result = model->model_->secondaryStatus();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setSecondaryStatus(Cbc_Model * model, int status)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setSecondaryStatus(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// cannot find this in Cbc, Osi, or OsiClp
//tbd  model->model_->setSecondaryStatus(status);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
COINLIBAPI double COINLINKAGE 
Cbc_optimizationDirection(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_optimizationDirection(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  result = model->model_->getObjSense();

  if (VERBOSE>0) printf("%s return %g\n",prefix, result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setOptimizationDirection(Cbc_Model * model, double value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setOptimizationDirection(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin, value = %g\n",prefix,value);

  model->model_->setObjSense(value);
//  model->model_->solver()->setObjSense(value);
  
  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Primal row solution */
COINLIBAPI double * COINLINKAGE 
Cbc_primalRowSolution(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_primalRowSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
  
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->primalRowSolution();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return result;
}
/* Primal column solution */
COINLIBAPI double * COINLINKAGE 
Cbc_primalColumnSolution(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_primalColumnSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//  result = model->model_->getColSolution();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return result;
}
/* Dual row solution */
COINLIBAPI double * COINLINKAGE 
Cbc_dualRowSolution(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dualRowSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->dualRowSolution();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Reduced costs */
COINLIBAPI double * COINLINKAGE 
Cbc_dualColumnSolution(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dualColumnSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->dualColumnSolution();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Row lower */
COINLIBAPI double* COINLINKAGE 
Cbc_rowLower(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_rowLower(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->rowLower();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Row upper  */
COINLIBAPI double* COINLINKAGE 
Cbc_rowUpper(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_rowUpper(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->rowUpper();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Objective Coefficients */
COINLIBAPI double * COINLINKAGE 
Cbc_objective(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_objective(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//  result = model->model_->objective();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Column Lower */
COINLIBAPI double * COINLINKAGE 
Cbc_columnLower(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_columnLower(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->columnLower();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Column Upper */
COINLIBAPI double * COINLINKAGE 
Cbc_columnUpper(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_columnUpper(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
// cannot find this in Cbc, Osi, or OsiClp
// may have to make it somehow
//tbd  return model->model_->columnUpper();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return NULL;
}
/* Number of elements in matrix */
COINLIBAPI int COINLINKAGE 
Cbc_getNumElements(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getNumElements(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
  result = model->model_->getNumElements();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}

// Column starts in matrix 
COINLIBAPI const CoinBigIndex * COINLINKAGE Cbc_getVectorStarts(Cbc_Model * model)
{
  const CoinPackedMatrix * matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  return (matrix == NULL) ? NULL : matrix->getVectorStarts();
}
// Row indices in matrix 
COINLIBAPI const int * COINLINKAGE 
Cbc_getIndices(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getIndices(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  const int * result = NULL;
  const CoinPackedMatrix * matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getIndices();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
} 

// Column vector lengths in matrix 
COINLIBAPI const int * COINLINKAGE 
Cbc_getVectorLengths(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getVectorLengths(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  const int * result = NULL;
  const CoinPackedMatrix * matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getVectorLengths();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
} 
  
// Element values in matrix 
COINLIBAPI const double * COINLINKAGE 
Cbc_getElements(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getElements(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  const double * result = NULL;
  const CoinPackedMatrix * matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getElements();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
// ======================================================================

/* Objective value */
COINLIBAPI double COINLINKAGE 
Cbc_objectiveValue(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_objectiveValue(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double result = 0.0;
  result = model->model_->getObjValue();

  if (VERBOSE>0) printf("%s return %g\n",prefix, result);
  return result;
}
/* Infeasibility/unbounded ray (NULL returned if none/wrong)
   Up to user to use delete [] on these arrays.  */
COINLIBAPI double * COINLINKAGE 
Cbc_infeasibilityRay(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_infeasibilityRay(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
  // lots of rays (probably too many) are available in
  // OsiClpSolverInterface::getDualRays()
  //
//tbd  result = model->model_->infeasibilityRay();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return result;
}
COINLIBAPI double * COINLINKAGE 
Cbc_unboundedRay(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_unboundedRay(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  double * result = NULL;
  // lots of rays (probably too many) are available in
  // OsiClpSolverInterface::getPrimalRays()
  //
//tbd  result = model->model_->unboundedRay();
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix, static_cast<void*>(result));
  return result;
}
/* See if status array exists (partly for OsiClp) */
COINLIBAPI int COINLINKAGE 
Cbc_statusExists(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_statusExists(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  int result = 0;
//tbd  result = model->model_->statusExists() ? 1 : 0;
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Return address of status array (char[numberRows+numberColumns]) */
COINLIBAPI void  COINLINKAGE 
Cbc_getBasisStatus(Cbc_Model * model, int * cstat, int * rstat)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getBasisStatus(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

// have to figure this out
//tbd  model->model_->solver()->getBasisStatus(cstat, rstat);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
}
/* Copy in status vector */
COINLIBAPI void COINLINKAGE 
setBasisStatus(Cbc_Model * model,  int * cstat, int * rstat)
{
  const char prefix[] = "Cbc_C_Interface::setBasisStatus(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

//  model->model_->solver()->setBasisStatus(cstat, rstat);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}

/* User pointer for whatever reason */
COINLIBAPI void COINLINKAGE 
Cbc_setUserPointer (Cbc_Model * model, void * pointer)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setUserPointer(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  // not sure what this is for
  //
//tbd  model->model_->setUserPointer(pointer);
  if (VERBOSE>0) printf("%s WARNING: NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
COINLIBAPI void * COINLINKAGE 
Cbc_getUserPointer (Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getUserPointer(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  void * result = NULL;
  // not sure what this is for
  //
//tbd result = model->model_->getUserPointer();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %p\n",prefix,result);
  return result;
}
/* Pass in Callback function */
COINLIBAPI void COINLINKAGE 
Cbc_registerCallBack(Cbc_Model * model, 
         cbc_callback userCallBack)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_registerCallBack(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  // Will be copy of users one
  delete model->handler_;
  model->handler_ = new Cbc_MessageHandler(*(model->model_->messageHandler()));
  model->handler_->setCallBack(userCallBack);
  model->handler_->setModel(model);
  model->model_->passInMessageHandler(model->handler_);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Unset Callback function */
COINLIBAPI void COINLINKAGE 
Cbc_clearCallBack(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_clearCallBack(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  delete model->handler_;
  model->handler_=NULL;

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Amount of print out:
   0 - none
   1 - just final
   2 - just factorizations
   3 - as 2 plus a bit more
   4 - verbose
   above that 8,16,32 etc just for selective debug
*/
COINLIBAPI void COINLINKAGE 
Cbc_setLogLevel(Cbc_Model * model, int value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setLogLevel(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  if (VERBOSE>1) printf("%s value = %i\n",prefix, value);
  
  model->model_->messageHandler()->setLogLevel(value);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
COINLIBAPI int COINLINKAGE 
Cbc_logLevel(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_logLevel(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->messageHandler()->logLevel();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* length of names (0 means no names0 */
COINLIBAPI int COINLINKAGE 
Cbc_lengthNames(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_lengthNames(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->lengthNames();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Fill in array (at least lengthNames+1 long) with a row name */
COINLIBAPI void COINLINKAGE 
Cbc_rowName(Cbc_Model * model, int iRow, char * name)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_rowName(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  sprintf(name,"ROW%5i",iRow);
// cannot find names in Cbc, Osi, or OsiClp
//tbd  std::string rowName=model->model_->rowName(iRow);
//tbd  strcpy(name,rowName.c_str());

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Fill in array (at least lengthNames+1 long) with a column name */
// cannot find names in Cbc, Osi, or OsiClp
COINLIBAPI void COINLINKAGE 
Cbc_columnName(Cbc_Model * model, int iColumn, char * name)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_columnName(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  sprintf(name,"COL%5i",iColumn);
//tbd  std::string columnName= model->model_->columnName(iColumn);
//tbd  strcpy(name,columnName.c_str());

  if (VERBOSE>0) printf("%s return\n",prefix);
}

/* General branch and bound solve algorithm which can do presolve.
   See  CbcSolve.hpp for options
*/
COINLIBAPI int COINLINKAGE 
Cbc_initialSolve(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_initialSolve(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  model->model_->initialSolve();
  result = model->model_->status();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* General solve algorithm which can do presolve.
   See  CbcModel.hpp for options
*/
COINLIBAPI int COINLINKAGE 
Cbc_branchAndBound(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_branchAndBound(): ";
//  const int  VERBOSE = 3;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  if (VERBOSE>2) Cbc_printModel(model,prefix);
  try {
    model->model_->branchAndBound();
    model->model_->solver()->resolve();
  } catch (CoinError e) {
    printf("%s ERROR: %s::%s, %s\n", prefix,
      e.className().c_str(), e.methodName().c_str(), e.message().c_str());
  }
  result = model->model_->status();

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) */
COINLIBAPI void COINLINKAGE 
Cbc_scaling(Cbc_Model * model, int mode)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_scaling(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  OsiSolverInterface * solver = model->model_->solver();
  solver->setHintParam(OsiDoScale,mode);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Gets scalingFlag */
COINLIBAPI int COINLINKAGE 
Cbc_scalingFlag(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_scalingFlag(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// try to use OsiSolverInterface::getHintParam(OsiDoScale, ???)
//tbd  result = model->model_->scalingFlag();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Crash - at present just aimed at dual, returns
   -2 if dual preferred and crash basis created
   -1 if dual preferred and all slack basis preferred
   0 if basis going in was not all slack
   1 if primal preferred and all slack basis preferred
   2 if primal preferred and crash basis created.
   
   if gap between bounds <="gap" variables can be flipped
   
   If "pivot" is
   0 No pivoting (so will just be choice of algorithm)
   1 Simple pivoting e.g. gub
   2 Mini iterations
*/
COINLIBAPI int COINLINKAGE 
Cbc_crash(Cbc_Model * model, double gap,int pivot)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_crash(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->crash(gap,pivot);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* If problem is primal feasible */
COINLIBAPI int COINLINKAGE 
Cbc_primalFeasible(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_primalFeasible(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isProvenPrimalInfeasible() ? 0 : 1;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* If problem is dual feasible */
COINLIBAPI int COINLINKAGE 
Cbc_dualFeasible(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dualFeasible(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isProvenDualInfeasible() ? 0 : 1;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Dual bound */
COINLIBAPI double COINLINKAGE 
Cbc_dualBound(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_dualBound(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
// cannot find in Cbc, Osi, or OsiClp
//tbd  result = model->model_->dualBound();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setDualBound(Cbc_Model * model, double value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setDualBound(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
// cannot find names in Cbc, Osi, or OsiClp
//tbd  model->model_->setDualBound(value);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Infeasibility cost */
COINLIBAPI double COINLINKAGE 
Cbc_infeasibilityCost(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_infeasibilityCost(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->solver()->infeasibilityCost();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setInfeasibilityCost(Cbc_Model * model, double value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setInfeasibilityCost(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
// cannot find names in Cbc, Osi, or OsiClp
//tbd  model->model_->setInfeasibilityCost(value);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Perturbation:
   50  - switch on perturbation
   100 - auto perturb if takes too long (1.0e-6 largest nonzero)
   101 - we are perturbed
   102 - don't try perturbing again
   default is 100
   others are for playing
*/
COINLIBAPI int COINLINKAGE 
Cbc_perturbation(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_perturbation(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->perturbation();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setPerturbation(Cbc_Model * model, int value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setPerturbation(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
// cannot find names in Cbc, Osi, or OsiClp
//tbd  model->model_->setPerturbation(value);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Current (or last) algorithm */
COINLIBAPI int COINLINKAGE 
Cbc_algorithm(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setPerturbation(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->algorithm();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Set algorithm */
COINLIBAPI void COINLINKAGE 
Cbc_setAlgorithm(Cbc_Model * model, int value)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setAlgorithm(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
// cannot find names in Cbc, Osi, or OsiClp
//tbd  model->model_->setAlgorithm(value);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return\n",prefix);
}
/* Sum of dual infeasibilities */
COINLIBAPI double COINLINKAGE 
Cbc_sumDualInfeasibilities(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_sumDualInfeasibilities(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->sumDualInfeasibilities();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
/* Number of dual infeasibilities */
COINLIBAPI int COINLINKAGE 
Cbc_numberDualInfeasibilities(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_numberDualInfeasibilities(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd  result = model->model_->numberDualInfeasibilities();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Sum of primal infeasibilities */
COINLIBAPI double COINLINKAGE 
Cbc_sumPrimalInfeasibilities(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_sumPrimalInfeasibilities(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
// cannot find names in Cbc, Osi, or OsiClp
//tbd result = model->model_->sumPrimalInfeasibilities();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
/* Number of primal infeasibilities */
COINLIBAPI int COINLINKAGE 
Cbc_numberPrimalInfeasibilities(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_numberPrimalInfeasibilities(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
//tbd  result = model->model_->getContinuousInfeasibilities();
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Save model to file, returns 0 if success.  This is designed for
   use outside algorithms so does not save iterating arrays etc.
   It does not save any messaging information. 
   Does not save scaling values.
   It does not know about all types of virtual functions.
*/
COINLIBAPI int COINLINKAGE 
Cbc_saveModel(Cbc_Model * model, const char * fileName)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_saveModel(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// there is a writeMPS method in Osi 
//tbd  result = model->model_->saveModel(fileName);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Restore model from file, returns 0 if success,
   deletes current model */
COINLIBAPI int COINLINKAGE 
Cbc_restoreModel(Cbc_Model * model, const char * fileName)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_restoreModel(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
// there is a readMPS method in Osi
//tbd  result = model->model_->restoreModel(fileName);
  if (VERBOSE>0) printf("%s WARNING:  NOT IMPLEMENTED\n",prefix);

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
  
/** Call this to really test if a valid solution can be feasible
    Solution is number columns in size.
    If fixVariables true then bounds of continuous solver updated.
    Returns objective value (worse than cutoff if not feasible)
*/
COINLIBAPI void COINLINKAGE 
Cbc_checkSolution(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_checkSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  // see CbcModel::checkSolution(double cutoff, const double * solution,
  //	       bool fixVariables);
//  model->model_->checkSolution();

  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
}
/* Number of rows */
COINLIBAPI int COINLINKAGE 
Cbc_getNumRows(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getNumRows(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->getNumRows();

  if (VERBOSE>0) printf("%s return %d\n",prefix, result);
  return result;
}
/* Number of columns */
COINLIBAPI int COINLINKAGE 
Cbc_getNumCols(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getNumCols(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->getNumCols();

  if (VERBOSE>0) printf("%s return %d\n",prefix, result);
  return result;
}
/* Number of iterations */
COINLIBAPI int COINLINKAGE 
Cbc_getIterationCount(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getIterationCount(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->getIterationCount();

  if (VERBOSE>0) printf("%s return %d\n",prefix, result);
  return result;
}
/* Are there a numerical difficulties? */
COINLIBAPI int COINLINKAGE 
Cbc_isAbandoned(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isAbandoned(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->isAbandoned() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Is optimality proven? */
COINLIBAPI int COINLINKAGE 
Cbc_isProvenOptimal(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isProvenOptimal(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->isProvenOptimal() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Is primal infeasiblity proven? */
COINLIBAPI int COINLINKAGE 
Cbc_isProvenPrimalInfeasible(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isProvenPrimalInfeasible(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isProvenPrimalInfeasible() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Is dual infeasiblity proven? */
COINLIBAPI int COINLINKAGE 
Cbc_isProvenDualInfeasible(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isProvenDualInfeasible(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isProvenDualInfeasible() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Is the given primal objective limit reached? */
COINLIBAPI int COINLINKAGE 
Cbc_isPrimalObjectiveLimitReached(Cbc_Model * model) 
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isPrimalObjectiveLimitReached(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isPrimalObjectiveLimitReached() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Is the given dual objective limit reached? */
COINLIBAPI int COINLINKAGE 
Cbc_isDualObjectiveLimitReached(Cbc_Model * model) 
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isDualObjectiveLimitReached(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isDualObjectiveLimitReached() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Iteration limit reached? */
COINLIBAPI int COINLINKAGE 
Cbc_isIterationLimitReached(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isIterationLimitReached(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  OsiSolverInterface * solver = model->model_->solver();
  result = solver->isIterationLimitReached() ? 1 : 0;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/* Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
COINLIBAPI double COINLINKAGE 
Cbc_getObjSense(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getObjSense(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
  result = model->model_->getObjSense();

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
/* Primal row solution */
COINLIBAPI const double * COINLINKAGE 
Cbc_getRowActivity(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getRowActivity(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getRowActivity();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Primal column solution */
COINLIBAPI const double * COINLINKAGE 
Cbc_getColSolution(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getColSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getColSolution();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
COINLIBAPI void COINLINKAGE 
Cbc_setColSolution(Cbc_Model * model, const double * input)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setColSolution(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  OsiSolverInterface * solver = model->model_->solver();
  solver->setColSolution(input);

  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
}
/* Dual row solution */
COINLIBAPI const double * COINLINKAGE 
Cbc_getRowPrice(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getRowPrice(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getRowPrice();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Reduced costs */
COINLIBAPI const double * COINLINKAGE 
Cbc_getReducedCost(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getReducedCost(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getReducedCost();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Row lower */
COINLIBAPI const double * COINLINKAGE 
Cbc_getRowLower(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getRowLower(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getRowLower();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Row upper  */
COINLIBAPI const double * COINLINKAGE 
Cbc_getRowUpper(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getRowUpper(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getRowUpper();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Objective Coefficients */
COINLIBAPI const double * COINLINKAGE 
Cbc_getObjCoefficients(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getObjCoefficients(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getObjCoefficients();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Column Lower */
COINLIBAPI const double * COINLINKAGE 
Cbc_getColLower(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getColLower(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getColLower();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Column Upper */
COINLIBAPI const double * COINLINKAGE 
Cbc_getColUpper(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getColUpper(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  const double * result = NULL;
  result = model->model_->getColUpper();

  if (VERBOSE>0)
    printf("%s return %p\n",prefix, static_cast<const void*>(result));
  return result;
}
/* Objective value */
COINLIBAPI double COINLINKAGE 
Cbc_getObjValue(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getObjValue(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
  result = model->model_->getObjValue();

  if (VERBOSE>0)
    printf("%s return %g\n",prefix, result);
  return result;
}
/* Print model */
COINLIBAPI void COINLINKAGE 
Cbc_printModel(Cbc_Model * model, const char * argPrefix)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_printModel(): ";
  const int  VERBOSE = 4;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  CbcModel *cbc_model = model->model_;
  int numrows    = cbc_model->getNumRows();
  int numcols    = cbc_model->getNumCols();
  int numelem    = cbc_model->getNumElements();
  const CoinPackedMatrix * matrix = cbc_model->solver()->getMatrixByCol();
  const CoinBigIndex     * start  = matrix->getVectorStarts();
  const int              * index  = matrix->getIndices();
  const double           * value  = matrix->getElements();
  const double           * collb  = cbc_model->getColLower();
  const double           * colub  = cbc_model->getColUpper();
  const double           * obj    = cbc_model->getObjCoefficients();
  const double           * rowlb  = cbc_model->getRowLower();
  const double           * rowub  = cbc_model->getRowUpper();

  printf("%s numcols = %i, numrows = %i, numelem = %i\n", 
    argPrefix, numcols, numrows, numelem);
  printf("%s model = %p, start = %p, index = %p, value = %p\n",
    argPrefix, static_cast<void*>(model), static_cast<const void*>(start),
	 static_cast<const void*>(index), static_cast<const void*>(value));
  matrix->dumpMatrix(NULL);
  {
    int i;
    for (i=0; i<=numcols; i++) 
      printf("%s start[%i] = %i\n", argPrefix, i, start[i]);
    for (i=0; i< numelem; i++)
      printf("%s index[%i] = %i, value[%i] = %g\n",
        argPrefix, i, index[i], i, value[i]);
  }
  
  printf("%s collb = %p, colub = %p, obj = %p, rowlb = %p, rowub = %p\n",
    argPrefix, static_cast<const void*>(collb),
	 static_cast<const void*>(colub), static_cast<const void*>(obj),
	 static_cast<const void*>(rowlb), static_cast<const void*>(rowub));
  printf("%s optimization direction = %g\n",argPrefix, Cbc_optimizationDirection(model));
  printf("  (1 - minimize, -1 - maximize, 0 - ignore)\n");
  {
    int i;
    for (i=0; i<numcols; i++) 
      printf("%s collb[%i] = %g, colub[%i] = %g, obj[%i] = %g\n",
        argPrefix, i, collb[i], i, colub[i], i, obj[i]);
    for (i=0; i< numrows; i++)
      printf("%s rowlb[%i] = %g, rowub[%i] = %g\n",
        argPrefix, i, rowlb[i], i, rowub[i]);
  }

  if (VERBOSE>0) printf("%s return\n",prefix);
}  // Cbc_printModel()

COINLIBAPI int COINLINKAGE 
Cbc_isInteger(Cbc_Model * model, int i)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isInteger(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  bool result = false;
  result = model->model_->isInteger(i);

  if (VERBOSE>0) printf("%s return %i\n",prefix, result);
  return (result) ? 1 : 0;
}

COINLIBAPI double COINLINKAGE 
Cbc_cpuTime(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_cpuTime(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  double result = 0;
  result = CoinCpuTime() ;

  if (VERBOSE>0) printf("%s return %g\n",prefix,result);
  return result;
}
/** Number of nodes explored in B&B tree */
COINLIBAPI int COINLINKAGE Cbc_getNodeCount(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getNodeCount(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  int result = 0;
  result = model->model_->getNodeCount() ;

  if (VERBOSE>0) printf("%s return %i\n",prefix,result);
  return result;
}
/** Return a copy of this model */
COINLIBAPI Cbc_Model * COINLINKAGE Cbc_clone(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_clone(): ";
  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  Cbc_Model * result = new Cbc_Model;
  result->model_     = new CbcModel(*(model->model_));
  result->solver_    = dynamic_cast< OsiClpSolverInterface*> (result->model_->solver());
  result->handler_   = NULL;

  if (VERBOSE>0) printf("%s return\n",prefix);
  return model;
}
/** Set this the variable to be continuous */
COINLIBAPI Cbc_Model * COINLINKAGE Cbc_setContinuous(Cbc_Model * model, int iColumn)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setContinuous(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  
  model->model_->solver()->setContinuous(iColumn);

  if (VERBOSE>0) printf("%s return\n",prefix);
  return model;
}
/* Delete all object information */
COINLIBAPI void  COINLINKAGE 
Cbc_deleteObjects(Cbc_Model * model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_deleteObjects(): ";
//  const int  VERBOSE = 2;
  if (VERBOSE>0) printf("%s begin\n",prefix);

  model->model_->deleteObjects();

  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
}
/* Add an SOS constraint to the model */
COINLIBAPI void  COINLINKAGE 
Cbc_addSOS_Dense(Cbc_Model * model, int numObjects, const int * len,
           const int ** which, const double * weights, const int type)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_addSOS_Dense(): ";
  const int  VERBOSE = 2;
  if (VERBOSE>0) printf("%sbegin\n",prefix);
  
  assert(1>0);// this is probably broken
  int i,j,l;
  const int *w;
  // I think this is a different model due to overriding = operator
  CbcModel m = *(model->model_); 
  
  CbcObject ** objects = new CbcObject * [numObjects];

  if (VERBOSE>1) printf("%s numObjects = %i\n",prefix,numObjects);
  for (i=0;i<numObjects;i++) {
    if (VERBOSE>1) {
      printf("%s len[%i] = %i, identifier = %i, type = %i\n",
        prefix,i,len[i], i, type); 
      fflush(stdout);
      for (j=0;j<len[i];j++) {
        if (VERBOSE>2 || j==0 || j==(len[i]-1)) {
          printf("%s which[%i][%i] = %d, weights[%i] = %g\n",
            prefix, i,j, which[i][j], j, weights[j]); 
          fflush(stdout);
        }
      }
    }
    
    // Make a CbcSOS and assign it to objects
    printf("%s len[%i] = %i\n",prefix, i, len[i]);
    l = len[i];
    w = which[i];
    if (VERBOSE>1) printf("%s new CbcSOS()\n",prefix);
    // ***
    objects[i] = new CbcSOS(model->model_, (int)(len[i]), 
      (const int*)which[i], (const double*)weights, (int)i, (int)type);
    // ***
    if (objects[i] == NULL) {
      printf("%s ERROR: objects[%i] == NULL\n",prefix,i); fflush(stdout);
      assert(objects[i] != NULL);
    }
  }
  if (VERBOSE>1) printf("%s calling addObjects()\n",prefix);
  fflush(stdout);
  model->model_->addObjects(numObjects, objects);
  if (VERBOSE>1) printf("%s finished addObjects()\n",prefix);
  
  for (i=0;i<numObjects;i++) delete objects[i]; 
  delete [] objects;

  if (VERBOSE>0) printf("%sreturn\n",prefix);
  return;
}
/* Add SOS constraints to the model using row-order matrix */
COINLIBAPI void  COINLINKAGE 
Cbc_addSOS_Sparse(Cbc_Model * model, const int * rowStarts,
           const int * rowIndices, const double * weights, const int type)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_addSOS_Sparse(): ";
//  const int  VERBOSE = 1;
  if (VERBOSE>0) printf("%sbegin\n",prefix);
  
  int numRows = Cbc_numberRows(model);
  if (VERBOSE>0) printf("%s numRows = %i\n",prefix,numRows);
  
  // The passed sparse matrix must have the same number of rows as the model
  assert(numRows==Cbc_numberRows(model)); 
  
  int row, i;
  const int *colIndex;
  const double *colWeight;
  
  // loop on rows and count number of objects according to numWeights>0
  int numObjects = 0;
  for (row=0; row<numRows; row++) {
    if (VERBOSE>2) {
      printf("%s row = %i\n",prefix,row);
      printf("%s rowStarts[%i] = %i\n",prefix,row,rowStarts[row]);
      printf("%s rowStarts[%i+1] = %i\n",prefix,row,rowStarts[row+1]);
      fflush(stdout);
    }
    const int numWeights = rowStarts[row+1]-rowStarts[row];
    if (VERBOSE>2) printf("%s  numWeights = %i\n",prefix,numWeights);
    if (numWeights>0) numObjects++;    
  }
  
  // make objects
  CbcObject ** objects = new CbcObject * [numObjects];
//  if (VERBOSE>1) printf("%s numObjects = %i, objects = %X\n",prefix,numObjects,objects);

  // loop on rows and make an object when numWeights>0
  int objNum = 0;
  for (row=0; row<numRows; row++) {
    if (VERBOSE>2) {
      printf("%s row = %i\n",prefix,row);
      printf("%s rowStarts[%i] = %i\n",prefix,row,rowStarts[row]);
      printf("%s rowStarts[%i+1] = %i\n",prefix,row,rowStarts[row+1]);
    }
    const int numWeights = rowStarts[row+1]-rowStarts[row];
    if (VERBOSE>2) printf("%s  numWeights = %i\n",prefix,numWeights);
    colIndex    = rowIndices+rowStarts[row];
    colWeight   = weights+rowStarts[row];
    if (numWeights>0) {
      // Make a CbcSOS and assign it to objects
      if (VERBOSE>3) {
        for (i=0;i<numWeights;i++) {
          printf("%s  colIndex [%i] = %i\n",prefix,i,colIndex[i]);
          printf("%s  colWeight[%i] = %f\n",prefix,i,colWeight[i]);
        }
        fflush(stdout);
      }
      objects[objNum] = new CbcSOS(model->model_, (int)(numWeights), 
        (const int*)colIndex, (const double*)colWeight, (int)objNum, (int)type);
//      if (VERBOSE>2) printf("%s objects[%i] = %X\n",prefix,objNum,objects[objNum]);
      if (objects[objNum] == NULL) {
        printf("%s ERROR: objects[%i] == NULL\n",prefix,objNum); fflush(stdout);
        assert(objects[objNum] != NULL);
      }
      objNum++;
    } 
  }
  if (VERBOSE>2) {
    printf("%s calling addObjects()\n",prefix);
/*
    printf("%s numObjects = %i, objects = %X\n",prefix,numObjects,objects);
    for (row=0; row<numObjects; row++)
      printf("%s  objects[%i] = %X\n",prefix,row,objects[row]);   
*/   
  }
  fflush(stdout);
  model->model_->addObjects(numObjects, objects);
  if (VERBOSE>1) printf("%s finished addObjects()\n",prefix);
  
  for (objNum=0;objNum<numObjects;objNum++) delete objects[objNum]; 
  delete [] objects;

  if (VERBOSE>0) printf("%sreturn\n",prefix);
  return;
}

/* Print the solution */
COINLIBAPI void  COINLINKAGE 
Cbc_printSolution(Cbc_Model * model)
{
  {
    //
    //  Now to print out row solution.  The methods used return const 
    //  pointers - which is of course much more virtuous.
    //  
    //  This version just does non-zero columns
    //
    
    // * Rows 
    
    int numberRows = Cbc_getNumRows(model);
    int iRow;
    
    
    const double * rowPrimal = Cbc_getRowActivity(model);
    // * Alternatively getReducedCost(model) 
    const double * rowDual = Cbc_getRowPrice(model);
    // * Alternatively getColLower(model) 
    const double * rowLower = Cbc_getRowLower(model);
    // * Alternatively getColUpper(model) 
    const double * rowUpper = Cbc_getRowUpper(model);
    printf("--------------------------------------\n");

    // * If we have not kept names (parameter to readMps) this will be 0 
  //    assert(Cbc_lengthNames(model));
    
    printf("                       Primal          Dual         Lower         Upper\n");
    for (iRow=0;iRow<numberRows;iRow++) {
      double value;
      value = rowPrimal[iRow];
      if (value>1.0e-8||value<-1.0e-8) {
      	char name[20];
  //      	Cbc_columnName(model,iColumn,name);
        sprintf(name,"ROW%5i",iRow);
      	printf("%6d %8s",iRow,name);
      	printf(" %13g",rowPrimal[iRow]);
      	printf(" %13g",rowDual[iRow]);
      	printf(" %13g",rowLower[iRow]);
      	printf(" %13g",rowUpper[iRow]);
      	printf("\n");
      }
    }
    printf("--------------------------------------\n");
  }
  {
    //
    //  Now to print out column solution.  The methods used return const 
    //  pointers - which is of course much more virtuous.
    //  
    //  This version just does non-zero columns
    //
    //
    
    // * Columns 
    
    int numberColumns = Cbc_numberColumns(model);
    int iColumn;
    
    
    // * Alternatively getColSolution(model) 
    const double * columnPrimal = Cbc_getColSolution(model);
    // * Alternatively getReducedCost(model) 
    const double * columnDual = Cbc_getReducedCost(model);
    // * Alternatively getColLower(model) 
    const double * columnLower = Cbc_getColLower(model);
    // * Alternatively getColUpper(model) 
    const double * columnUpper = Cbc_getColUpper(model);
    // * Alternatively getObjCoefficients(model) 
    const double * columnObjective = Cbc_getObjCoefficients(model);

    const char * isInteger = Cbc_integerInformation(model);

    printf("--------------------------------------\n");

    // * If we have not kept names (parameter to readMps) this will be 0 
//    assert(Cbc_lengthNames(model));
    
    printf("                       Primal          Dual         Lower         Upper          Cost     isInteger\n");
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value;
      value = columnPrimal[iColumn];
      if (value>1.0e-8||value<-1.0e-8) {
      	char name[20];
//      	Cbc_columnName(model,iColumn,name);
        sprintf(name,"COL%5i",iColumn);
      	printf("%6d %8s",iColumn,name);
      	printf(" %13g",columnPrimal[iColumn]);
      	printf(" %13g",columnDual[iColumn]);
      	printf(" %13g",columnLower[iColumn]);
      	printf(" %13g",columnUpper[iColumn]);
      	printf(" %13g",columnObjective[iColumn]);
      	printf(" %13i",isInteger[iColumn]);
      	printf("\n");
      }
    }
    printf("--------------------------------------\n");
  }
  if (0) Cbc_printModel(model, "cbc::main(): ");
  return;
}
#if defined(__MWERKS__) 
#pragma export off
#endif
