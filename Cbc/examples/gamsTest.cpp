//Gams tests - provided by Stefan Vigerske
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <iostream>
using namespace std;
#include "CoinHelperFunctions.hpp"
#include "CoinError.hpp"
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp" //for CbcSOS
#include "CbcBranchLotsize.hpp" //for CbcLotsize
#include "OsiClpSolverInterface.hpp"
#define testtol 1e-6
/** model sos1a from the GAMS test library
 * http://www.gams.com/testlib/libhtml/sos1a.htm */
void sos1a(int& error_count, int& warning_count);
/** model semicon1 from the GAMS test library
 * http://www.gams.com/testlib/libhtml/semicon1.htm */
void semicon1(int& error_count, int& warning_count);
/** model semiint1 from the GAMS test library
 * http://www.gams.com/testlib/libhtml/semiint1.htm */
void semiint1(int& error_count, int& warning_count);
int main (int argc, const char *argv[])
{
   // only in CoinUtils/trunk: WindowsErrorPopupBlocker();
   int error_count = 0;
   int warning_count = 0;
   
   sos1a(error_count, warning_count);
   cout << "\n***********************\n" << endl;
   semicon1(error_count, warning_count);
   cout << "\n***********************\n" << endl;
   semiint1(error_count, warning_count);
   
   cout << endl << "Finished - there have been " << error_count << " errors and " << warning_count << " warnings." << endl;
   return error_count;
}
void sos1a(int& error_count, int& warning_count) {
   OsiClpSolverInterface solver1;
   
   int numcols = 3;
   int numrows = 1;
   int nnz = 3;
   CoinBigIndex *start = new int[numcols+1];
   int* index = new int[nnz];
   double* value = new double[nnz];
   double *collb = new double[numcols];
   double *colub = new double[numcols];
   double *obj = new double[numcols];
   double *rowlb = new double[numrows];
   double *rowub = new double[numrows];
   // objective
   obj[0] = .9;  obj[1] = 1.;  obj[2] = 1.1;
   
   // column bounds
   collb[0] = 0.;  colub[0] = .8;
   collb[1] = 0.;  colub[1] = .6;
   collb[2] = 0.;  colub[2] = .6;
   // matrix
   start[0] = 0;  index[0] = 0;  value[0] = 1.;
   start[1] = 1;  index[1] = 0;  value[1] = 1.;
   start[2] = 2;  index[2] = 0;  value[2] = 1.;
   start[3] = 3;
   
   // row bounds
   rowlb[0] = -solver1.getInfinity(); rowub[0] = 1.;
   solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
   solver1.setObjSense(-1);
   
  CbcModel model(solver1);
  CbcMain0(model);
  
  int which[3] = { 0, 1, 2 };
   CbcObject* sosobject = new CbcSOS(&model, 3, which, NULL, 0, 1); 
  model.addObjects(1, &sosobject);
  delete sosobject;
  
  const char * argv2[]={"gamstest_sos1a","-solve","-quit"};
  CbcMain1(3,argv2,model);
  cout << endl;
  if (!model.isProvenOptimal()) {
     cerr << "Error: Model sos1a not solved to optimality." << endl;
     ++error_count;
     return; // other tests make no sense ---- memory leak here
  }
  
  OsiSolverInterface* solver = model.solver();
  assert(solver);
  
  cout << "Objective value model: " << model.getObjValue()
    << "\t solver: " << solver->getObjValue()
    << "\t expected: 0.72" << endl;
  if (CoinAbs(model.getObjValue()-0.72)>testtol || CoinAbs(solver->getObjValue()-0.72)>testtol) {
     cerr << "Error: Objective value incorrect." << endl;
     ++error_count;
  }
  
  cout << "Primal value variable 0 in model: " << model.bestSolution()[0]
    << "\t in solver: " << solver->getColSolution()[0]
    << "\t expected: 0.8" << endl;
  if (CoinAbs(model.bestSolution()[0]-0.8)>testtol || CoinAbs(solver->getColSolution()[0]-0.8)>testtol) {
     cerr << "Error: Primal value incorrect." << endl;
     ++error_count;
  }
  cout << "Primal value variable 1 in model: " << model.bestSolution()[1]
    << "\t in solver: " << solver->getColSolution()[1]
    << "\t expected: 0.0" << endl;
  if (CoinAbs(model.bestSolution()[1])>testtol || CoinAbs(solver->getColSolution()[1])>testtol) {
     cerr << "Error: Primal value incorrect." << endl;
     ++error_count;
  }
  cout << "Primal value variable 2 in model: " << model.bestSolution()[2]
    << "\t in solver: " << solver->getColSolution()[2]
    << "\t expected: 0.0" << endl;
  if (CoinAbs(model.bestSolution()[2])>testtol || CoinAbs(solver->getColSolution()[2])>testtol) {
     cerr << "Error: Primal value incorrect." << endl;
     ++error_count;
  }
  delete[] start;
   delete[] index;
   delete[] value;
   delete[] collb;
   delete[] colub;
   delete[] obj;
   delete[] rowlb;
   delete[] rowub;
}
void semicon1(int& error_count, int& warning_count) {
   OsiClpSolverInterface solver1;
   
   int numcols = 4; // s, pup, plo, x
   int numrows = 3; // bigx, smallx, f
   int nnz = 6;
   CoinBigIndex *start = new int[numcols+1];
   int* index = new int[nnz];
   double* value = new double[nnz];
   double *collb = new double[numcols];
   double *colub = new double[numcols];
   double *obj = new double[numcols];
   double *rowlb = new double[numrows];
   double *rowub = new double[numrows];
   // objective
   obj[0] = 0;  obj[1] = 1.;  obj[2] = 1;   obj[3] = 0;
   
   // column bounds
   collb[0] = 0.;  colub[0] = 10.;
   collb[1] = 0.;  colub[1] = solver1.getInfinity();
   collb[2] = 0.;  colub[2] = solver1.getInfinity();
   collb[3] = 0.;  colub[3] = solver1.getInfinity();
   // matrix
   start[0] = 0;  index[0] = 2;  value[0] =  1.;
   start[1] = 1;  index[1] = 0;  value[1] = -1.;
   start[2] = 2;  index[2] = 1;  value[2] =  1.;
   start[3] = 3;  index[3] = 0;  value[3] =  1.;
                  index[4] = 1;  value[4] =  1.;
                  index[5] = 2;  value[5] =  1.;
   start[4] = nnz;
   
   // row bounds
   rowlb[0] = -solver1.getInfinity(); rowub[0] = 8.9;
   rowlb[1] = 8.9;  rowub[1] = solver1.getInfinity();
   rowlb[2] = 10.;  rowub[2] = 10.;
   solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
   for (int testcase = 0; testcase < 5; ++testcase) {
      CbcModel model(solver1);
      CbcMain0(model);
      
      double points[4] = { 0., 0., 0., 10. };
      double objval;
      double primalval[4];
      double redcost[4];
      double row2marg;
      redcost[1] = 1.0;
      redcost[2] = 1.0;
      redcost[3] = 0.0;
      switch(testcase) {
         case 0:
            points[2] = 0.;
            objval = 0.;
            primalval[0] = 1.1;
            primalval[1] = 0.0;
            primalval[2] = 0.0;
            primalval[3] = 8.9;
            redcost[0] = 0.0;
            row2marg = 0.0; 
            break;
         case 1:
            points[2] = 1.;
            objval = 0.;
            primalval[0] = 1.1;
            primalval[1] = 0.0;
            primalval[2] = 0.0;
            primalval[3] = 8.9;
            redcost[0] = 0.0;
            row2marg = 0.0; 
            break;
         case 2:
            points[2] = 1.5;
            objval = 0.4;
            primalval[0] = 1.5;
            primalval[1] = 0.0;
            primalval[2] = 0.4;
            primalval[3] = 8.5;
            redcost[0] = 1.0;
            row2marg = -1.0; 
            break;
         case 3:
            points[2] = 2.1;
            objval = 1.0;
            primalval[0] = 2.1;
            primalval[1] = 0.0;
            primalval[2] = 1.0;
            primalval[3] = 7.9;
            redcost[0] = 1.0;
            row2marg = -1.0; 
            break;
         case 4:
            points[2] = 2.8;
            objval = 1.1;
            primalval[0] = 0.0;
            primalval[1] = 1.1;
            primalval[2] = 0.0;
            primalval[3] = 10.0;
            redcost[0] = -1.0;
            row2marg = 1.0; 
            break;
      }
      
      CbcObject* semiconobject = new CbcLotsize(&model, 0, 2, points, true); 
      model.addObjects(1, &semiconobject);
      delete semiconobject;
      
      cout << "\nSolving semicon1 model for lotsize variable being either 0 or between " << points[2] << " and 10.\n" << endl; 
      const char * argv2[]={"gamstest_semicon1","-solve","-quit"};
      CbcMain1(3,argv2,model);
      cout << endl;
      if (!model.isProvenOptimal()) {
         cerr << "Error: Model semicon1 not solved to optimality." << endl;
         ++error_count;
         continue; // other tests make no sense
      }
      OsiSolverInterface* solver = model.solver();
      assert(solver);
      cout << "Objective value in model: " << model.getObjValue()
      << "\t in solver: " << solver->getObjValue()
      << "\t expected: " << objval << endl;
      if (CoinAbs(model.getObjValue()-objval)>testtol || CoinAbs(solver->getObjValue()-objval)>testtol) {
         cerr << "Error: Objective value incorrect." << endl;
         ++error_count;
      }
      for (int i=0; i<numcols; ++i) {
         cout << "Primal value variable " << i << " in model: " << model.bestSolution()[i]
        << "\t in solver: " << solver->getColSolution()[i]
        << "\t expected: " << primalval[i]
        << endl;
      if (CoinAbs(model.bestSolution()[i]-primalval[i])>testtol || CoinAbs(solver->getColSolution()[i]-primalval[i])>testtol) {
         cerr << "Error: Primal value incorrect." << endl;
         ++error_count;
      }
      }
      for (int i=0; i<numcols; ++i) {
         cout << "Reduced cost variable " << i << " in model: " << model.getReducedCost()[i]
        << "\t in solver: " << solver->getReducedCost()[i]
        << "\t expected: " << redcost[i]
        << endl;
      if (CoinAbs(model.getReducedCost()[i]-redcost[i])>testtol || CoinAbs(solver->getReducedCost()[i]-redcost[i])>testtol) {
         cerr << "Warning: Reduced cost incorrect." << endl;
         ++warning_count;
      }
      }
      
      cout << "Row 2 marginal (price) in model: " << model.getRowPrice()[2] 
      << "\t in solver: " << solver->getRowPrice()[2]
        << "\t expected: " << row2marg << endl;
    if (CoinAbs(model.getRowPrice()[2]-row2marg)>testtol || CoinAbs(solver->getRowPrice()[2]-row2marg)>testtol) {
       cerr << "Warning: Row price incorrect." << endl;
       ++warning_count;
    }
   
   }
   
   delete[] start;
   delete[] index;
   delete[] value;
   delete[] collb;
   delete[] colub;
   delete[] obj;
   delete[] rowlb;
   delete[] rowub;
}
void semiint1(int& error_count, int& warning_count) {
   OsiClpSolverInterface solver1;
   
   int numcols = 4; // s, pup, plo, x
   int numrows = 3; // bigx, smallx, f
   int nnz = 6;
   CoinBigIndex *start = new int[numcols+1];
   int* index = new int[nnz];
   double* value = new double[nnz];
   double *collb = new double[numcols];
   double *colub = new double[numcols];
   double *obj = new double[numcols];
   double *rowlb = new double[numrows];
   double *rowub = new double[numrows];
   // objective
   obj[0] = 0;  obj[1] = 1.;  obj[2] = 1;   obj[3] = 0;
   
   // column bounds
   collb[0] = 0.;  colub[0] = 10.;
   collb[1] = 0.;  colub[1] = solver1.getInfinity();
   collb[2] = 0.;  colub[2] = solver1.getInfinity();
   collb[3] = 0.;  colub[3] = solver1.getInfinity();
   // matrix
   start[0] = 0;  index[0] = 2;  value[0] =  1.;
   start[1] = 1;  index[1] = 0;  value[1] = -1.;
   start[2] = 2;  index[2] = 1;  value[2] =  1.;
   start[3] = 3;  index[3] = 0;  value[3] =  1.;
                  index[4] = 1;  value[4] =  1.;
                  index[5] = 2;  value[5] =  1.;
   start[4] = nnz;
   
   // row bounds
   rowlb[0] = -solver1.getInfinity(); rowub[0] = 7.9;
   rowlb[1] = 7.9;  rowub[1] = solver1.getInfinity();
   rowlb[2] = 10.;  rowub[2] = 10.;
   solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
   solver1.setInteger(0); 
   
   for (int testcase = 0; testcase < 6; ++testcase) {
      CbcModel model(solver1);
      CbcMain0(model);
      
      double points[4] = { 0., 0., 0., 10. };
      double objval;
      double primalval[4];
      double redcost[4];
      double row2marg;
      redcost[2] = 1.0;
      redcost[3] = 0.0;
      switch(testcase) {
         case 0:
            points[2] = 0.;
            objval = 0.;
            primalval[0] = 2.0;
            primalval[1] = 0.0;
            primalval[2] = 0.0;
            primalval[3] = 8.0;
            redcost[0] = 0.0;
            redcost[1] = 1.0;
            row2marg = 0.0; 
            break;
         case 1:
            points[2] = 1.;
            objval = 0.1;
            primalval[0] = 2.0;
            primalval[1] = 0.1;
            primalval[2] = 0.0;
            primalval[3] = 8.0;
            redcost[0] = -1.0;
            redcost[1] =  0.0;
            row2marg = 1.0; 
            break;
         case 2:
            points[2] = 2.;
            objval = 0.1;
            primalval[0] = 2.0;
            primalval[1] = 0.1;
            primalval[2] = 0.0;
            primalval[3] = 8.0;
            redcost[0] = -1.0;
            redcost[1] =  0.0;
            row2marg = 1.0; 
            break;
         case 3:
            points[2] = 3.;
            objval = 0.9;
            primalval[0] = 3.0;
            primalval[1] = 0.0;
            primalval[2] = 0.9;
            primalval[3] = 7.0;
            redcost[0] = 1.0;
            redcost[1] = 1.0;
            row2marg = -1.0; 
            break;
         case 4:
            points[2] = 4.;
            objval = 1.9;
            primalval[0] = 4.0;
            primalval[1] = 0.0;
            primalval[2] = 1.9;
            primalval[3] = 6.0;
            redcost[0] = 1.0;
            redcost[1] = 1.0;
            row2marg = -1.0; 
            break;
         case 5:
            points[2] = 5.;
            objval = 2.1;
            primalval[0] = 0.0;
            primalval[1] = 2.1;
            primalval[2] = 0.0;
            primalval[3] = 10.0;
            redcost[0] = -1.0;
            redcost[1] =  0.0;
            row2marg = 1.0; 
            break;
      }
      if (testcase!=3) {
	CbcObject* semiintobject = new CbcLotsize(&model, 0, 2, points, true); 
	model.addObjects(1, &semiintobject);
	delete semiintobject;
      } else {
	// Do alternate way
	double points2[]={0.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	int numberPoints = (int) sizeof(points2)/sizeof(double);
	CbcObject* semiintobject = new CbcLotsize(&model, 0, numberPoints,
						  points2); 
	model.addObjects(1, &semiintobject);
	delete semiintobject;
      }
      
      cout << "\nSolving semiint1 model for integer lotsize variable being either 0 or between " << points[2] << " and 10.\n" << endl; 
      const char * argv2[]={"gamstest_semiint1","-solve","-quit"};
      CbcMain1(3,argv2,model);
      cout << endl;
      if (!model.isProvenOptimal()) {
         cerr << "Error: Model semiint1 not solved to optimality." << endl;
         ++error_count;
         continue; // other tests make no sense
      }
      OsiSolverInterface* solver = model.solver();
      assert(solver);
      cout << "Objective value in model: " << model.getObjValue()
      << "\t in solver: " << solver->getObjValue()
      << "\t expected: " << objval << endl;
      if (CoinAbs(model.getObjValue()-objval)>testtol || CoinAbs(solver->getObjValue()-objval)>testtol) {
         cerr << "Error: Objective value incorrect." << endl;
         ++error_count;
      }
      for (int i=0; i<numcols; ++i) {
         cout << "Primal value variable " << i << " in model: " << model.bestSolution()[i]
        << "\t in solver: " << solver->getColSolution()[i]
        << "\t expected: " << primalval[i]
        << endl;
      if (CoinAbs(model.bestSolution()[i]-primalval[i])>testtol || CoinAbs(solver->getColSolution()[i]-primalval[i])>testtol) {
         cerr << "Error: Primal value incorrect." << endl;
         ++error_count;
      }
      }
      for (int i=0; i<numcols; ++i) {
         cout << "Reduced cost variable " << i << " in model: " << model.getReducedCost()[i]
        << "\t in solver: " << solver->getReducedCost()[i]
        << "\t expected: " << redcost[i]
        << endl;
      if (CoinAbs(model.getReducedCost()[i]-redcost[i])>testtol || CoinAbs(solver->getReducedCost()[i]-redcost[i])>testtol) {
         cerr << "Warning: Reduced cost incorrect." << endl;
         ++warning_count;
      }
      }
      
      cout << "Row 2 marginal (price) in model: " << model.getRowPrice()[2] 
      << "\t in solver: " << solver->getRowPrice()[2]
        << "\t expected: " << row2marg << endl;
    if (CoinAbs(model.getRowPrice()[2]-row2marg)>testtol || CoinAbs(solver->getRowPrice()[2]-row2marg)>testtol) {
       cerr << "Warning: Row price incorrect." << endl;
       ++warning_count;
    }
   
   }
   
   delete[] start;
   delete[] index;
   delete[] value;
   delete[] collb;
   delete[] colub;
   delete[] obj;
   delete[] rowlb;
   delete[] rowub;
}
