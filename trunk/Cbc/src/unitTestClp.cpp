// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cstdio>
#include <string>
#include <iostream>

#include "CoinTime.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "OsiClpSolverInterface.hpp"

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

//#############################################################################

// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
  std::cout <<msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}

//#############################################################################

void CbcClpUnitTest (const CbcModel & saveModel, std::string& miplibDir,
		     bool unitTestOnly)
{
  unsigned int m ;
  // See if files exist
  FILE * fp;
  bool doTest=false;

  // Set directory containing miplib data files.
  std::string test1 = miplibDir +"p0033";
  fp=fopen(test1.c_str(),"r");
  if (fp) {
    doTest=true;
    fclose(fp);
  }
#ifdef COIN_HAS_ZLIB
  test1 += ".gz";
  fp=fopen(test1.c_str(),"r");
  if (fp) {
    doTest=true;
    fclose(fp);
  }
#endif
  if (!doTest) {
    printf("Not doing miplib run as can't find mps files - ? .gz without libz\n");
    return;
  }
  /*
    Vectors to hold test problem names and characteristics. The objective value
    after optimization (objValue) must agree to the specified tolerance
    (objValueTol).
  */
  std::vector<std::string> mpsName ;
  std::vector<int> nRows ;
  std::vector<int> nCols ;
  std::vector<double> objValueC ;
  std::vector<double> objValue ;
  std::vector<int> strategy ;
  /*
    And a macro to make the vector creation marginally readable.
  */
#define PUSH_MPS(zz_mpsName_zz,\
		 zz_nRows_zz,zz_nCols_zz,zz_objValue_zz,zz_objValueC_zz, \
                 zz_strategy_zz) \
  mpsName.push_back(zz_mpsName_zz) ; \
  nRows.push_back(zz_nRows_zz) ; \
  nCols.push_back(zz_nCols_zz) ; \
  objValueC.push_back(zz_objValueC_zz) ; \
  strategy.push_back(zz_strategy_zz) ; \
  objValue.push_back(zz_objValue_zz) ;

  if (unitTestOnly) {
    PUSH_MPS("p0033",16,33,3089,2520.57,7);
    PUSH_MPS("p0201",133,201,7615,6875.0,7);
    PUSH_MPS("p0548",176,548,8691,315.29,7);
  } else {
    /*
      Load up the problem vector. Note that the row counts here include the
      objective function.
    */
    // 0 for no test, 1 for some, 2 for many, 3 for all
#define HOWMANY 2
#if HOWMANY
#if HOWMANY>1
    PUSH_MPS("10teams",230,2025,924,917,7);
#endif
    PUSH_MPS("air03",124,10757,340160,338864.25,7);
#if HOWMANY==3
    PUSH_MPS("air04",823,8904,56137,55535.436,8);
    PUSH_MPS("air05",426,7195,26374,25877.609,8);
#endif
  //    PUSH_MPS("arki001",1048,1388,7580813.0459,7579599.80787,7);
    PUSH_MPS("bell3a",123,133,878430.32,862578.64,7);
#if HOWMANY>1
    PUSH_MPS("bell5",91,104,8966406.49,8608417.95,7);
#endif
    PUSH_MPS("blend2",274,353,7.598985,6.9156751140,7);
#if HOWMANY>1
    PUSH_MPS("cap6000",2176,6000,-2451377,-2451537.325,7);
#endif
    //    PUSH_MPS("dano3mip",3202,13873,728.1111,576.23162474,7);
    //PUSH_MPS("danoint",664,521,65.67,62.637280418,7);
    PUSH_MPS("dcmulti",290,548,188182,183975.5397,7);
    PUSH_MPS("dsbmip",1182,1886,-305.19817501,-305.19817501,7);
    PUSH_MPS("egout",98,141,568.101,149.589,7);
    PUSH_MPS("enigma",21,100,0.0,0.0,7);
#if HOWMANY==3
    PUSH_MPS("fast0507",507,63009,174,172.14556668,7);
#endif
    PUSH_MPS("fiber",363,1298,405935.18000,156082.51759,7);
#if HOWMANY>1
    PUSH_MPS("fixnet6",478,878,3983,1200.88,7);
#endif
    PUSH_MPS("flugpl",18,18,1201500,1167185.7,7);
    PUSH_MPS("gen",780,870,112313,112130.0,7);
#if HOWMANY>1
    PUSH_MPS("gesa2",1392,1224,25779856.372,25476489.678,7);
    PUSH_MPS("gesa2_o",1248,1224,25779856.372,25476489.678,7);
#endif
    PUSH_MPS("gesa3",1368,1152,27991042.648,27833632.451,7);
    PUSH_MPS("gesa3_o",1224,1152,27991042.648,27833632.451,7);
    PUSH_MPS("gt2",29,188,21166.000,13460.233074,7);
#if HOWMANY==3
    PUSH_MPS("harp2",112,2993,-73899798.00,-74353341.502,7);
#endif
    PUSH_MPS("khb05250",101,1350,106940226,95919464.0,7);
#if HOWMANY>1
    PUSH_MPS("l152lav",97,1989,4722,4656.36,7);
#endif
    PUSH_MPS("lseu",28,89,1120,834.68,7);
    PUSH_MPS("misc03",96,160,3360,1910.,7);
    PUSH_MPS("misc06",820,1808,12850.8607,12841.6,7);
#if HOWMANY>1
    PUSH_MPS("misc07",212,260,2810,1415.0,7);
    PUSH_MPS("mitre",2054,10724,115155,114740.5184,7);
#endif
    PUSH_MPS("mod008",6,319,307,290.9,7);
    PUSH_MPS("mod010",146,2655,6548,6532.08,7);
#if HOWMANY==3
    PUSH_MPS("mod011",4480,10958,-54558535,-62121982.55,7);
    PUSH_MPS("modglob",291,422,20740508,20430947.,7);
    PUSH_MPS("noswot",182,128,-43,-43.0,7);
#endif
#if HOWMANY>1
    PUSH_MPS("nw04",36,87482,16862,16310.66667,7);
#endif
    PUSH_MPS("p0033",16,33,3089,2520.57,7);
    PUSH_MPS("p0201",133,201,7615,6875.0,7);
    PUSH_MPS("p0282",241,282,258411,176867.50,7);
    PUSH_MPS("p0548",176,548,8691,315.29,7);
    PUSH_MPS("p2756",755,2756,3124,2688.75,7);
#if HOWMANY==3
    PUSH_MPS("pk1",45,86,11.0,0.0,7);
#endif
#if HOWMANY>1
    PUSH_MPS("pp08a",136,240,7350.0,2748.3452381,7);
    PUSH_MPS("pp08aCUTS",246,240,7350.0,5480.6061563,7);
#endif
#if HOWMANY==3
    PUSH_MPS("qiu",1192,840,-132.873137,-931.638857,7);
#endif
    PUSH_MPS("qnet1",503,1541,16029.692681,14274.102667,7);
    PUSH_MPS("qnet1_o",456,1541,16029.692681,12095.571667,7);
    PUSH_MPS("rentacar",6803,9557,30356761,28806137.644,7);
    PUSH_MPS("rgn",24,180,82.1999,48.7999,7);
#if HOWMANY==3
    PUSH_MPS("rout",291,556,1077.56,981.86428571,7);
    PUSH_MPS("set1ch",492,712,54537.75,32007.73,7);
#endif
    //    PUSH_MPS("seymour",4944,1372,423,403.84647413,7);
    PUSH_MPS("stein27",118,27,18,13.0,7);
#if HOWMANY>1
    PUSH_MPS("stein45",331,45,30,22.0,7);
#endif
    PUSH_MPS("vpm1",234,378,20,15.4167,7);
    PUSH_MPS("vpm2",234,378,13.75,9.8892645972,7);
#endif
  }
#undef PUSH_MPS
    
  int numProbSolved = 0;
  double timeTaken=0.0;
  
  /*
    Open the main loop to step through the MPS problems.
  */
  for (m = 0 ; m < mpsName.size() ; m++) {
    std::cout << "  processing mps file: " << mpsName[m] 
              << " (" << m+1 << " out of " << mpsName.size() << ")\n";
    /*
      Stage 1: Read the MPS
      and make sure the size of the constraint matrix is correct.
    */
    CbcModel * model = new CbcModel(saveModel);
      
    std::string fn = miplibDir+mpsName[m] ;
    model->solver()->readMps(fn.c_str(),"") ;
    assert(model->getNumRows() == nRows[m]) ;
    assert(model->getNumCols() == nCols[m]) ;

    /*
      Stage 2: Call solver to solve the problem.
      then check the return code and objective.
    */

    double startTime = CoinCpuTime()+CoinCpuTimeJustChildren();
    if (model->getMaximumNodes()>200000) {
      model->setMaximumNodes(200000);
    }
    OsiClpSolverInterface * si =
      dynamic_cast<OsiClpSolverInterface *>(model->solver()) ;
    assert (si != NULL);
    // get clp itself
    ClpSimplex * modelC = si->getModelPtr();
    modelC->tightenPrimalBounds();
    model->initialSolve();
    if (modelC->dualBound()==1.0e10) {
      // user did not set - so modify
      // get largest scaled away from bound
      double largest=1.0e-12;
      int numberRows = modelC->numberRows();
      const double * rowPrimal = modelC->primalRowSolution();
      const double * rowLower = modelC->rowLower();
      const double * rowUpper = modelC->rowUpper();
      const double * rowScale = modelC->rowScale();
      int iRow;
      for (iRow=0;iRow<numberRows;iRow++) {
	double value = rowPrimal[iRow];
	double above = value-rowLower[iRow];
	double below = rowUpper[iRow]-value;
	if (rowScale) {
	  double multiplier = rowScale[iRow];
	  above *= multiplier;
	  below *= multiplier;
	}
	if (above<1.0e12) {
	  largest = CoinMax(largest,above);
	}
	if (below<1.0e12) {
	  largest = CoinMax(largest,below);
	}
      }
      
      int numberColumns = modelC->numberColumns();
      const double * columnPrimal = modelC->primalColumnSolution();
      const double * columnLower = modelC->columnLower();
      const double * columnUpper = modelC->columnUpper();
      const double * columnScale = modelC->columnScale();
      int iColumn;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double value = columnPrimal[iColumn];
	double above = value-columnLower[iColumn];
	double below = columnUpper[iColumn]-value;
	if (columnScale) {
	  double multiplier = 1.0/columnScale[iColumn];
	  above *= multiplier;
	  below *= multiplier;
	}
	if (above<1.0e12) {
	  largest = CoinMax(largest,above);
	}
	if (below<1.0e12) {
	  largest = CoinMax(largest,below);
	}
      }
      //std::cout<<"Largest (scaled) away from bound "<<largest<<std::endl;
      modelC->setDualBound(CoinMax(1.0001e8,
				   CoinMin(1000.0*largest,1.00001e10)));
    }
    model->setMinimumDrop(min(5.0e-2,
			      fabs(model->getMinimizationObjValue())*1.0e-3+1.0e-4));
    if (model->getNumCols()<500) {
      model->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
    } else if (model->getNumCols()<5000) {
      model->setMaximumCutPassesAtRoot(100); // use minimum drop
    } else {
      model->setMaximumCutPassesAtRoot(20);
    }
    // If defaults then increase trust for small models
    if (model->numberStrong()==5&&model->numberBeforeTrust()==10) {
      int numberColumns = model->getNumCols();
      if (numberColumns<=50) {
	model->setNumberBeforeTrust(1000);
      } else if (numberColumns<=100) {
	model->setNumberBeforeTrust(100);
      } else if (numberColumns<=300) {
	model->setNumberBeforeTrust(50);
      }
    }
    model->branchAndBound();
      
    double timeOfSolution = CoinCpuTime()+CoinCpuTimeJustChildren()-startTime;
    // Print more statistics
    std::cout<<"Cuts at root node changed objective from "<<model->getContinuousObjective()
	     <<" to "<<model->rootObjectiveAfterCuts()<<std::endl;
    int numberGenerators = model->numberCutGenerators();
    for (int iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CbcCutGenerator * generator = model->cutGenerator(iGenerator);
      std::cout<<generator->cutGeneratorName()<<" was tried "
	       <<generator->numberTimesEntered()<<" times and created "
	       <<generator->numberCutsInTotal()<<" cuts of which "
	       <<generator->numberCutsActive()<<" were active after adding rounds of cuts";
      if (generator->timing())
	std::cout<<" ( "<<generator->timeInCutGenerator()<<" seconds)"<<std::endl;
      else
	std::cout<<std::endl;
    }
    if (!model->status()) { 
      double soln = model->getObjValue();       
      CoinRelFltEq eq(1.0e-3) ;
      if (eq(soln,objValue[m])) { 
        std::cout 
          <<"cbc_clp"<<" "
          << soln << " = " << objValue[m] << " ; okay";
        numProbSolved++;
      } else  { 
        std::cout <<"cbc_clp" <<" " <<soln << " != " <<objValue[m]
		  << "; error=" << fabs(objValue[m] - soln); 
      }
    } else {
      std::cout << "error; too many nodes" ;
    }
    std::cout<<" - took " <<timeOfSolution<<" seconds.\n";
    timeTaken += timeOfSolution;
    delete model;
  }
  std::cout 
    <<"cbc_clp" 
    <<" solved " 
    <<numProbSolved
    <<" out of "
    <<objValue.size()
    <<" and took "
    <<timeTaken
    <<" seconds."
    <<std::endl;
}
