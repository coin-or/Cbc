// Copyright (C) 2000-2009, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// UnitTest for CglGomory adapted for lift-and-project

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "CoinPragma.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglLandP.hpp"

void
CglLandPUnitTest(
    OsiSolverInterface * si,
    const std::string &mpsDir)
{
    CoinRelFltEq eq(1e-05);
    // Test default constructor
    {
        CglLandP aGenerator;
        assert(aGenerator.parameter().pivotLimit==20);
        assert(aGenerator.parameter().maxCutPerRound==5000);
        assert(aGenerator.parameter().failedPivotLimit==1);
        assert(aGenerator.parameter().degeneratePivotLimit==0);
        assert(eq(aGenerator.parameter().pivotTol, 1e-04));
        assert(eq(aGenerator.parameter().away, 5e-04));
        assert(eq(aGenerator.parameter().timeLimit, COIN_DBL_MAX));
        assert(eq(aGenerator.parameter().singleCutTimeLimit, COIN_DBL_MAX));
        assert(aGenerator.parameter().useTableauRow==true);
        assert(aGenerator.parameter().modularize==false);
        assert(aGenerator.parameter().strengthen==true);
        assert(aGenerator.parameter().perturb==true);
        assert(aGenerator.parameter().pivotSelection==CglLandP::mostNegativeRc);
    }


    // Test copy constructor
    {
        CglLandP a;
        {
            CglLandP b;
            b.parameter().pivotLimit = 100;
            b.parameter().maxCutPerRound = 100;
            b.parameter().failedPivotLimit = 10;
            b.parameter().degeneratePivotLimit = 10;
            b.parameter().pivotTol = 1e-07;
            b.parameter().away = 1e-10;
            b.parameter().timeLimit = 120;
            b.parameter().singleCutTimeLimit = 15;
            b.parameter().useTableauRow = true;
            b.parameter().modularize = true;
            b.parameter().strengthen = false;
            b.parameter().perturb = false;
            b.parameter().pivotSelection=CglLandP::bestPivot;
            //Test Copy
            CglLandP c(b);
            assert(c.parameter().pivotLimit == 100);
            assert(c.parameter().maxCutPerRound == 100);
            assert(c.parameter().failedPivotLimit == 10);
            assert(c.parameter().degeneratePivotLimit == 10);
            assert(c.parameter().pivotTol == 1e-07);
            assert(c.parameter().away == 1e-10);
            assert(c.parameter().timeLimit == 120);
            assert(c.parameter().singleCutTimeLimit == 15);
            assert(c.parameter().useTableauRow == true);
            assert(c.parameter().modularize == true);
            assert(c.parameter().strengthen == false);
            assert(c.parameter().perturb == false);
            assert(c.parameter().pivotSelection == CglLandP::bestPivot);
            a=b;
            assert(a.parameter().pivotLimit == 100);
            assert(a.parameter().maxCutPerRound == 100);
            assert(a.parameter().failedPivotLimit == 10);
            assert(a.parameter().degeneratePivotLimit == 10);
            assert(a.parameter().pivotTol == 1e-07);
            assert(a.parameter().away == 1e-10);
            assert(a.parameter().timeLimit == 120);
            assert(a.parameter().singleCutTimeLimit == 15);
            assert(a.parameter().useTableauRow == true);
            assert(a.parameter().modularize == true);
            assert(a.parameter().strengthen == false);
            assert(a.parameter().perturb == false);
            assert(a.parameter().pivotSelection == CglLandP::bestPivot);
        }
    }

    {
        //  Maximize  2 x2
        // s.t.
        //    2x1 +  2x2 <= 3
        //   -2x1 +  2x2 <= 1
        //    7x1 +  4x2 <= 8
        //   -7x1 +  4x2 <= 1
        //     x1, x2 >= 0 and x1, x2 integer
        // Slacks are s1, s2, s3, s4



        //Test that problem is correct
        // Optimal Basis is x1, x2, s3, s4 with tableau
        //    x1            0.25 s1  -0.25 s2             =  0.5
        //           x2     0.25 s1   0.25 s2             =  1
        //                 -2.75 s1   0.75 s2    s3       =  0.5
        //                  0.75 s1  -2.75 s2        s4   =  0.5
        // z=              -0.25 s1  -0.25 s2             =  -1
        // Gomory cut from variable x1 is x2 <= 0.5
        // Can be improved by first pivoting s2 in and s4 out, then s1 in and s3 out
        // to x2 <= 0.25
        {
            CoinBigIndex start[2] = {0,4};
            int length[2] = {4,4};
            int rows[8] = {0,1,2,3,0,1,2,3};
            double elements[8] = {2.0,-2.0,7.0,-7.0,2.0,2.0,4.0,4.0};
            CoinPackedMatrix  columnCopy(true,4,2,8,elements,rows,start,length);

            double rowLower[4]={-COIN_DBL_MAX,-COIN_DBL_MAX,
                                -COIN_DBL_MAX,-COIN_DBL_MAX};
            double rowUpper[4]={3.,1.,8.,1.};
            double colLower[2]={0.0,0.0};
            double colUpper[2]={1.0,1.0};
            double obj[2]={-1,-1};
            int intVar[2]={0,1};

            OsiSolverInterface  * siP = si->clone();
            siP->loadProblem(columnCopy, colLower, colUpper, obj, rowLower, rowUpper);
            siP->setInteger(intVar,2);
            CglLandP test;
            test.setLogLevel(2);
            test.parameter().sepSpace = CglLandP::Full;
            siP->resolve();
            // Test generateCuts method
            if (0) {
                OsiCuts cuts;
                test.generateCuts(*siP,cuts);
                cuts.printCuts();
                assert(cuts.sizeRowCuts()==1);
                OsiRowCut aCut = cuts.rowCut(0);
                assert(eq(aCut.lb(), -.0714286));
                CoinPackedVector row = aCut.row();
                if (row.getNumElements() == 1)
                {
                    assert(row.getIndices()[0]==1);
                    assert(eq(row.getElements()[0], -4*.0714286));
                }
                else if (row.getNumElements() == 2)
                {
                    assert(row.getIndices()[0]==0);
                    assert(eq(row.getElements()[0], 0.));
                    assert(row.getIndices()[1]==1);
                    assert(eq(row.getElements()[1], -1));
                }
                OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);

                siP->resolve();
            }
            if (0)
            {
                OsiCuts cuts;
                test.generateCuts(*siP,cuts);
                cuts.printCuts();
                assert(cuts.sizeRowCuts()==1);
                OsiRowCut aCut = cuts.rowCut(0);
                CoinPackedVector row = aCut.row();
                if (row.getNumElements() == 1)
                {
                    assert(row.getIndices()[0]==1);
                    assert(eq(row.getElements()[0], -1));
                }
                else if (row.getNumElements() == 2)
                {
                    assert(row.getIndices()[0]==0);
                    assert(eq(row.getElements()[0], 0.));
                    assert(row.getIndices()[1]==1);
                    assert(eq(row.getElements()[1], -1));
                }
                assert(eq(aCut.lb(), 0.));
                OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);

                siP->resolve();
            }
            delete siP;
        }
    }

    if (1)  //Test on p0033
    {
        // Setup
        OsiSolverInterface  * siP = si->clone();
        std::string fn(mpsDir+"p0033");
        siP->readMps(fn.c_str(),"mps");
        siP->activateRowCutDebugger("p0033");
        CglLandP test;

        // Solve the LP relaxation of the model and
        // print out ofv for sake of comparison
        siP->initialSolve();
        double lpRelaxBefore=siP->getObjValue();
        assert( eq(lpRelaxBefore, 2520.5717391304347) );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif

        OsiCuts cuts;

        // Test generateCuts method
        test.generateCuts(*siP,cuts);
        OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);

        siP->resolve();
        double lpRelaxAfter=siP->getObjValue();
        //assert( eq(lpRelaxAfter, 2592.1908295194507) );

        std::cout<<"Relaxation after "<<lpRelaxAfter<<std::endl;
        assert( lpRelaxAfter> 2840. );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
        printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
        assert( lpRelaxBefore < lpRelaxAfter );

        delete siP;
    }
    if (1)  //test again with modularization
    {
        // Setup
        OsiSolverInterface  * siP = si->clone();
        std::string fn(mpsDir+"p0033");
        siP->readMps(fn.c_str(),"mps");
        siP->activateRowCutDebugger("p0033");
        CglLandP test;
        test.parameter().modularize = true;
        // Solve the LP relaxation of the model and
        // print out ofv for sake of comparison
        siP->initialSolve();
        double lpRelaxBefore=siP->getObjValue();
        assert( eq(lpRelaxBefore, 2520.5717391304347) );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif

        OsiCuts cuts;

        // Test generateCuts method
        test.generateCuts(*siP,cuts);
        OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);

        siP->resolve();
        double lpRelaxAfter=siP->getObjValue();
        //assert( eq(lpRelaxAfter, 2592.1908295194507) );

        std::cout<<"Relaxation after "<<lpRelaxAfter<<std::endl;
        assert( lpRelaxAfter> 2840. );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
        printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
        assert( lpRelaxBefore < lpRelaxAfter );

        delete siP;
    }
    if (1)  //test again with alternate pivoting rule
    {
        // Setup
        OsiSolverInterface  * siP = si->clone();
        std::string fn(mpsDir+"p0033");
        siP->readMps(fn.c_str(),"mps");
        siP->activateRowCutDebugger("p0033");
        CglLandP test;
        test.parameter().pivotSelection = CglLandP::bestPivot;
        // Solve the LP relaxation of the model and
        // print out ofv for sake of comparison
        siP->initialSolve();
        double lpRelaxBefore=siP->getObjValue();
        assert( eq(lpRelaxBefore, 2520.5717391304347) );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif

        OsiCuts cuts;

        // Test generateCuts method
        test.generateCuts(*siP,cuts);
        OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);

        siP->resolve();
        double lpRelaxAfter=siP->getObjValue();
        //assert( eq(lpRelaxAfter, 2592.1908295194507) );

        std::cout<<"Relaxation after "<<lpRelaxAfter<<std::endl;
        assert( lpRelaxAfter> 2840. );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
        printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
        assert( lpRelaxBefore < lpRelaxAfter );

        delete siP;
    }

    if (1)  //Finally test code in documentation
    {
        // Setup
        OsiSolverInterface  * siP = si->clone();
        std::string fn(mpsDir+"p0033");
        siP->readMps(fn.c_str(),"mps");
        siP->activateRowCutDebugger("p0033");
        CglLandP landpGen;

        landpGen.parameter().timeLimit = 10.;
        landpGen.parameter().pivotLimit = 2;


        // Solve the LP relaxation of the model and
        // print out ofv for sake of comparison
        siP->initialSolve();
        double lpRelaxBefore=siP->getObjValue();
        assert( eq(lpRelaxBefore, 2520.5717391304347) );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif

        OsiCuts cuts;

        // Test generateCuts method
        landpGen.generateCuts(*siP, cuts);
        OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);

        siP->resolve();
        double lpRelaxAfter=siP->getObjValue();
        //assert( eq(lpRelaxAfter, 2592.1908295194507) );

        std::cout<<"Relaxation after "<<lpRelaxAfter<<std::endl;
        assert( lpRelaxAfter> 2840. );
#ifdef CGL_DEBUG
        printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
        printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
        assert( lpRelaxBefore < lpRelaxAfter );

        delete siP;
    }
}
