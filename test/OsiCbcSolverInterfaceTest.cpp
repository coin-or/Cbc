// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"

//#include <cassert>
//#include <cstdlib>
//#include <cstdio>
//#include <iostream>

#include "OsiCbcSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiUnitTests.hpp"
#include "CoinMessage.hpp"
#include "CoinModel.hpp"

//#############################################################################

namespace {
CoinPackedMatrix &BuildExmip1Mtx()
/*
  Simple function to build a packed matrix for the exmip1 example used in
  tests. The function exists solely to hide the intermediate variables.
  Probably could be written as an initialised declaration.
  See COIN/Mps/Sample/exmip1.mps for a human-readable presentation.

  Ordered triples seem easiest. They're listed in row-major order.
*/

{
  int rowndxs[] = { 0, 0, 0, 0, 0,
    1, 1,
    2, 2,
    3, 3,
    4, 4, 4 };
  int colndxs[] = { 0, 1, 3, 4, 7,
    1, 2,
    2, 5,
    3, 6,
    0, 4, 7 };
  double coeffs[] = { 3.0, 1.0, -2.0, -1.0, -1.0,
    2.0, 1.1,
    1.0, 1.0,
    2.8, -1.2,
    5.6, 1.0, 1.9 };

  static CoinPackedMatrix exmip1mtx = CoinPackedMatrix(true, &rowndxs[0], &colndxs[0], &coeffs[0], 14);

  return (exmip1mtx);
}
}

//--------------------------------------------------------------------------
// test solution methods.
void OsiCbcSolverInterfaceUnitTest(const std::string &mpsDir, const std::string &netlibDir)
{
  {
    CoinRelFltEq eq;
    OsiCbcSolverInterface m;
    std::string fn = mpsDir + "exmip1";
    m.readMps(fn.c_str(), "mps");

    {
      OsiCbcSolverInterface im;
      OSIUNITTEST_ASSERT_ERROR(im.getNumCols() == 0, {}, "cbc", "default constructor");
      OSIUNITTEST_ASSERT_ERROR(im.getModelPtr() != NULL, {}, "cbc", "default constructor");
    }

    // Test copy constructor and assignment operator
    {
      OsiCbcSolverInterface lhs;
      {
        OsiCbcSolverInterface im(m);

        OsiCbcSolverInterface imC1(im);
        OSIUNITTEST_ASSERT_ERROR(imC1.getModelPtr() != im.getModelPtr(), {}, "cbc", "copy constructor");
        OSIUNITTEST_ASSERT_ERROR(imC1.getNumCols() == im.getNumCols(), {}, "cbc", "copy constructor");
        OSIUNITTEST_ASSERT_ERROR(imC1.getNumRows() == im.getNumRows(), {}, "cbc", "copy constructor");

        OsiCbcSolverInterface imC2(im);
        OSIUNITTEST_ASSERT_ERROR(imC2.getModelPtr() != im.getModelPtr(), {}, "cbc", "copy constructor");
        OSIUNITTEST_ASSERT_ERROR(imC2.getNumCols() == im.getNumCols(), {}, "cbc", "copy constructor");
        OSIUNITTEST_ASSERT_ERROR(imC2.getNumRows() == im.getNumRows(), {}, "cbc", "copy constructor");

        OSIUNITTEST_ASSERT_ERROR(imC1.getModelPtr() != imC2.getModelPtr(), {}, "cbc", "copy constructor");

        lhs = imC2;
      }

      // Test that lhs has correct values even though rhs has gone out of scope
      OSIUNITTEST_ASSERT_ERROR(lhs.getModelPtr() != m.getModelPtr(), {}, "cbc", "assignment operator");
      OSIUNITTEST_ASSERT_ERROR(lhs.getNumCols() == m.getNumCols(), {}, "cbc", "copy constructor");
      OSIUNITTEST_ASSERT_ERROR(lhs.getNumRows() == m.getNumRows(), {}, "cbc", "copy constructor");
    }

    // Test clone
    {
      OsiCbcSolverInterface cbcSi(m);
      OsiSolverInterface *siPtr = &cbcSi;
      OsiSolverInterface *siClone = siPtr->clone();
      OsiCbcSolverInterface *cbcClone = dynamic_cast< OsiCbcSolverInterface * >(siClone);

      OSIUNITTEST_ASSERT_ERROR(cbcClone != NULL, {}, "cbc", "clone");
      OSIUNITTEST_ASSERT_ERROR(cbcClone->getModelPtr() != cbcSi.getModelPtr(), {}, "cbc", "clone");
      OSIUNITTEST_ASSERT_ERROR(cbcClone->getNumRows() == cbcSi.getNumRows(), {}, "cbc", "clone");
      OSIUNITTEST_ASSERT_ERROR(cbcClone->getNumCols() == m.getNumCols(), {}, "cbc", "clone");

      delete siClone;
    }

    // test infinity
    {
      OsiCbcSolverInterface si;
      OSIUNITTEST_ASSERT_ERROR(si.getInfinity() == OsiCbcInfinity, {}, "cbc", "infinity");
    }

    // Test some catches
    if (!OsiCbcHasNDEBUG()) {
      OsiCbcSolverInterface solver;
      try {
        solver.setObjCoeff(0, 0.0);
        OSIUNITTEST_ADD_OUTCOME("cbc", "setObjCoeff on empty model", "should throw exception", OsiUnitTest::TestOutcome::ERROR, false);
      } catch (CoinError e) {
        if (OsiUnitTest::verbosity >= 1)
          std::cout << "Correct throw from setObjCoeff on empty model" << std::endl;
      }

      std::string fn = mpsDir + "exmip1";
      solver.readMps(fn.c_str(), "mps");
      OSIUNITTEST_CATCH_ERROR(solver.setObjCoeff(0, 0.0), {}, "cbc", "setObjCoeff on nonempty model");

      try {
        int index[] = { 0, 20 };
        double value[] = { 0.0, 0.0, 0.0, 0.0 };
        solver.setColSetBounds(index, index + 2, value);
        OSIUNITTEST_ADD_OUTCOME("cbc", "setColSetBounds on cols not in model", "should throw exception", OsiUnitTest::TestOutcome::ERROR, false);
      } catch (CoinError e) {
        if (OsiUnitTest::verbosity >= 1)
          std::cout << "Correct throw from setObjCoeff on empty model" << std::endl;
      }
    }

    {
      OsiCbcSolverInterface cbcSi(m);
      int nc = cbcSi.getNumCols();
      int nr = cbcSi.getNumRows();
      const double *cl = cbcSi.getColLower();
      const double *cu = cbcSi.getColUpper();
      const double *rl = cbcSi.getRowLower();
      const double *ru = cbcSi.getRowUpper();
      OSIUNITTEST_ASSERT_ERROR(nc == 8, return, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(nr == 5, return, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cl[0], 2.5), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cl[1], 0.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cu[1], 4.1), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cu[2], 1.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(rl[0], 2.5), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(rl[4], 3.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(ru[1], 2.1), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(ru[4], 15.), {}, "cbc", "read and copy exmip1");

      const double *cs = cbcSi.getColSolution();
      OSIUNITTEST_ASSERT_ERROR(eq(cs[0], 2.5), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cs[7], 0.0), {}, "cbc", "read and copy exmip1");

      OSIUNITTEST_ASSERT_ERROR(!eq(cl[3], 1.2345), {}, "cbc", "set col lower");
      cbcSi.setColLower(3, 1.2345);
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getColLower()[3], 1.2345), {}, "cbc", "set col lower");

      OSIUNITTEST_ASSERT_ERROR(!eq(cbcSi.getColUpper()[4], 10.2345), {}, "cbc", "set col upper");
      cbcSi.setColUpper(4, 10.2345);
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getColUpper()[4], 10.2345), {}, "cbc", "set col upper");

      // LH: Objective will depend on how underlying solver constructs and maintains initial solution
      double objValue = cbcSi.getObjValue();
      OSIUNITTEST_ASSERT_ERROR(eq(objValue, 3.5) || eq(objValue, 10.5), {}, "cbc", "getObjValue() before solve");

      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[0], 1.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[1], 0.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[2], 0.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[3], 0.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[4], 2.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[5], 0.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[6], 0.0), {}, "cbc", "read and copy exmip1");
      OSIUNITTEST_ASSERT_ERROR(eq(cbcSi.getObjCoefficients()[7], -1.0), {}, "cbc", "read and copy exmip1");
    }

    // Test matrixByRow method
    {
      const OsiCbcSolverInterface si(m);
      const CoinPackedMatrix *smP = si.getMatrixByRow();

      OSIUNITTEST_ASSERT_ERROR(smP->getMajorDim() == 5, return, "cbc", "getMatrixByRow: major dim");
      OSIUNITTEST_ASSERT_ERROR(smP->getMinorDim() == 8, return, "cbc", "getMatrixByRow: major dim");
      OSIUNITTEST_ASSERT_ERROR(smP->getNumElements() == 14, return, "cbc", "getMatrixByRow: num elements");
      OSIUNITTEST_ASSERT_ERROR(smP->getSizeVectorStarts() == 6, return, "cbc", "getMatrixByRow: num elements");

#ifdef OSICBC_TEST_MTX_STRUCTURE
      CoinRelFltEq eq;
      const double *ev = smP->getElements();
      OSIUNITTEST_ASSERT_ERROR(eq(ev[0], 3.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[1], 1.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[2], -2.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[3], -1.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[4], -1.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[5], 2.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[6], 1.1), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[7], 1.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[8], 1.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[9], 2.8), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[10], -1.2), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[11], 5.6), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[12], 1.0), {}, "cbc", "getMatrixByRow: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[13], 1.9), {}, "cbc", "getMatrixByRow: elements");

      const int *mi = smP->getVectorStarts();
      OSIUNITTEST_ASSERT_ERROR(mi[0] == 0, {}, "cbc", "getMatrixByRow: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[1] == 5, {}, "cbc", "getMatrixByRow: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[2] == 7, {}, "cbc", "getMatrixByRow: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[3] == 9, {}, "cbc", "getMatrixByRow: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[4] == 11, {}, "cbc", "getMatrixByRow: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[5] == 14, {}, "cbc", "getMatrixByRow: vector starts");

      const int *ei = smP->getIndices();
      OSIUNITTEST_ASSERT_ERROR(ei[0] == 0, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[1] == 1, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[2] == 3, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[3] == 4, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[4] == 7, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[5] == 1, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[6] == 2, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[7] == 2, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[8] == 5, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[9] == 3, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[10] == 6, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[11] == 0, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[12] == 4, {}, "cbc", "getMatrixByRow: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[13] == 7, {}, "cbc", "getMatrixByRow: indices");
#else // OSICBC_TEST_MTX_STRUCTURE

      CoinPackedMatrix exmip1Mtx;
      exmip1Mtx.reverseOrderedCopyOf(BuildExmip1Mtx());
      OSIUNITTEST_ASSERT_ERROR(exmip1Mtx.isEquivalent(*smP), {}, "cbc", "getMatrixByRow");
#endif // OSICBC_TEST_MTX_STRUCTURE
    }

    // Test adding several cuts, and handling of a coefficient of infinity
    // in the constraint matrix.
    {
      OsiCbcSolverInterface fim;
      std::string fn = mpsDir + "exmip1";
      fim.readMps(fn.c_str(), "mps");
      // exmip1.mps has 2 integer variables with index 2 & 3
      fim.initialSolve();
      OsiRowCut cuts[3];

      // Generate one ineffective cut plus two trivial cuts
      int c;
      int nc = fim.getNumCols();
      int *inx = new int[nc];
      for (c = 0; c < nc; c++)
        inx[c] = c;
      double *el = new double[nc];
      for (c = 0; c < nc; c++)
        el[c] = 1.0e-50 + ((double)c) * ((double)c);

      cuts[0].setRow(nc, inx, el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      el[4] = 0.0; // to get inf later

      for (c = 2; c < 4; c++) {
        el[0] = 1.0;
        inx[0] = c;
        cuts[c - 1].setRow(1, inx, el);
        cuts[c - 1].setLb(1.);
        cuts[c - 1].setUb(100.);
        cuts[c - 1].setEffectiveness(c);
      }
      fim.writeMps("x1.mps");
      fim.applyRowCuts(3, cuts);
      fim.writeMps("x2.mps");
      // resolve - should get message about zero elements
      fim.resolve();
      fim.writeMps("x3.mps");
      // check integer solution
      const double *cs = fim.getColSolution();
      CoinRelFltEq eq;
      OSIUNITTEST_ASSERT_ERROR(eq(cs[2], 1.0), {}, "cbc", "add cuts");
      OSIUNITTEST_ASSERT_ERROR(eq(cs[3], 1.0), {}, "cbc", "add cuts");
      // check will find invalid matrix
      el[0] = 1.0 / el[4];
      inx[0] = 0;
      cuts[0].setRow(nc, inx, el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      fim.applyRowCut(cuts[0]);
      // resolve - should get message about zero elements
      fim.resolve();
      OSIUNITTEST_ASSERT_WARNING(fim.isAbandoned(), {}, "cbc", "add cuts");
      delete[] el;
      delete[] inx;
    }

    // Test matrixByCol method
    {
      const OsiCbcSolverInterface si(m);
      const CoinPackedMatrix *smP = si.getMatrixByCol();

      OSIUNITTEST_ASSERT_ERROR(smP->getMajorDim() == 8, return, "cbc", "getMatrixByCol: major dim");
      OSIUNITTEST_ASSERT_ERROR(smP->getMinorDim() == 5, return, "cbc", "getMatrixByCol: minor dim");
      OSIUNITTEST_ASSERT_ERROR(smP->getNumElements() == 14, return, "cbc", "getMatrixByCol: number of elements");
      OSIUNITTEST_ASSERT_ERROR(smP->getSizeVectorStarts() == 9, return, "cbc", "getMatrixByCol: vector starts size");

#ifdef OSICBC_TEST_MTX_STRUCTURE
      CoinRelFltEq eq;
      const double *ev = smP->getElements();
      OSIUNITTEST_ASSERT_ERROR(eq(ev[0], 3.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[1], 5.6), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[2], 1.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[3], 2.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[4], 1.1), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[5], 1.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[6], -2.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[7], 2.8), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[8], -1.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[9], 1.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[10], 1.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[11], -1.2), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[12], -1.0), {}, "cbc", "getMatrixByCol: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[13], 1.9), {}, "cbc", "getMatrixByCol: elements");

      const CoinBigIndex *mi = smP->getVectorStarts();
      OSIUNITTEST_ASSERT_ERROR(mi[0] == 0, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[1] == 2, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[2] == 4, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[3] == 6, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[4] == 8, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[5] == 10, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[6] == 11, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[7] == 12, {}, "cbc", "getMatrixByCol: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[8] == 14, {}, "cbc", "getMatrixByCol: vector starts");

      const int *ei = smP->getIndices();
      OSIUNITTEST_ASSERT_ERROR(ei[0] == 0, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[1] == 4, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[2] == 0, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[3] == 1, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[4] == 1, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[5] == 2, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[6] == 0, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[7] == 3, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[8] == 0, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[9] == 4, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[10] == 2, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[11] == 3, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[12] == 0, {}, "cbc", "getMatrixByCol: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[13] == 4, {}, "cbc", "getMatrixByCol: indices");
#else // OSICBC_TEST_MTX_STRUCTURE

      CoinPackedMatrix &exmip1Mtx = BuildExmip1Mtx();
      OSIUNITTEST_ASSERT_ERROR(exmip1Mtx.isEquivalent(*smP), {}, "cbc", "getMatrixByCol");
#endif // OSICBC_TEST_MTX_STRUCTURE
    }

    //--------------
    // Test rowsense, rhs, rowrange, matrixByRow, solver assignment
    {
      OsiCbcSolverInterface lhs;
      {
        OsiCbcSolverInterface siC1(m);

        const char *siC1rs = siC1.getRowSense();
        OSIUNITTEST_ASSERT_ERROR(siC1rs[0] == 'G', {}, "cbc", "row sense");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[1] == 'L', {}, "cbc", "row sense");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[2] == 'E', {}, "cbc", "row sense");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[3] == 'R', {}, "cbc", "row sense");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[4] == 'R', {}, "cbc", "row sense");

        const double *siC1rhs = siC1.getRightHandSide();
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[0], 2.5), {}, "cbc", "right hand side");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[1], 2.1), {}, "cbc", "right hand side");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[2], 4.0), {}, "cbc", "right hand side");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[3], 5.0), {}, "cbc", "right hand side");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[4], 15.), {}, "cbc", "right hand side");

        const double *siC1rr = siC1.getRowRange();
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[0], 0.0), {}, "cbc", "row range");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[1], 0.0), {}, "cbc", "row range");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[2], 0.0), {}, "cbc", "row range");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[3], 5.0 - 1.8), {}, "cbc", "row range");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[4], 15.0 - 3.0), {}, "cbc", "row range");

        const CoinPackedMatrix *siC1mbr = siC1.getMatrixByRow();
        OSIUNITTEST_ASSERT_ERROR(siC1mbr != NULL, {}, "cbc", "matrix by row");
        OSIUNITTEST_ASSERT_ERROR(siC1mbr->getMajorDim() == 5, return, "cbc", "matrix by row: major dim");
        OSIUNITTEST_ASSERT_ERROR(siC1mbr->getMinorDim() == 8, return, "cbc", "matrix by row: major dim");
        OSIUNITTEST_ASSERT_ERROR(siC1mbr->getNumElements() == 14, return, "cbc", "matrix by row: num elements");
        OSIUNITTEST_ASSERT_ERROR(siC1mbr->getSizeVectorStarts() == 6, return, "cbc", "matrix by row: num elements");

#ifdef OSICBC_TEST_MTX_STRUCTURE
        const double *ev = siC1mbr->getElements();
        OSIUNITTEST_ASSERT_ERROR(eq(ev[0], 3.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[1], 1.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[2], -2.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[3], -1.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[4], -1.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[5], 2.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[6], 1.1), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[7], 1.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[8], 1.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[9], 2.8), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[10], -1.2), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[11], 5.6), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[12], 1.0), {}, "cbc", "matrix by row: elements");
        OSIUNITTEST_ASSERT_ERROR(eq(ev[13], 1.9), {}, "cbc", "matrix by row: elements");

        const CoinBigIndex *mi = siC1mbr->getVectorStarts();
        OSIUNITTEST_ASSERT_ERROR(mi[0] == 0, {}, "cbc", "matrix by row: vector starts");
        OSIUNITTEST_ASSERT_ERROR(mi[1] == 5, {}, "cbc", "matrix by row: vector starts");
        OSIUNITTEST_ASSERT_ERROR(mi[2] == 7, {}, "cbc", "matrix by row: vector starts");
        OSIUNITTEST_ASSERT_ERROR(mi[3] == 9, {}, "cbc", "matrix by row: vector starts");
        OSIUNITTEST_ASSERT_ERROR(mi[4] == 11, {}, "cbc", "matrix by row: vector starts");
        OSIUNITTEST_ASSERT_ERROR(mi[5] == 14, {}, "cbc", "matrix by row: vector starts");

        const int *ei = siC1mbr->getIndices();
        OSIUNITTEST_ASSERT_ERROR(ei[0] == 0, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[1] == 1, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[2] == 3, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[3] == 4, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[4] == 7, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[5] == 1, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[6] == 2, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[7] == 2, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[8] == 5, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[9] == 3, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[10] == 6, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[11] == 0, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[12] == 4, {}, "cbc", "matrix by row: indices");
        OSIUNITTEST_ASSERT_ERROR(ei[13] == 7, {}, "cbc", "matrix by row: indices");
#else // OSICBC_TEST_MTX_STRUCTURE

        CoinPackedMatrix exmip1Mtx;
        exmip1Mtx.reverseOrderedCopyOf(BuildExmip1Mtx());
        OSIUNITTEST_ASSERT_ERROR(exmip1Mtx.isEquivalent(*siC1mbr), {}, "cbc", "matrix by row");
#endif // OSICBC_TEST_MTX_STRUCTURE

        OSIUNITTEST_ASSERT_WARNING(siC1rs == siC1.getRowSense(), {}, "cbc", "row sense");
        OSIUNITTEST_ASSERT_WARNING(siC1rhs == siC1.getRightHandSide(), {}, "cbc", "right hand side");
        OSIUNITTEST_ASSERT_WARNING(siC1rr == siC1.getRowRange(), {}, "cbc", "row range");

        // Change CBC Model by adding free row
        OsiRowCut rc;
        rc.setLb(-COIN_DBL_MAX);
        rc.setUb(COIN_DBL_MAX);
        OsiCuts cuts;
        cuts.insert(rc);
        siC1.applyCuts(cuts);

        siC1rs = siC1.getRowSense();
        OSIUNITTEST_ASSERT_ERROR(siC1rs[0] == 'G', {}, "cbc", "row sense after adding row");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[1] == 'L', {}, "cbc", "row sense after adding row");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[2] == 'E', {}, "cbc", "row sense after adding row");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[3] == 'R', {}, "cbc", "row sense after adding row");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[4] == 'R', {}, "cbc", "row sense after adding row");
        OSIUNITTEST_ASSERT_ERROR(siC1rs[5] == 'N', {}, "cbc", "row sense after adding row");

        siC1rhs = siC1.getRightHandSide();
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[0], 2.5), {}, "cbc", "right hand side after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[1], 2.1), {}, "cbc", "right hand side after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[2], 4.0), {}, "cbc", "right hand side after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[3], 5.0), {}, "cbc", "right hand side after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[4], 15.), {}, "cbc", "right hand side after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rhs[5], 0.0), {}, "cbc", "right hand side after adding row");

        siC1rr = siC1.getRowRange();
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[0], 0.0), {}, "cbc", "row range after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[1], 0.0), {}, "cbc", "row range after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[2], 0.0), {}, "cbc", "row range after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[3], 5.0 - 1.8), {}, "cbc", "row range after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[4], 15.0 - 3.0), {}, "cbc", "row range after adding row");
        OSIUNITTEST_ASSERT_ERROR(eq(siC1rr[5], 0.0), {}, "cbc", "row range after adding row");

        lhs = siC1;
      }

      // Test that lhs has correct values even though siC1 has gone out of scope
      const char *lhsrs = lhs.getRowSense();
      OSIUNITTEST_ASSERT_ERROR(lhsrs[0] == 'G', {}, "cbc", "row sense after assignment");
      OSIUNITTEST_ASSERT_ERROR(lhsrs[1] == 'L', {}, "cbc", "row sense after assignment");
      OSIUNITTEST_ASSERT_ERROR(lhsrs[2] == 'E', {}, "cbc", "row sense after assignment");
      OSIUNITTEST_ASSERT_ERROR(lhsrs[3] == 'R', {}, "cbc", "row sense after assignment");
      OSIUNITTEST_ASSERT_ERROR(lhsrs[4] == 'R', {}, "cbc", "row sense after assignment");
      OSIUNITTEST_ASSERT_ERROR(lhsrs[5] == 'N', {}, "cbc", "row sense after assignment");

      const double *lhsrhs = lhs.getRightHandSide();
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrhs[0], 2.5), {}, "cbc", "right hand side after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrhs[1], 2.1), {}, "cbc", "right hand side after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrhs[2], 4.0), {}, "cbc", "right hand side after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrhs[3], 5.0), {}, "cbc", "right hand side after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrhs[4], 15.), {}, "cbc", "right hand side after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrhs[5], 0.0), {}, "cbc", "right hand side after assignment");

      const double *lhsrr = lhs.getRowRange();
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrr[0], 0.0), {}, "cbc", "row range after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrr[1], 0.0), {}, "cbc", "row range after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrr[2], 0.0), {}, "cbc", "row range after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrr[3], 5.0 - 1.8), {}, "cbc", "row range after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrr[4], 15.0 - 3.0), {}, "cbc", "row range after assignment");
      OSIUNITTEST_ASSERT_ERROR(eq(lhsrr[5], 0.0), {}, "cbc", "row range after assignment");

      const CoinPackedMatrix *lhsmbr = lhs.getMatrixByRow();
      OSIUNITTEST_ASSERT_ERROR(lhsmbr != NULL, {}, "cbc", "matrix by row after assignment");
      OSIUNITTEST_ASSERT_ERROR(lhsmbr->getMajorDim() == 6, return, "cbc", "matrix by row after assignment: major dim");
      OSIUNITTEST_ASSERT_ERROR(lhsmbr->getNumElements() == 14, return, "cbc", "matrix by row after assignment: num elements");

#ifdef OSICBC_TEST_MTX_STRUCTURE
      const double *ev = lhsmbr->getElements();
      OSIUNITTEST_ASSERT_ERROR(eq(ev[0], 3.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[1], 1.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[2], -2.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[3], -1.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[4], -1.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[5], 2.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[6], 1.1), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[7], 1.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[8], 1.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[9], 2.8), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[10], -1.2), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[11], 5.6), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[12], 1.0), {}, "cbc", "matrix by row after assignment: elements");
      OSIUNITTEST_ASSERT_ERROR(eq(ev[13], 1.9), {}, "cbc", "matrix by row after assignment: elements");

      const CoinBigIndex *mi = lhsmbr->getVectorStarts();
      OSIUNITTEST_ASSERT_ERROR(mi[0] == 0, {}, "cbc", "matrix by row after assignment: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[1] == 5, {}, "cbc", "matrix by row after assignment: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[2] == 7, {}, "cbc", "matrix by row after assignment: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[3] == 9, {}, "cbc", "matrix by row after assignment: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[4] == 11, {}, "cbc", "matrix by row after assignment: vector starts");
      OSIUNITTEST_ASSERT_ERROR(mi[5] == 14, {}, "cbc", "matrix by row after assignment: vector starts");

      const int *ei = lhsmbr->getIndices();
      OSIUNITTEST_ASSERT_ERROR(ei[0] == 0, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[1] == 1, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[2] == 3, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[3] == 4, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[4] == 7, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[5] == 1, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[6] == 2, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[7] == 2, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[8] == 5, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[9] == 3, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[10] == 6, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[11] == 0, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[12] == 4, {}, "cbc", "matrix by row after assignment: indices");
      OSIUNITTEST_ASSERT_ERROR(ei[13] == 7, {}, "cbc", "matrix by row after assignment: indices");
#else // OSICBC_TEST_MTX_STRUCTURE

      /*
  This admittedly looks bogus, but it's the equivalent operation on the matrix
  for inserting a cut of the form -Inf <= +Inf (i.e., a cut with no
  coefficients).
*/
      CoinPackedMatrix exmip1Mtx;
      exmip1Mtx.reverseOrderedCopyOf(BuildExmip1Mtx());
      CoinPackedVector freeRow;
      exmip1Mtx.appendRow(freeRow);
      OSIUNITTEST_ASSERT_ERROR(exmip1Mtx.isEquivalent(*lhsmbr), {}, "cbc", "matrix by row after assignment");
#endif // OSICBC_TEST_MTX_STRUCTURE
    }
  }

  // Test add/delete columns
  {
    OsiCbcSolverInterface m;
    std::string fn = mpsDir + "p0033";
    m.readMps(fn.c_str(), "mps");
    double inf = m.getInfinity();

    CoinPackedVector c0;
    c0.insert(0, 4);
    c0.insert(1, 1);
    m.addCol(c0, 0, inf, 3);
    m.initialSolve();
    double objValue = m.getObjValue();
    CoinRelFltEq eq(1.0e-2);
    OSIUNITTEST_ASSERT_ERROR(eq(objValue, 2520.57), {}, "cbc", "objvalue after adding col");

    // Try deleting first column that's nonbasic at lower bound (0).
    int *d = new int[1];
    CoinWarmStartBasis *cwsb = dynamic_cast< CoinWarmStartBasis * >(m.getWarmStart());
    OSIUNITTEST_ASSERT_ERROR(cwsb != NULL, {}, "cbc", "get warmstart basis");
    CoinWarmStartBasis::Status stati;
    int iCol;
    for (iCol = 0; iCol < cwsb->getNumStructural(); iCol++) {
      stati = cwsb->getStructStatus(iCol);
      if (stati == CoinWarmStartBasis::atLowerBound)
        break;
    }
    d[0] = iCol;
    m.deleteCols(1, d);
    delete[] d;
    delete cwsb;
    d = NULL;
    m.resolve();
    objValue = m.getObjValue();
    OSIUNITTEST_ASSERT_ERROR(eq(objValue, 2520.57), {}, "clp", "objvalue after deleting first col");

    // Try deleting column we added. If basic, go to initialSolve as deleting
    // basic variable trashes basis required for warm start.
    iCol = m.getNumCols() - 1;
    cwsb = dynamic_cast< CoinWarmStartBasis * >(m.getWarmStart());
    stati = cwsb->getStructStatus(iCol);
    delete cwsb;
    m.deleteCols(1, &iCol);
    if (stati == CoinWarmStartBasis::basic) {
      m.initialSolve();
    } else {
      m.resolve();
    }
    objValue = m.getObjValue();
    OSIUNITTEST_ASSERT_ERROR(eq(objValue, 2520.57), {}, "clp", "objvalue after deleting added col");
  }

  // Build a model
  {
    OsiCbcSolverInterface model;
    std::string fn = mpsDir + "p0033";
    model.readMps(fn.c_str(), "mps");
    // Point to data
    int numberRows = model.getNumRows();
    const double *rowLower = model.getRowLower();
    const double *rowUpper = model.getRowUpper();
    int numberColumns = model.getNumCols();
    const double *columnLower = model.getColLower();
    const double *columnUpper = model.getColUpper();
    const double *columnObjective = model.getObjCoefficients();
    // get row copy
    CoinPackedMatrix rowCopy = *model.getMatrixByRow();
    const int *column = rowCopy.getIndices();
    const int *rowLength = rowCopy.getVectorLengths();
    const CoinBigIndex *rowStart = rowCopy.getVectorStarts();
    const double *element = rowCopy.getElements();

    // solve
    model.initialSolve();
    // Now build new model
    CoinModel build;
    // Row bounds
    int iRow;
    for (iRow = 0; iRow < numberRows; iRow++) {
      build.setRowBounds(iRow, rowLower[iRow], rowUpper[iRow]);
    }
    // Column bounds and objective
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      build.setColumnLower(iColumn, columnLower[iColumn]);
      build.setColumnUpper(iColumn, columnUpper[iColumn]);
      build.setObjective(iColumn, columnObjective[iColumn]);
    }
    // Adds elements one by one by row (backwards by row)
    for (iRow = numberRows - 1; iRow >= 0; iRow--) {
      int start = rowStart[iRow];
      for (int j = start; j < start + rowLength[iRow]; j++)
        build(iRow, column[j], element[j]);
    }
    // Now create Model
    OsiCbcSolverInterface model2;
    model2.loadFromCoinModel(build);
    model2.initialSolve();
    // Save - should be continuous
    model2.writeMps("continuous");
    int *whichInteger = new int[numberColumns];
    for (iColumn = 0; iColumn < numberColumns; iColumn++)
      whichInteger[iColumn] = iColumn;
    // mark as integer
    model2.setInteger(whichInteger, numberColumns);
    delete[] whichInteger;
    // save - should be integer
    model2.writeMps("integer");

    // Now do with strings attached
    // Save build to show how to go over rows
    CoinModel saveBuild = build;
    build = CoinModel();
    // Column bounds
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      build.setColumnLower(iColumn, columnLower[iColumn]);
      build.setColumnUpper(iColumn, columnUpper[iColumn]);
    }
    // Objective - half the columns as is and half with multiplier of "1.0+multiplier"
    // Pick up from saveBuild (for no reason at all)
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = saveBuild.objective(iColumn);
      if (iColumn * 2 < numberColumns) {
        build.setObjective(iColumn, columnObjective[iColumn]);
      } else {
        // create as string
        char temp[100];
        sprintf(temp, "%g + abs(%g*multiplier)", value, value);
        build.setObjective(iColumn, temp);
      }
    }
    // It then adds rows one by one but for half the rows sets their values
    //      with multiplier of "1.0+1.5*multiplier"
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (iRow * 2 < numberRows) {
        // add row in simple way
        int start = rowStart[iRow];
        build.addRow(rowLength[iRow], column + start, element + start,
          rowLower[iRow], rowUpper[iRow]);
      } else {
        // As we have to add one by one let's get from saveBuild
        CoinModelLink triple = saveBuild.firstInRow(iRow);
        while (triple.column() >= 0) {
          int iColumn = triple.column();
          if (iColumn * 2 < numberColumns) {
            // just value as normal
            build(iRow, triple.column(), triple.value());
          } else {
            // create as string
            char temp[100];
            sprintf(temp, "%g + (1.5*%g*multiplier)", triple.value(), triple.value());
            build(iRow, iColumn, temp);
          }
          triple = saveBuild.next(triple);
        }
        // but remember to do rhs
        build.setRowLower(iRow, rowLower[iRow]);
        build.setRowUpper(iRow, rowUpper[iRow]);
      }
    }
    // If small switch on error printing
    if (numberColumns < 50)
      build.setLogLevel(1);
    // should fail as we never set multiplier
    OSIUNITTEST_ASSERT_ERROR(model2.loadFromCoinModel(build) != 0, {}, "cbc", "build model with missing multipliers");
    build.associateElement("multiplier", 0.0);
    OSIUNITTEST_ASSERT_ERROR(model2.loadFromCoinModel(build) == 0, {}, "cbc", "build model");
    model2.initialSolve();
    // It then loops with multiplier going from 0.0 to 2.0 in increments of 0.1
    for (double multiplier = 0.0; multiplier < 2.0; multiplier += 0.1) {
      build.associateElement("multiplier", multiplier);
      OSIUNITTEST_ASSERT_ERROR(model2.loadFromCoinModel(build, true) == 0, {}, "cbc", "build model with increasing multiplier");
      model2.resolve();
    }
  }

  // branch and bound
  {
    OsiCbcSolverInterface m;
    std::string fn = mpsDir + "p0033";
    m.readMps(fn.c_str(), "mps");
    m.initialSolve();
    //m.messageHandler()->setLogLevel(0);
    m.getModelPtr()->messageHandler()->setLogLevel(0);
    m.branchAndBound();
  }

  // branch and bound using CbcModel!!!!!!!
  {
    OsiCbcSolverInterface mm;
    OsiCbcSolverInterface m(&mm);
    std::string fn = mpsDir + "p0033";
    m.readMps(fn.c_str(), "mps");
    m.initialSolve();
    m.branchAndBound();
  }

  // Do common solverInterface testing
  {
    OsiCbcSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir, netlibDir);
  }
  {
    OsiCbcSolverInterface mm;
    OsiCbcSolverInterface m(&mm);
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir, netlibDir);
  }
}
