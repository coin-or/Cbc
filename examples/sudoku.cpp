// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CoinModel.hpp"
// For all different
#include "CbcBranchCut.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchAllDifferent.hpp"
#include "CbcCutGenerator.hpp"
#include "CglAllDifferent.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CglStored.hpp"

#include "CoinTime.hpp"

/************************************************************************

This shows how we can define a new branching method to solve problems with
all different constraints.

We are going to solve a sudoku problem such as

1, , ,4, , ,7, ,
 ,2, , ,5, , ,8,
8,7,3, , ,6, , ,9
4, , ,7, , ,1, ,
 ,5, , ,8, , ,2,
 , ,6, ,4,9, , ,3
7, , ,1, , ,4, ,
 ,8, ,6,2, , ,5,
 , ,9, ,7,3, ,1,6

The input should be exported from spreadsheet as a csv file where cells are
empty unless value is known.

We set up a fake objective and simple constraints to say sum must be 45

and then we add all different branching (CbcBranchAllDifferent) 
and all different cuts (to fix variables) (CglAllDifferent).

CbcBranchAllDifferent is really just an example of a cut branch.  If we wish to stop x==y
then we can have two branches - one x <= y-1 and the other x >= y+1.  It should be easy for the user to 
make up similar cut branches for other uses.

Note - this is all we need to solve most 9 x 9 puzzles because they seem to solve
at root node or easily.  To solve 16 x 16 puzzles we need more.  All different cuts
need general integer variables to be fixed while we can branch so they are just at bounds.  
To get round that we can introduce extra 0-1 variables such that general integer x = sum j * delta j
and then do N way branching on these (CbcNWay) so that we fix one delta j to 1.  At the same time we use the
new class CbcConsequence (used in N way branching) which when delta j goes to 1 fixes other variables.
So it will fix x to the correct value and while we are at it we can fix some delta variables in other
sets to zero (as per all different rules).  Finally as well as storing the instructions which say if 
delta 11 is 1 then delta 21 is 0 we can also add this in as a cut using new trivial cut class
CglStored.

************************************************************************/

int main(int argc, const char *argv[])
{
  // Get data
  std::string fileName = "./sudoku_sample.csv";
  if (argc >= 2)
    fileName = argv[1];
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    printf("Unable to open file %s\n", fileName.c_str());
    exit(0);
  }
#define MAX_SIZE 16
  int valueOffset = 1;
  double lo[MAX_SIZE * MAX_SIZE], up[MAX_SIZE * MAX_SIZE];
  char line[80];
  int row, column;
  /*************************************** 
     Read .csv file and see if 9 or 16 Su Doku
  ***************************************/
  int size = 9;
  for (row = 0; row < size; row++) {
    fgets(line, 80, fp);
    // Get size of sudoku puzzle (9 or 16)
    if (!row) {
      int get = 0;
      size = 1;
      while (line[get] >= 32) {
        if (line[get] == ',')
          size++;
        get++;
      }
      assert(size == 9 || size == 16);
      printf("Solving Su Doku of size %d\n", size);
      if (size == 16)
        valueOffset = 0;
    }
    int get = 0;
    for (column = 0; column < size; column++) {
      lo[size * row + column] = valueOffset;
      up[size * row + column] = valueOffset - 1 + size;
      if (line[get] != ',' && line[get] >= 32) {
        // skip blanks
        if (line[get] == ' ') {
          get++;
          continue;
        }
        int value = line[get] - '0';
        if (size == 9) {
          assert(value >= 1 && value <= 9);
        } else {
          assert(size == 16);
          if (value < 0 || value > 9) {
            if (line[get] == '"') {
              get++;
              value = 10 + line[get] - 'A';
              if (value < 10 || value > 15) {
                value = 10 + line[get] - 'a';
              }
              get++;
            } else {
              value = 10 + line[get] - 'A';
              if (value < 10 || value > 15) {
                value = 10 + line[get] - 'a';
              }
            }
          }
          assert(value >= 0 && value <= 15);
        }
        lo[size * row + column] = value;
        up[size * row + column] = value;
        get++;
      }
      get++;
    }
  }
  int block_size = (int)sqrt((double)size);
  /*************************************** 
     Now build rules for all different
     3*9 or 3*16 sets of variables
     Number variables by row*size+column
  ***************************************/
  int starts[3 * MAX_SIZE + 1];
  int which[3 * MAX_SIZE * MAX_SIZE];
  int put = 0;
  int set = 0;
  starts[0] = 0;
  // By row
  for (row = 0; row < size; row++) {
    for (column = 0; column < size; column++)
      which[put++] = row * size + column;
    starts[set + 1] = put;
    set++;
  }
  // By column
  for (column = 0; column < size; column++) {
    for (row = 0; row < size; row++)
      which[put++] = row * size + column;
    starts[set + 1] = put;
    set++;
  }
  // By box
  for (row = 0; row < size; row += block_size) {
    for (column = 0; column < size; column += block_size) {
      for (int row2 = row; row2 < row + block_size; row2++) {
        for (int column2 = column; column2 < column + block_size; column2++)
          which[put++] = row2 * size + column2;
      }
      starts[set + 1] = put;
      set++;
    }
  }
  OsiClpSolverInterface solver1;

  /*************************************** 
     Create model
     Set variables to be general integer variables although
     priorities probably mean that it won't matter
  ***************************************/
  CoinModel build;
  // Columns
  char name[4];
  for (row = 0; row < size; row++) {
    for (column = 0; column < size; column++) {
      if (row < 10) {
        if (column < 10)
          sprintf(name, "X%d%d", row, column);
        else
          sprintf(name, "X%d%c", row, 'A' + (column - 10));
      } else {
        if (column < 10)
          sprintf(name, "X%c%d", 'A' + (row - 10), column);
        else
          sprintf(name, "X%c%c", 'A' + (row - 10), 'A' + (column - 10));
      }
      double value = CoinDrand48() * 100.0;
      build.addColumn(0, NULL, NULL, lo[size * row + column],
        up[size * row + column], value, name, true);
    }
  }
  /*************************************** 
     Now add in extra variables for N way branching
  ***************************************/
  for (row = 0; row < size; row++) {
    for (column = 0; column < size; column++) {
      int iColumn = size * row + column;
      double value = lo[iColumn];
      if (value < up[iColumn]) {
        for (int i = 0; i < size; i++)
          build.addColumn(0, NULL, NULL, 0.0, 1.0, 0.0);
      } else {
        // fixed
        // obviously could do better  if we missed out variables
        int which = ((int)value) - valueOffset;
        for (int i = 0; i < size; i++) {
          if (i != which)
            build.addColumn(0, NULL, NULL, 0.0, 0.0, 0.0);
          else
            build.addColumn(0, NULL, NULL, 0.0, 1.0, 0.0);
        }
      }
    }
  }

  /*************************************** 
     Now rows
  ***************************************/
  double values[] = { 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0 };
  int indices[MAX_SIZE + 1];
  double rhs = size == 9 ? 45.0 : 120.0;
  for (row = 0; row < 3 * size; row++) {
    int iStart = starts[row];
    for (column = 0; column < size; column++)
      indices[column] = which[column + iStart];
    build.addRow(size, indices, values, rhs, rhs);
  }
  double values2[MAX_SIZE + 1];
  values2[0] = -1.0;
  for (row = 0; row < size; row++)
    values2[row + 1] = row + valueOffset;
  // Now add rows for extra variables
  for (row = 0; row < size; row++) {
    for (column = 0; column < size; column++) {
      int iColumn = row * size + column;
      int base = size * size + iColumn * size;
      indices[0] = iColumn;
      for (int i = 0; i < size; i++)
        indices[i + 1] = base + i;
      build.addRow(size + 1, indices, values2, 0.0, 0.0);
    }
  }
  solver1.loadFromCoinModel(build);
  build.writeMps("xx.mps");

  double time1 = CoinCpuTime();
  solver1.initialSolve();
  CbcModel model(solver1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  model.solver()->setHintParam(OsiDoScale, false, OsiHintTry);
  /*************************************** 
    Add in All different cut generator and All different branching
    So we will have integers then cut branching then N way branching
    in reverse priority order
  ***************************************/

  // Cut generator
  CglAllDifferent allDifferent(3 * size, starts, which);
  model.addCutGenerator(&allDifferent, -99, "allDifferent");
  model.cutGenerator(0)->setWhatDepth(5);
  CbcObject **objects = new CbcObject *[4 * size * size];
  int nObj = 0;
  for (row = 0; row < 3 * size; row++) {
    int iStart = starts[row];
    objects[row] = new CbcBranchAllDifferent(&model, size, which + iStart);
    objects[row]->setPriority(2000 + nObj); // do after rest satisfied
    nObj++;
  }
  /*************************************** 
     Add in N way branching and while we are at it add in cuts
  ***************************************/
  CglStored stored;
  for (row = 0; row < size; row++) {
    for (column = 0; column < size; column++) {
      int iColumn = row * size + column;
      int base = size * size + iColumn * size;
      int i;
      for (i = 0; i < size; i++)
        indices[i] = base + i;
      CbcNWay *obj = new CbcNWay(&model, size, indices, nObj);
      int seq[200];
      int newUpper[200];
      memset(newUpper, 0, sizeof(newUpper));
      for (i = 0; i < size; i++) {
        int state = 9999;
        int one = 1;
        int nFix = 1;
        // Fix real variable
        seq[0] = iColumn;
        newUpper[0] = valueOffset + i;
        int kColumn = base + i;
        int j;
        // same row
        for (j = 0; j < size; j++) {
          int jColumn = row * size + j;
          int jjColumn = size * size + jColumn * size + i;
          if (jjColumn != kColumn) {
            seq[nFix++] = jjColumn; // must be zero
          }
        }
        // same column
        for (j = 0; j < size; j++) {
          int jColumn = j * size + column;
          int jjColumn = size * size + jColumn * size + i;
          if (jjColumn != kColumn) {
            seq[nFix++] = jjColumn; // must be zero
          }
        }
        // same block
        int kRow = row / block_size;
        kRow *= block_size;
        int kCol = column / block_size;
        kCol *= block_size;
        for (j = kRow; j < kRow + block_size; j++) {
          for (int jc = kCol; jc < kCol + block_size; jc++) {
            int jColumn = j * size + jc;
            int jjColumn = size * size + jColumn * size + i;
            if (jjColumn != kColumn) {
              seq[nFix++] = jjColumn; // must be zero
            }
          }
        }
        // seem to need following?
        const int *upperAddress = newUpper;
        const int *seqAddress = seq;
        CbcFixVariable fix(1, &state, &one, &upperAddress, &seqAddress, &nFix, &upperAddress, &seqAddress);
        obj->setConsequence(indices[i], fix);
        // Now do as cuts
        for (int kk = 1; kk < nFix; kk++) {
          int jColumn = seq[kk];
          int cutInd[2];
          cutInd[0] = kColumn;
          if (jColumn > kColumn) {
            cutInd[1] = jColumn;
            stored.addCut(-COIN_DBL_MAX, 1.0, 2, cutInd, values);
          }
        }
      }
      objects[nObj] = obj;
      objects[nObj]->setPriority(nObj);
      nObj++;
    }
  }
  model.addObjects(nObj, objects);
  for (row = 0; row < nObj; row++)
    delete objects[row];
  delete[] objects;
  model.messageHandler()->setLogLevel(1);
  model.addCutGenerator(&stored, 1, "stored");
  // Say we want timings
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }
  // Set this to get all solutions (all ones in newspapers should only have one)
  //model.setCutoffIncrement(-1.0e6);
  /*************************************** 
     Do branch and bound
  ***************************************/
  // Do complete search
  model.branchAndBound();
  std::cout << "took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  /*************************************** 
     Print solution and check it is feasible
     We could modify output so could be imported by spreadsheet
  ***************************************/
  if (model.getMinimizationObjValue() < 1.0e50) {

    const double *solution = model.bestSolution();
    int put = 0;
    for (row = 0; row < size; row++) {
      for (column = 0; column < size; column++) {
        int value = (int)floor(solution[row * size + column] + 0.5);
        assert(value >= lo[put] && value <= up[put]);
        // save for later test
        lo[put++] = value;
        printf("%d ", value);
      }
      printf("\n");
    }
    // check valid
    bool valid = true;
    // By row
    for (row = 0; row < size; row++) {
      put = 0;
      for (column = 0; column < size; column++)
        which[put++] = row * size + column;
      assert(put == size);
      int i;
      for (i = 0; i < put; i++)
        which[i] = (int)lo[which[i]];
      std::sort(which, which + put);
      int last = valueOffset - 1;
      for (i = 0; i < put; i++) {
        int value = which[i];
        if (value != last + 1)
          valid = false;
        last = value;
      }
    }
    // By column
    for (column = 0; column < size; column++) {
      put = 0;
      for (row = 0; row < size; row++)
        which[put++] = row * size + column;
      assert(put == size);
      int i;
      for (i = 0; i < put; i++)
        which[i] = (int)lo[which[i]];
      std::sort(which, which + put);
      int last = valueOffset - 1;
      for (i = 0; i < put; i++) {
        int value = which[i];
        if (value != last + 1)
          valid = false;
        last = value;
      }
    }
    // By box
    for (row = 0; row < size; row += block_size) {
      for (column = 0; column < size; column += block_size) {
        put = 0;
        for (int row2 = row; row2 < row + block_size; row2++) {
          for (int column2 = column; column2 < column + block_size; column2++)
            which[put++] = row2 * size + column2;
        }
        assert(put == size);
        int i;
        for (i = 0; i < put; i++)
          which[i] = (int)lo[which[i]];
        std::sort(which, which + put);
        int last = valueOffset - 1;
        for (i = 0; i < put; i++) {
          int value = which[i];
          if (value != last + 1)
            valid = false;
          last = value;
        }
      }
    }
    if (valid) {
      printf("solution is valid\n");
    } else {
      printf("solution is not valid\n");
      abort();
    }
  }
  return 0;
}
