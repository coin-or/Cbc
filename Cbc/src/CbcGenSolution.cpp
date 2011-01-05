/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/


#include <string>
#include <cstdio>
#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "CoinFileIO.hpp"

#include "CbcGenCtlBlk.hpp"
#include "CbcGenParam.hpp"

namespace {

char svnid[] = "$Id: CbcGenSolution.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

namespace {

/*
  Helper routine to generate masks for selecting names to print. Returns true
  if masks are generated without error, false otherwise.

  This is John Forrest's code, shamelessly stolen from CoinSolve and tweaked
  just enough to allow it to be yanked out of the main body of CoinSolve.

  Returns the number of generated masks, -1 on error.
*/

int generateMasks (std::string proto, int longestName,
                   int *&maskStarts, char **&finalMasks)

{
    int nAst = 0 ;
    const char *pMask2 = proto.c_str() ;
    char pMask[100] ;
    int iChar ;
    int lengthMask = strlen(pMask2) ;
    assert (lengthMask < 100) ;

    int maxMasks = -1 ;

    maskStarts = 0 ;
    finalMasks = 0 ;
    /*
      Remove surrounding matched quotes if present. Abort if unmatched.
    */
    if (*pMask2 == '"') {
        if (pMask2[lengthMask-1] != '"') {
            printf("generateMasks: Mismatched \" in mask %s\n", pMask2) ;
            return (-1) ;
        } else {
            strcpy(pMask, pMask2 + 1) ;
            *strchr(pMask, '"') = '\0' ;
        }
    } else if (*pMask2 == '\'') {
        if (pMask2[lengthMask-1] != '\'') {
            printf("mismatched ' in mask %s\n", pMask2) ;
            return (maxMasks) ;
        } else {
            strcpy(pMask, pMask2 + 1) ;
            *strchr(pMask, '\'') = '\0' ;
        }
    } else {
        strcpy(pMask, pMask2) ;
    }
    /*
      Mask should not be longer than longest name.
    */
    if (lengthMask > longestName) {
        printf("mask %s too long - skipping\n", pMask) ;
        return (maxMasks) ;
    }
    /*
      Expand `*' to multiple masks with varying number of `?' characters.
    */
    maxMasks = 1 ;
    for (iChar = 0; iChar < lengthMask; iChar++) {
        if (pMask[iChar] == '*') {
            nAst++ ;
            maxMasks *= (longestName + 1) ;
        }
    }
    int nEntries = 1 ;
    maskStarts = new int[longestName+2] ;
    char ** masks = new char * [maxMasks] ;
    char ** newMasks = new char * [maxMasks] ;
    int i ;
    for (i = 0; i < maxMasks; i++) {
        masks[i] = new char[longestName+1] ;
        newMasks[i] = new char[longestName+1] ;
    }
    strcpy(masks[0], pMask) ;
    for (int iAst = 0; iAst < nAst; iAst++) {
        int nOldEntries = nEntries ;
        nEntries = 0 ;
        for (int iEntry = 0; iEntry < nOldEntries; iEntry++) {
            char * oldMask = masks[iEntry] ;
            char * ast = strchr(oldMask, '*') ;
            assert (ast) ;
            int length = strlen(oldMask) - 1 ;
            int nBefore = ast - oldMask ;
            int nAfter = length - nBefore ;
            // and add null
            nAfter++ ;
            for (int i = 0; i <= longestName - length; i++) {
                char * maskOut = newMasks[nEntries] ;
                memcpy(maskOut, oldMask, nBefore) ;
                for (int k = 0; k < i; k++)
                    maskOut[k+nBefore] = '?' ;
                memcpy(maskOut + nBefore + i, ast + 1, nAfter) ;
                nEntries++ ;
                assert (nEntries <= maxMasks) ;
            }
        }
        char ** temp = masks ;
        masks = newMasks ;
        newMasks = temp ;
    }
    /*
      Trim trailing blanks and record final length.
    */
    int * sort = new int[nEntries] ;
    for (i = 0; i < nEntries; i++) {
        char * maskThis = masks[i] ;
        int length = strlen(maskThis) ;
        while (maskThis[length-1] == ' ')
            length-- ;
        maskThis[length] = '\0' ;
        sort[i] = length ;
    }
    /*
      Sort by length.
    */
    CoinSort_2(sort, sort + nEntries, masks) ;
    int lastLength = -1 ;
    for (i = 0; i < nEntries; i++) {
        int length = sort[i] ;
        while (length > lastLength)
            maskStarts[++lastLength] = i ;
    }
    maskStarts[++lastLength] = nEntries ;
    delete [] sort ;
    for (i = 0; i < maxMasks; i++)
        delete [] newMasks[i] ;
    delete [] newMasks ;

    finalMasks = masks ;

    return (maxMasks) ;
}

/*
  Utility routine to check a string against the array of masks. Borrowed
  from CoinSolve.
*/

bool maskMatches (const int *starts, char **masks, const char *checkC)

{
    int length = strlen(checkC);
    while (checkC[length-1] == ' ')
        length--;
    for (int i = starts[length]; i < starts[length+1]; i++) {
        char * thisMask = masks[i];
        int k;
        for ( k = 0; k < length; k++) {
            if (thisMask[k] != '?' && thisMask[k] != checkC[k])
                break;
        }
        if (k == length)
            return true;
    }
    return (false) ;
}

}  // end unnamed namespace

/*
  Routine to write out the solution. Minimally adapted from John's code in
  CoinSolve to break out a few subroutines and use generic OSI functions.

  The print mode is established by the printingOptions command, and the integer
  coding is established by the order of declaration of the keyword parameters.
  As of 060920, known modes are

    normal (0)	  print nonzero primal variables
    integer (1)	  print nonzero integer primal variables
    special (2)	  print in a format suitable for OsiRowCutDebugger
    rows (3)	  `normal', plus rows with nonzero row activity
    all (4)	  all primal variables and row activities
*/

int CbcGenParamUtils::doSolutionParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;
    CbcModel *model = ctlBlk->model_ ;
    assert (model != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      It's hard to print a solution we don't have.
    */
    if (ctlBlk->bab_.haveAnswer_ == false) {
        std::cout << "There is no solution available to print." << std::endl ;
        return (retval) ;
    }
    OsiSolverInterface *osi = ctlBlk->bab_.answerSolver_ ;
    assert (osi != 0) ;
    /*
      Figure out where we're going to write the solution. As special cases,
      `$' says `use the previous output file' and `-' says `use stdout'.

      cbc will also accept `no string value' as stdout, but that'd be a real pain
      in this architecture.
    */
    std::string field = genParam->strVal() ;
    std::string fileName ;
    if (field == "$") {
        fileName = ctlBlk->lastSolnOut_ ;
        field = fileName ;
    }
    if (field == "-") {
        fileName = "stdout" ;
        field = fileName ;
    } else {
        fileName = field ;
    }
    if (!(fileName == "stdout" || fileName == "stderr")) {
        if (fileName[0] == '~') {
            char dirSep = CoinFindDirSeparator() ;
            if (fileName[1] == dirSep) {
                char *environVar = getenv("HOME") ;
                if (environVar) {
                    std::string home(environVar) ;
                    fileName = home + fileName.substr(1) ;
                }
            }
        }
        if (!(fileAbsPath(fileName) || fileName.substr(0, 2) == "./")) {
            fileName = ctlBlk->dfltDirectory_ + fileName ;
        }
    }
    /*
      See if we can open the file. Bail out if the open fails.
    */
    FILE *fp ;
    if (fileName == "stdout") {
        fp = stdout ;
    } else if (fileName == "stderr") {
        fp = stderr ;
    } else {
        fp = fopen(fileName.c_str(), "w") ;
        fp = fopen(fileName.c_str(), "w") ;
    }
    if (!fp) {
        std::cout
            << "Unable to open file `" << fileName
            << "', original name '" << genParam->strVal() << "'." << std::endl ;
        return (retval) ;
    } else {
        std::cout
            << "Writing solution to `" << fileName << "'." << std::endl ;
    }

    int m = osi->getNumRows() ;
    int n = osi->getNumCols() ;
    int iColumn ;
    const double *primalColSolution = osi->getColSolution() ;
    /*
      Separate printing a solution for humans from printing a solution for the
      row cut debugger. For the row cut debugger, we want to produce C++ code
      that can be pasted into the debugger's set of known problems.
    */
    if (ctlBlk->printMode_ == 2) {
        int k = 0 ;
        bool newLine = true ;
        bool comma = false ;
        fprintf(fp, "  int intIndicesV[] = {") ;
        for (iColumn = 0 ; iColumn < n ; iColumn++ ) {
            if (fabs(primalColSolution[iColumn]) > 0.5 && osi->isInteger(iColumn)) {
                if (newLine) {
                    fprintf(fp, "\n\t") ;
                    newLine = false ;
                } else {
                    fprintf(fp, ", ") ;
                }
                fprintf(fp, "%d", iColumn) ;
                if (++k == 10) {
                    k = 0 ;
                    newLine = true ;
                }
            }
        }
        fprintf(fp, "\n      } ;\n") ;
        k = 0 ;
        newLine = true ;
        fprintf(fp, "  double intSolnV[] = {") ;
        for (iColumn = 0 ; iColumn < n ; iColumn++) {
            double value = primalColSolution[iColumn] ;
            if (fabs(value) > 0.5 && osi->isInteger(iColumn)) {
                if (newLine) {
                    fprintf(fp, "\n\t") ;
                    newLine = false ;
                } else {
                    fprintf(fp, ", ") ;
                }
                if (value > 0) {
                    value = floor(value + .5)  ;
                } else {
                    value = ceil(value - .5) ;
                }
                int ivalue = static_cast<int>(value) ;
                fprintf(fp, "%d.0", ivalue) ;
                if (++k == 10) {
                    k = 0 ;
                    newLine = true ;
                }
            }
        }
        fprintf(fp, "\n      } ;\n") ;
        return (0) ;
    }
    /*
      Begin the code to generate output meant for a human.  What's our longest
      name? Scan the names we're going to print. printMode_ of 3 or 4 requires we
      scan the row names too. Force between 8 and 20 characters in any event.
    */
    int longestName = 0 ;
    for (int j = 0 ; j < n ; j++) {
        int len = osi->getColName(j).length() ;
        longestName = CoinMax(longestName, len) ;
    }
    if (ctlBlk->printMode_ >= 3) {
        for (int i = 0 ; i < m ; i++) {
            int len = osi->getRowName(i).length() ;
            longestName = CoinMax(longestName, len) ;
        }
    }
    /*
      Generate masks if we need to do so.
    */
    bool doMask = ctlBlk->printMask_ != "" ;
    int *maskStarts = NULL ;
    int maxMasks = 0 ;
    char **masks = NULL ;
    if (doMask) {
        maxMasks = generateMasks(ctlBlk->printMask_, longestName, maskStarts, masks) ;
        if (maxMasks < 0) {
            return (retval) ;
        }
    }
    /*
      Force the space allocated to names to be between 8 and 20 characters.
    */
    if (longestName < 8) {
        longestName = 8 ;
    } else if (longestName > 20) {
        longestName = 20 ;
    }
    /*
      Print requested components of the row solution. Only modes 3 (rows) and 4
      (all) will print row information. For the rows that we print, print both
      the row activity and the value of the associated dual.

      Which to print? Violated constraints will always be flagged to print.
      Otherwise, if m < 50 or all rows are requested, print all rows.  Otherwise,
      print tight constraints (non-zero dual).

      All of this is filtered through printMask, if specified.
    */
    double primalTolerance ;
    osi->getDblParam(OsiPrimalTolerance, primalTolerance) ;

    int iRow ;
    if (ctlBlk->printMode_ >= 3) {
        const double *dualRowSolution = osi->getRowPrice() ;
        const double *primalRowSolution = osi->getRowActivity() ;
        const double *rowLower = osi->getRowLower() ;
        const double *rowUpper = osi->getRowUpper() ;

        fprintf(fp, "\n   %7s %-*s%15s%15s\n\n",
                "Index", longestName, "Row", "Activity", "Dual") ;

        for (iRow = 0 ; iRow < m ; iRow++) {
            bool violated = false ;
            bool print = false ;
            if (primalRowSolution[iRow] > rowUpper[iRow] + primalTolerance ||
                    primalRowSolution[iRow] < rowLower[iRow] - primalTolerance) {
                violated = true ;
                print = true ;
            } else {
                if (m < 50 || ctlBlk->printMode_ >= 4) {
                    print = true ;
                } else if (fabs(dualRowSolution[iRow]) > 1.0e-8) {
                    print = true ;
                }
            }
            const char *name = osi->getRowName(iRow).c_str() ;
            if (doMask && !maskMatches(maskStarts, masks, name)) {
                print = false ;
            }
            if (print) {
                if (violated) {
                    fprintf(fp, "** ") ;
                } else {
                    fprintf(fp, "%3s", " ") ;
                }
                fprintf(fp, "%7d %-*s%15.8g%15.8g\n", iRow, longestName, name,
                        primalRowSolution[iRow], dualRowSolution[iRow]) ;
            }
        }
        fprintf(fp, "\n") ;
    }
    /*
      Now do the columns. This first block handles all modes except 2 (special).
      Out-of-bounds variables are flagged with `**'. If there are less than 50
      variables, all are printed. All of this is filtered through `integer only'
      and can be further filtered using printMask.
    */
    if (ctlBlk->printMode_ != 2) {
        const double *columnLower = osi->getColLower() ;
        const double *columnUpper = osi->getColUpper() ;
        const double *dualColSolution = osi->getReducedCost() ;

        fprintf(fp, "\n   %7s %-*s%15s%15s\n\n",
                "Index", longestName, "Column", "Value", "Reduced Cost") ;

        for (iColumn = 0 ; iColumn < n ; iColumn++) {
            bool violated = false ;
            bool print = false ;
            if (primalColSolution[iColumn] > columnUpper[iColumn] + primalTolerance ||
                    primalColSolution[iColumn] < columnLower[iColumn] - primalTolerance) {
                violated = true ;
                print = true ;
            } else {
                if (n < 50 || ctlBlk->printMode_ == 4) {
                    print = true ;
                } else if (fabs(primalColSolution[iColumn]) > 1.0e-8) {
                    if (ctlBlk->printMode_ == 1) {
                        print = osi->isInteger(iColumn) ;
                    } else {
                        print = true ;
                    }
                }
            }
            const char *name = osi->getColName(iColumn).c_str() ;
            if (doMask && !maskMatches(maskStarts, masks, name)) {
                print = false ;
            }
            if (print) {
                if (violated) {
                    fprintf(fp, "** ") ;
                } else {
                    fprintf(fp, "%3s", " ") ;
                }
                fprintf(fp, "%7d %-*s%15.8g%15.8g\n", iColumn, longestName, name,
                        primalColSolution[iColumn], dualColSolution[iColumn]) ;
            }
        }
    }
    /*
      Close out the file, but don't close stdout. Delete any masks.
    */
    if (fp != stdout) {
        fclose(fp) ;
    }

    if (masks) {
        delete [] maskStarts ;
        for (int i = 0 ; i < maxMasks ; i++)
            delete [] masks[i] ;
        delete [] masks ;
    }

    return (0) ;
}


/*
  Routine to do initial verification of a print mask. We don't generate the
  full set of print masks here, but we do some verification checks to make sure
  it's valid.
*/

int CbcGenParamUtils::doPrintMaskParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Now do a bit of verification of the mask. It should be non-null and, if
      quoted, the quotes should be matched. Aribtrarily put the absolute maximum
      length at 50 characters. If we have a model loaded, that'll be tightened to
      the length of the longest name.
    */
    std::string maskProto = param->strVal() ;
    int maskLen = maskProto.length() ;
    if (maskLen <= 0 || maskLen > 50) {
        std::cerr
            << "Mask |" << maskProto
            << "| is " << maskLen << " characters; should be between "
            << 0 << " and " << 50 << "." << std::endl ;
        return (retval) ;
    }
    /*
      Remove surrounding matched quotes if present. Abort if unmatched.
    */
    if (maskProto[0] == '"' || maskProto[0] == '\'') {
        char quoteChar = maskProto[0] ;
        if (maskProto[maskLen-1] != quoteChar) {
            std::cerr
                << "Mismatched quotes around mask |" << maskProto
                << "|." << std::endl ;
            return (retval) ;
        } else {
            maskProto = maskProto.substr(1, maskLen - 2) ;
        }
    }
    /*
      Mask should not be longer than longest name. Of course, if we don't have a
      model, we can't do this check.
    */
    if (ctlBlk->goodModel_) {
        CbcModel *model = ctlBlk->model_ ;
        assert (model != 0) ;
        OsiSolverInterface *osi = model->solver() ;
        assert (osi != 0) ;
        int longestName = 0 ;
        int n = osi->getNumCols() ;
        for (int j = 0 ; j < n ; j++) {
            int len = osi->getColName(j).length() ;
            longestName = CoinMax(longestName, len) ;
        }
        int m = osi->getNumRows() ;
        for (int i = 0 ; i < m ; i++) {
            int len = osi->getRowName(i).length() ;
            longestName = CoinMax(longestName, len) ;
        }
        if (maskLen > longestName) {
            std::cerr
                << "Mask |" << maskProto << "| has " << maskLen << " chars; this"
                << " is longer than the longest name (" << longestName
                << " chars)." << std::endl ;
            return (retval) ;
        }
    }

    ctlBlk->printMask_ = maskProto ;

    return (0) ;
}

