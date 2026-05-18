// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
#include "CoinUtilsConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinModel.hpp"
#include "CoinMessage.hpp"
#include "CoinSort.hpp"
#include "CoinMpsIO.hpp"
#include "CoinFloatEqual.hpp"
#include "CoinPackedMatrix.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CoinBaseModel::CoinBaseModel()
  : numberRows_(0)
  , numberColumns_(0)
  , optimizationDirection_(1.0)
  , objectiveOffset_(0.0)
  , handler_(NULL)
  , logLevel_(0)
{
  messages_ = CoinMessage();
  handler_ = new CoinMessageHandler();
  problemName_ = "";
  rowBlockName_ = "row_master";
  columnBlockName_ = "column_master";
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CoinBaseModel::CoinBaseModel(const CoinBaseModel &rhs)
  : numberRows_(rhs.numberRows_)
  , numberColumns_(rhs.numberColumns_)
  , optimizationDirection_(rhs.optimizationDirection_)
  , objectiveOffset_(rhs.objectiveOffset_)
  , logLevel_(rhs.logLevel_)
{
  problemName_ = rhs.problemName_;
  rowBlockName_ = rhs.rowBlockName_;
  columnBlockName_ = rhs.columnBlockName_;
  if (rhs.handler_ != NULL)
    handler_ = new CoinMessageHandler(*rhs.handler_) ;
  else
    handler_ = NULL ;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CoinBaseModel::~CoinBaseModel()
{
  delete handler_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CoinBaseModel &
CoinBaseModel::operator=(const CoinBaseModel &rhs)
{
  if (this != &rhs) {
    problemName_ = rhs.problemName_;
    rowBlockName_ = rhs.rowBlockName_;
    columnBlockName_ = rhs.columnBlockName_;
    numberRows_ = rhs.numberRows_;
    numberColumns_ = rhs.numberColumns_;
    optimizationDirection_ = rhs.optimizationDirection_;
    objectiveOffset_ = rhs.objectiveOffset_;
    delete handler_;
    if (rhs.handler_ != NULL)
      handler_ = new CoinMessageHandler(*rhs.handler_) ;
    else
      handler_ = NULL ;
    logLevel_ = rhs.logLevel_;
  }
  return *this;
}
void CoinBaseModel::setLogLevel(int value)
{
  if (value >= 0 && value < 3)
    logLevel_ = value;
}
void CoinBaseModel::setProblemName(const char *name)
{
  if (name)
    problemName_ = name;
  else
    problemName_ = "";
}
// Pass in message handler
void CoinBaseModel::setMessageHandler(CoinMessageHandler *handler)
{
  handler_ = handler;
  if (handler)
    logLevel_ = -1;
  else
    logLevel_ = std::max(0, logLevel_);
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CoinModel::CoinModel()
  : CoinBaseModel()
  , maximumRows_(0)
  , maximumColumns_(0)
  , numberElements_(0)
  , maximumElements_(0)
  , numberQuadraticElements_(0)
  , maximumQuadraticElements_(0)
  , rowLower_(NULL)
  , rowUpper_(NULL)
  , rowType_(NULL)
  , objective_(NULL)
  , columnLower_(NULL)
  , columnUpper_(NULL)
  , integerType_(NULL)
  , columnType_(NULL)
  , start_(NULL)
  , elements_(NULL)
  , packedMatrix_(NULL)
  , quadraticElements_(NULL)
  , sortIndices_(NULL)
  , sortElements_(NULL)
  , sortSize_(0)
  , sizeAssociated_(0)
  , associated_(NULL)
  , numberSOS_(0)
  , startSOS_(NULL)
  , memberSOS_(NULL)
  , typeSOS_(NULL)
  , prioritySOS_(NULL)
  , referenceSOS_(NULL)
  , priority_(NULL)
  , cut_(NULL)
  , moreInfo_(NULL)
  , type_(-1)
  , noNames_(false)
  , links_(0)
{
}
/* Constructor with sizes. */
CoinModel::CoinModel(int firstRows, int firstColumns,
  CoinBigIndex firstElements, bool noNames)
  : CoinBaseModel()
  , maximumRows_(0)
  , maximumColumns_(0)
  , numberElements_(0)
  , maximumElements_(0)
  , numberQuadraticElements_(0)
  , maximumQuadraticElements_(0)
  , rowLower_(NULL)
  , rowUpper_(NULL)
  , rowType_(NULL)
  , objective_(NULL)
  , columnLower_(NULL)
  , columnUpper_(NULL)
  , integerType_(NULL)
  , columnType_(NULL)
  , start_(NULL)
  , elements_(NULL)
  , packedMatrix_(NULL)
  , quadraticElements_(NULL)
  , sortIndices_(NULL)
  , sortElements_(NULL)
  , sortSize_(0)
  , sizeAssociated_(0)
  , associated_(NULL)
  , numberSOS_(0)
  , startSOS_(NULL)
  , memberSOS_(NULL)
  , typeSOS_(NULL)
  , prioritySOS_(NULL)
  , referenceSOS_(NULL)
  , priority_(NULL)
  , cut_(NULL)
  , moreInfo_(NULL)
  , type_(-1)
  , noNames_(noNames)
  , links_(0)
{
  if (!firstRows) {
    if (firstColumns) {
      type_ = 1;
      resize(0, firstColumns, firstElements);
    }
  } else {
    type_ = 0;
    resize(firstRows, 0, firstElements);
    if (firstColumns) {
      // mixed - do linked lists for columns
      //createList(2);
    }
  }
}
/* Read a problem in MPS or GAMS format from the given filename.
 */
CoinModel::CoinModel(const char *fileName, int allowStrings)
  : CoinBaseModel()
  , maximumRows_(0)
  , maximumColumns_(0)
  , numberElements_(0)
  , maximumElements_(0)
  , numberQuadraticElements_(0)
  , maximumQuadraticElements_(0)
  , rowLower_(NULL)
  , rowUpper_(NULL)
  , rowType_(NULL)
  , objective_(NULL)
  , columnLower_(NULL)
  , columnUpper_(NULL)
  , integerType_(NULL)
  , columnType_(NULL)
  , start_(NULL)
  , elements_(NULL)
  , packedMatrix_(NULL)
  , quadraticElements_(NULL)
  , sortIndices_(NULL)
  , sortElements_(NULL)
  , sortSize_(0)
  , sizeAssociated_(0)
  , associated_(NULL)
  , numberSOS_(0)
  , startSOS_(NULL)
  , memberSOS_(NULL)
  , typeSOS_(NULL)
  , prioritySOS_(NULL)
  , referenceSOS_(NULL)
  , priority_(NULL)
  , cut_(NULL)
  , moreInfo_(NULL)
  , type_(-1)
  , noNames_(false)
  , links_(0)
{
  rowBlockName_ = "row_master";
  columnBlockName_ = "column_master";
  int status = 0;
  if (!strcmp(fileName, "-") || !strcmp(fileName, "stdin")) {
    // stdin
  } else {
    std::string name = fileName;
    bool readable = fileCoinReadable(name);
    if (!readable) {
      std::cerr << "Unable to open file "
                << fileName << std::endl;
      status = -1;
    }
  }
  CoinMpsIO m;
  m.setAllowStringElements(allowStrings);
  m.setConvertObjective(true);
  if (!status) {
    try {
      status = m.readMps(fileName, "");
    } catch (CoinError &e) {
      e.print();
      status = -1;
    }
  }
  if (!status) {
    // set problem name
    problemName_ = m.getProblemName();
    objectiveOffset_ = m.objectiveOffset();
    // build model
    int numberRows = m.getNumRows();
    int numberColumns = m.getNumCols();

    // Build by row from scratch
    CoinPackedMatrix matrixByRow = *m.getMatrixByRow();
    const double *element = matrixByRow.getElements();
    const int *column = matrixByRow.getIndices();
    const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
    const int *rowLength = matrixByRow.getVectorLengths();
    const double *rowLower = m.getRowLower();
    const double *rowUpper = m.getRowUpper();
    const double *columnLower = m.getColLower();
    const double *columnUpper = m.getColUpper();
    const double *objective = m.getObjCoefficients();
    int i;
    for (i = 0; i < numberRows; i++) {
      addRow(rowLength[i], column + rowStart[i],
        element + rowStart[i], rowLower[i], rowUpper[i], m.rowName(i));
    }
    int numberIntegers = 0;
    // Now do column part
    for (i = 0; i < numberColumns; i++) {
      setColumnBounds(i, columnLower[i], columnUpper[i]);
      setColumnObjective(i, objective[i]);
      if (m.isInteger(i)) {
        setColumnIsInteger(i, true);
        ;
        numberIntegers++;
      }
    }
    bool quadraticInteger = (numberIntegers != 0) && m.reader()->whichSection() == COIN_QUAD_SECTION;
    // do names
    int iRow;
    for (iRow = 0; iRow < numberRows_; iRow++) {
      const char *name = m.rowName(iRow);
      setRowName(iRow, name);
    }
    bool ifStrings = (m.numberStringElements() != 0);
    int nChanged = 0;
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      // Replace - + or * if strings
      if (!ifStrings && !quadraticInteger) {
        const char *name = m.columnName(iColumn);
        setColumnName(iColumn, name);
      } else {
        assert(strlen(m.columnName(iColumn)) < 100);
        char temp[100];
        strcpy(temp, m.columnName(iColumn));
        int n = CoinStrlenAsInt(temp);
        bool changed = false;
        for (int i = 0; i < n; i++) {
          if (temp[i] == '-') {
            temp[i] = '_';
            changed = true;
          } else if (temp[i] == '+') {
            temp[i] = '$';
            changed = true;
          } else if (temp[i] == '*') {
            temp[i] = '&';
            changed = true;
          }
        }
        if (changed)
          nChanged++;
        setColumnName(iColumn, temp);
      }
    }
    if (nChanged)
      printf("%d column names changed to eliminate - + or *\n", nChanged);
    if (ifStrings) {
      // add in
      int numberElements = m.numberStringElements();
      for (int i = 0; i < numberElements; i++) {
        const char *line = m.stringElement(i);
        int iRow;
        int iColumn;
        sscanf(line, "%d,%d,", &iRow, &iColumn);
        assert(iRow >= 0 && iRow <= numberRows_ + 2);
        assert(iColumn >= 0 && iColumn <= numberColumns_);
        const char *pos = strchr(line, ',');
        assert(pos);
        pos = strchr(pos + 1, ',');
        assert(pos);
        pos++;
        if (iRow < numberRows_ && iColumn < numberColumns_) {
          // element
          setElement(iRow, iColumn, pos);
        } else {
          fprintf(stderr, "code CoinModel strings for rim\n");
          abort();
        }
      }
    }
    // get quadratic part
    if (m.reader()->whichSection() == COIN_QUAD_SECTION) {
      CoinBigIndex *start = NULL;
      int *column = NULL;
      double *element = NULL;
      status = m.readQuadraticMps(NULL, start, column, element, 2);
      if (!status) {
        // If strings allowed 13 then just for Hans convert to constraint
        int objRow = -1;
        if (allowStrings == 13) {
          int objColumn = numberColumns_;
          objRow = numberRows_;
          // leave linear part in objective
          addColumn(0, NULL, NULL, -COIN_DBL_MAX, COIN_DBL_MAX, 1.0, "obj");
          double minusOne = -1.0;
          addRow(1, &objColumn, &minusOne, -COIN_DBL_MAX, 0.0, "objrow");
        }
        if (!ifStrings && !numberIntegers) {
          // no strings - add to quadratic (not done yet)
          for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
            for (CoinBigIndex j = start[iColumn]; j < start[iColumn + 1]; j++) {
              int jColumn = column[j];
              double value = element[j];
              // what about diagonal etc
              if (jColumn == iColumn) {
                printf("diag %d %d %g\n", iColumn, jColumn, value);
                setQuadraticElement(iColumn, jColumn, 0.5 * value);
              } else if (jColumn > iColumn) {
                printf("above diag %d %d %g\n", iColumn, jColumn, value);
              } else if (jColumn < iColumn) {
                printf("below diag %d %d %g\n", iColumn, jColumn, value);
                setQuadraticElement(iColumn, jColumn, value);
              }
            }
          }
        } else {
          // add in as strings
          for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
            char temp[20000];
            temp[0] = '\0';
            int put = 0;
            int n = 0;
            bool ifFirst = true;
            double value = getColumnObjective(iColumn);
            if (value && objRow < 0) {
              sprintf(temp, "%g", value);
              ifFirst = false;
              /* static cast is safe, temp is at most 20000 chars */
              put = CoinStrlenAsInt(temp);
            }
            for (CoinBigIndex j = start[iColumn]; j < start[iColumn + 1]; j++) {
              int jColumn = column[j];
              double value = element[j];
              // what about diagonal etc
              if (jColumn == iColumn) {
                //printf("diag %d %d %g\n",iColumn,jColumn,value);
                value *= 0.5;
              } else if (jColumn > iColumn) {
                //printf("above diag %d %d %g\n",iColumn,jColumn,value);
              } else if (jColumn < iColumn) {
                //printf("below diag %d %d %g\n",iColumn,jColumn,value);
                value = 0.0;
              }
              if (value) {
                n++;
                const char *name = columnName(jColumn);
                if (value == 1.0) {
                  sprintf(temp + put, "%s%s", ifFirst ? "" : "+", name);
                } else {
                  if (ifFirst || value < 0.0)
                    sprintf(temp + put, "%g*%s", value, name);
                  else
                    sprintf(temp + put, "+%g*%s", value, name);
                }
                put += CoinStrlenAsInt(temp + put);
                assert(put < 20000);
                ifFirst = false;
              }
            }
            if (n) {
              if (objRow < 0)
                setObjective(iColumn, temp);
              else
                setElement(objRow, iColumn, temp);
              //printf("el for objective column c%7.7d is %s\n",iColumn,temp);
            }
          }
        }
      }
      delete[] start;
      delete[] column;
      delete[] element;
    }
  }
}
// From arrays
CoinModel::CoinModel(int numberRows, int numberColumns,
  const CoinPackedMatrix *matrix,
  const double *rowLower, const double *rowUpper,
  const double *columnLower, const double *columnUpper,
  const double *objective)
  : CoinBaseModel()
  , maximumRows_(numberRows)
  , maximumColumns_(numberColumns)
  , numberElements_(matrix->getNumElements())
  , maximumElements_(matrix->getNumElements())
  , numberQuadraticElements_(0)
  , maximumQuadraticElements_(0)
  , rowType_(NULL)
  , integerType_(NULL)
  , columnType_(NULL)
  , start_(NULL)
  , elements_(NULL)
  , packedMatrix_(NULL)
  , quadraticElements_(NULL)
  , sortIndices_(NULL)
  , sortElements_(NULL)
  , sortSize_(0)
  , sizeAssociated_(0)
  , associated_(NULL)
  , numberSOS_(0)
  , startSOS_(NULL)
  , memberSOS_(NULL)
  , typeSOS_(NULL)
  , prioritySOS_(NULL)
  , referenceSOS_(NULL)
  , priority_(NULL)
  , cut_(NULL)
  , moreInfo_(NULL)
  , type_(-1)
  , noNames_(false)
  , links_(0)
{
  numberRows_ = numberRows;
  numberColumns_ = numberColumns;
  assert(numberRows_ >= matrix->getNumRows());
  assert(numberColumns_ >= matrix->getNumCols());
  type_ = 3;
  packedMatrix_ = new CoinPackedMatrix(*matrix);
  rowLower_ = CoinCopyOfArray(rowLower, numberRows_);
  rowUpper_ = CoinCopyOfArray(rowUpper, numberRows_);
  objective_ = CoinCopyOfArray(objective, numberColumns_);
  columnLower_ = CoinCopyOfArray(columnLower, numberColumns_);
  columnUpper_ = CoinCopyOfArray(columnUpper, numberColumns_);
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CoinModel::CoinModel(const CoinModel &rhs)
  : CoinBaseModel(rhs)
  , maximumRows_(rhs.maximumRows_)
  , maximumColumns_(rhs.maximumColumns_)
  , numberElements_(rhs.numberElements_)
  , maximumElements_(rhs.maximumElements_)
  , numberQuadraticElements_(rhs.numberQuadraticElements_)
  , maximumQuadraticElements_(rhs.maximumQuadraticElements_)
  , rowName_(rhs.rowName_)
  , columnName_(rhs.columnName_)
  , string_(rhs.string_)
  , hashElements_(rhs.hashElements_)
  , rowList_(rhs.rowList_)
  , columnList_(rhs.columnList_)
  , hashQuadraticElements_(rhs.hashQuadraticElements_)
  , sortSize_(rhs.sortSize_)
  , quadraticRowList_(rhs.quadraticRowList_)
  , quadraticColumnList_(rhs.quadraticColumnList_)
  , sizeAssociated_(rhs.sizeAssociated_)
  , numberSOS_(rhs.numberSOS_)
  , type_(rhs.type_)
  , noNames_(rhs.noNames_)
  , links_(rhs.links_)
{
  rowLower_ = CoinCopyOfArray(rhs.rowLower_, maximumRows_);
  rowUpper_ = CoinCopyOfArray(rhs.rowUpper_, maximumRows_);
  rowType_ = CoinCopyOfArray(rhs.rowType_, maximumRows_);
  objective_ = CoinCopyOfArray(rhs.objective_, maximumColumns_);
  columnLower_ = CoinCopyOfArray(rhs.columnLower_, maximumColumns_);
  columnUpper_ = CoinCopyOfArray(rhs.columnUpper_, maximumColumns_);
  integerType_ = CoinCopyOfArray(rhs.integerType_, maximumColumns_);
  columnType_ = CoinCopyOfArray(rhs.columnType_, maximumColumns_);
  sortIndices_ = CoinCopyOfArray(rhs.sortIndices_, sortSize_);
  sortElements_ = CoinCopyOfArray(rhs.sortElements_, sortSize_);
  associated_ = CoinCopyOfArray(rhs.associated_, sizeAssociated_);
  priority_ = CoinCopyOfArray(rhs.priority_, maximumColumns_);
  cut_ = CoinCopyOfArray(rhs.cut_, maximumRows_);
  moreInfo_ = rhs.moreInfo_;
  if (rhs.packedMatrix_)
    packedMatrix_ = new CoinPackedMatrix(*rhs.packedMatrix_);
  else
    packedMatrix_ = NULL;
  if (numberSOS_) {
    startSOS_ = CoinCopyOfArray(rhs.startSOS_, numberSOS_ + 1);
    int numberMembers = startSOS_[numberSOS_];
    memberSOS_ = CoinCopyOfArray(rhs.memberSOS_, numberMembers);
    typeSOS_ = CoinCopyOfArray(rhs.typeSOS_, numberSOS_);
    prioritySOS_ = CoinCopyOfArray(rhs.prioritySOS_, numberSOS_);
    referenceSOS_ = CoinCopyOfArray(rhs.referenceSOS_, numberMembers);
  } else {
    startSOS_ = NULL;
    memberSOS_ = NULL;
    typeSOS_ = NULL;
    prioritySOS_ = NULL;
    referenceSOS_ = NULL;
  }
  if (type_ == 0) {
    start_ = CoinCopyOfArray(rhs.start_, maximumRows_ + 1);
  } else if (type_ == 1) {
    start_ = CoinCopyOfArray(rhs.start_, maximumColumns_ + 1);
  } else {
    start_ = NULL;
  }
  elements_ = CoinCopyOfArray(rhs.elements_, maximumElements_);
  quadraticElements_ = CoinCopyOfArray(rhs.quadraticElements_, maximumQuadraticElements_);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CoinModel::~CoinModel()
{
  delete[] rowLower_;
  delete[] rowUpper_;
  delete[] rowType_;
  delete[] objective_;
  delete[] columnLower_;
  delete[] columnUpper_;
  delete[] integerType_;
  delete[] columnType_;
  delete[] start_;
  delete[] elements_;
  delete[] quadraticElements_;
  delete[] sortIndices_;
  delete[] sortElements_;
  delete[] associated_;
  delete[] startSOS_;
  delete[] memberSOS_;
  delete[] typeSOS_;
  delete[] prioritySOS_;
  delete[] referenceSOS_;
  delete[] priority_;
  delete[] cut_;
  delete packedMatrix_;
}
// Clone
CoinBaseModel *
CoinModel::clone() const
{
  return new CoinModel(*this);
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CoinModel &
CoinModel::operator=(const CoinModel &rhs)
{
  if (this != &rhs) {
    CoinBaseModel::operator=(rhs);
    delete[] rowLower_;
    delete[] rowUpper_;
    delete[] rowType_;
    delete[] objective_;
    delete[] columnLower_;
    delete[] columnUpper_;
    delete[] integerType_;
    delete[] columnType_;
    delete[] start_;
    delete[] elements_;
    delete[] quadraticElements_;
    delete[] sortIndices_;
    delete[] sortElements_;
    delete[] associated_;
    delete[] startSOS_;
    delete[] memberSOS_;
    delete[] typeSOS_;
    delete[] prioritySOS_;
    delete[] referenceSOS_;
    delete[] priority_;
    delete[] cut_;
    delete packedMatrix_;
    maximumRows_ = rhs.maximumRows_;
    maximumColumns_ = rhs.maximumColumns_;
    numberElements_ = rhs.numberElements_;
    maximumElements_ = rhs.maximumElements_;
    numberQuadraticElements_ = rhs.numberQuadraticElements_;
    maximumQuadraticElements_ = rhs.maximumQuadraticElements_;
    sortSize_ = rhs.sortSize_;
    rowName_ = rhs.rowName_;
    columnName_ = rhs.columnName_;
    string_ = rhs.string_;
    hashElements_ = rhs.hashElements_;
    hashQuadraticElements_ = rhs.hashQuadraticElements_;
    rowList_ = rhs.rowList_;
    quadraticColumnList_ = rhs.quadraticColumnList_;
    quadraticRowList_ = rhs.quadraticRowList_;
    columnList_ = rhs.columnList_;
    sizeAssociated_ = rhs.sizeAssociated_;
    numberSOS_ = rhs.numberSOS_;
    type_ = rhs.type_;
    noNames_ = rhs.noNames_;
    links_ = rhs.links_;
    rowLower_ = CoinCopyOfArray(rhs.rowLower_, maximumRows_);
    rowUpper_ = CoinCopyOfArray(rhs.rowUpper_, maximumRows_);
    rowType_ = CoinCopyOfArray(rhs.rowType_, maximumRows_);
    objective_ = CoinCopyOfArray(rhs.objective_, maximumColumns_);
    columnLower_ = CoinCopyOfArray(rhs.columnLower_, maximumColumns_);
    columnUpper_ = CoinCopyOfArray(rhs.columnUpper_, maximumColumns_);
    integerType_ = CoinCopyOfArray(rhs.integerType_, maximumColumns_);
    columnType_ = CoinCopyOfArray(rhs.columnType_, maximumColumns_);
    priority_ = CoinCopyOfArray(rhs.priority_, maximumColumns_);
    cut_ = CoinCopyOfArray(rhs.cut_, maximumRows_);
    moreInfo_ = rhs.moreInfo_;
    if (rhs.packedMatrix_)
      packedMatrix_ = new CoinPackedMatrix(*rhs.packedMatrix_);
    else
      packedMatrix_ = NULL;
    if (numberSOS_) {
      startSOS_ = CoinCopyOfArray(rhs.startSOS_, numberSOS_ + 1);
      int numberMembers = startSOS_[numberSOS_];
      memberSOS_ = CoinCopyOfArray(rhs.memberSOS_, numberMembers);
      typeSOS_ = CoinCopyOfArray(rhs.typeSOS_, numberSOS_);
      prioritySOS_ = CoinCopyOfArray(rhs.prioritySOS_, numberSOS_);
      referenceSOS_ = CoinCopyOfArray(rhs.referenceSOS_, numberMembers);
    } else {
      startSOS_ = NULL;
      memberSOS_ = NULL;
      typeSOS_ = NULL;
      prioritySOS_ = NULL;
      referenceSOS_ = NULL;
    }
    if (type_ == 0) {
      start_ = CoinCopyOfArray(rhs.start_, maximumRows_ + 1);
    } else if (type_ == 1) {
      start_ = CoinCopyOfArray(rhs.start_, maximumColumns_ + 1);
    } else {
      start_ = NULL;
    }
    elements_ = CoinCopyOfArray(rhs.elements_, maximumElements_);
    quadraticElements_ = CoinCopyOfArray(rhs.quadraticElements_, maximumQuadraticElements_);
    sortIndices_ = CoinCopyOfArray(rhs.sortIndices_, sortSize_);
    sortElements_ = CoinCopyOfArray(rhs.sortElements_, sortSize_);
    associated_ = CoinCopyOfArray(rhs.associated_, sizeAssociated_);
  }
  return *this;
}
/* add a row -  numberInRow may be zero */
void CoinModel::addRow(int numberInRow, const int *columns,
  const double *elements, double rowLower,
  double rowUpper, const char *name)
{
  if (type_ == -1) {
    // initial
    type_ = 0;
    resize(100, 0, 1000);
  } else if (type_ == 1) {
    // mixed - do linked lists for rows
    createList(1);
  } else if (type_ == 3) {
    badType();
  }
  int newColumn = -1;
  if (numberInRow > 0) {
    // Move and sort
    if (numberInRow > sortSize_) {
      delete[] sortIndices_;
      delete[] sortElements_;
      sortSize_ = numberInRow + 100;
      sortIndices_ = new int[sortSize_];
      sortElements_ = new double[sortSize_];
    }
    bool sorted = true;
    int last = -1;
    int i;
    for (i = 0; i < numberInRow; i++) {
      int k = columns[i];
      if (k <= last)
        sorted = false;
      last = k;
      sortIndices_[i] = k;
      sortElements_[i] = elements[i];
    }
    if (!sorted) {
      CoinSort_2(sortIndices_, sortIndices_ + numberInRow, sortElements_);
    }
    // check for duplicates etc
    if (sortIndices_[0] < 0) {
      printf("bad index %d\n", sortIndices_[0]);
      // clean up
      abort();
    }
    last = -1;
    bool duplicate = false;
    for (i = 0; i < numberInRow; i++) {
      int k = sortIndices_[i];
      if (k == last)
        duplicate = true;
      last = k;
    }
    if (duplicate) {
      printf("duplicates - what do we want\n");
      abort();
    }
    newColumn = std::max(newColumn, last);
  }
  int newRow = 0;
  CoinBigIndex newElement = 0;
  if (numberElements_ + numberInRow > maximumElements_) {
    newElement = (3 * (numberElements_ + numberInRow) / 2) + 1000;
    if (numberRows_ * 10 > maximumRows_ * 9)
      newRow = (maximumRows_ * 3) / 2 + 100;
  }
  if (numberRows_ == maximumRows_)
    newRow = (maximumRows_ * 3) / 2 + 100;
  if (newRow || newColumn >= maximumColumns_ || newElement) {
    if (newColumn < maximumColumns_) {
      // columns okay
      resize(newRow, 0, newElement);
    } else {
      // newColumn will be new numberColumns_
      resize(newRow, (3 * newColumn) / 2 + 100, newElement);
    }
  }
  // If rows extended - take care of that
  fillRows(numberRows_, false, true);
  // Do name
  if (name) {
    rowName_.addHash(numberRows_, name);
  } else if (!noNames_) {
    char name[9];
    Coin8CharacterName('r',numberRows_,name);
    rowName_.addHash(numberRows_, name);
  }
  rowLower_[numberRows_] = rowLower;
  rowUpper_[numberRows_] = rowUpper;
  // If columns extended - take care of that
  fillColumns(newColumn, false);
  if (type_ == 0) {
    // can do simply
    CoinBigIndex put = start_[numberRows_];
    assert(put == numberElements_);
    bool doHash = hashElements_.numberItems() != 0;
    for (int i = 0; i < numberInRow; i++) {
      setRowAndStringInTriple(elements_[put], numberRows_, false);
      //elements_[put].row=static_cast<unsigned int>(numberRows_);
      //elements_[put].string=0;
      elements_[put].column = sortIndices_[i];
      elements_[put].value = sortElements_[i];
      if (doHash)
        hashElements_.addHash(put, numberRows_, sortIndices_[i], elements_);
      put++;
    }
    start_[numberRows_ + 1] = put;
    numberElements_ += numberInRow;
  } else {
    if (numberInRow) {
      // must update at least one link
      assert(links_);
      if (links_ == 1 || links_ == 3) {
        CoinBigIndex first = rowList_.addEasy(numberRows_, numberInRow, sortIndices_, sortElements_, elements_,
          hashElements_);
        if (links_ == 3)
          columnList_.addHard(first, elements_, rowList_.firstFree(), rowList_.lastFree(),
            rowList_.next());
        numberElements_ = std::max(numberElements_, rowList_.numberElements());
        if (links_ == 3)
          assert(columnList_.numberElements() == rowList_.numberElements());
      } else if (links_ == 2) {
        columnList_.addHard(numberRows_, numberInRow, sortIndices_, sortElements_, elements_,
          hashElements_);
        numberElements_ = std::max(numberElements_, columnList_.numberElements());
      }
    }
    numberElements_ = std::max(numberElements_, hashElements_.numberItems());
  }
  numberRows_++;
}
// add a column - numberInColumn may be zero */
void CoinModel::addColumn(int numberInColumn, const int *rows,
  const double *elements,
  double columnLower,
  double columnUpper, double objectiveValue,
  const char *name, bool isInteger)
{
  if (type_ == -1) {
    // initial
    type_ = 1;
    resize(0, 100, 1000);
  } else if (type_ == 0) {
    // mixed - do linked lists for columns
    createList(2);
  } else if (type_ == 3) {
    badType();
  }
  int newRow = -1;
  if (numberInColumn > 0) {
    // Move and sort
    if (numberInColumn > sortSize_) {
      delete[] sortIndices_;
      delete[] sortElements_;
      sortSize_ = numberInColumn + 100;
      sortIndices_ = new int[sortSize_];
      sortElements_ = new double[sortSize_];
    }
    bool sorted = true;
    int last = -1;
    int i;
    for (i = 0; i < numberInColumn; i++) {
      int k = rows[i];
      if (k <= last)
        sorted = false;
      last = k;
      sortIndices_[i] = k;
      sortElements_[i] = elements[i];
    }
    if (!sorted) {
      CoinSort_2(sortIndices_, sortIndices_ + numberInColumn, sortElements_);
    }
    // check for duplicates etc
    if (sortIndices_[0] < 0) {
      printf("bad index %d\n", sortIndices_[0]);
      // clean up
      abort();
    }
    last = -1;
    bool duplicate = false;
    for (i = 0; i < numberInColumn; i++) {
      int k = sortIndices_[i];
      if (k == last)
        duplicate = true;
      last = k;
    }
    if (duplicate) {
      printf("duplicates - what do we want\n");
      abort();
    }
    newRow = std::max(newRow, last);
  }
  int newColumn = 0;
  CoinBigIndex newElement = 0;
  if (numberElements_ + numberInColumn > maximumElements_) {
    newElement = (3 * (numberElements_ + numberInColumn) / 2) + 1000;
    if (numberColumns_ * 10 > maximumColumns_ * 9)
      newColumn = (maximumColumns_ * 3) / 2 + 100;
  }
  if (numberColumns_ == maximumColumns_)
    newColumn = (maximumColumns_ * 3) / 2 + 100;
  if (newColumn || newRow >= maximumRows_ || newElement) {
    if (newRow < maximumRows_) {
      // rows okay
      resize(0, newColumn, newElement);
    } else {
      // newRow will be new numberRows_
      resize((3 * newRow) / 2 + 100, newColumn, newElement);
    }
  }
  // If columns extended - take care of that
  fillColumns(numberColumns_, false, true);
  // Do name
  if (name) {
    columnName_.addHash(numberColumns_, name);
  } else if (!noNames_) {
    char name[9];
    Coin8CharacterName('c',numberColumns_,name);
    columnName_.addHash(numberColumns_, name);
  }
  columnLower_[numberColumns_] = columnLower;
  columnUpper_[numberColumns_] = columnUpper;
  objective_[numberColumns_] = objectiveValue;
  if (isInteger)
    integerType_[numberColumns_] = 1;
  else
    integerType_[numberColumns_] = 0;
  // If rows extended - take care of that
  fillRows(newRow, false);
  if (type_ == 1) {
    // can do simply
    CoinBigIndex put = start_[numberColumns_];
    assert(put == numberElements_);
    bool doHash = hashElements_.numberItems() != 0;
    for (int i = 0; i < numberInColumn; i++) {
      elements_[put].column = numberColumns_;
      setRowAndStringInTriple(elements_[put], sortIndices_[i], false);
      //elements_[put].string=0;
      //elements_[put].row=static_cast<unsigned int>(sortIndices_[i]);
      elements_[put].value = sortElements_[i];
      if (doHash)
        hashElements_.addHash(put, sortIndices_[i], numberColumns_, elements_);
      put++;
    }
    start_[numberColumns_ + 1] = put;
    numberElements_ += numberInColumn;
  } else {
    if (numberInColumn) {
      // must update at least one link
      assert(links_);
      if (links_ == 2 || links_ == 3) {
        CoinBigIndex first = columnList_.addEasy(numberColumns_, numberInColumn, sortIndices_, sortElements_, elements_,
          hashElements_);
        if (links_ == 3)
          rowList_.addHard(first, elements_, columnList_.firstFree(), columnList_.lastFree(),
            columnList_.next());
        numberElements_ = std::max(numberElements_, columnList_.numberElements());
        if (links_ == 3)
          assert(columnList_.numberElements() == rowList_.numberElements());
      } else if (links_ == 1) {
        rowList_.addHard(numberColumns_, numberInColumn, sortIndices_, sortElements_, elements_,
          hashElements_);
        numberElements_ = std::max(numberElements_, rowList_.numberElements());
      }
    }
  }
  numberColumns_++;
}
// Sets value for row i and column j
void CoinModel::setElement(int i, int j, double value)
{
  if (type_ == -1) {
    // initial
    type_ = 0;
    resize(100, 100, 1000);
    createList(2);
  } else if (type_ == 3) {
    badType();
  } else if (!links_) {
    if (type_ == 0 || type_ == 2) {
      createList(1);
    } else if (type_ == 1) {
      createList(2);
    }
  }
  if (!hashElements_.maximumItems()) {
    hashElements_.resize(maximumElements_, elements_);
  }
  CoinBigIndex position = hashElements_.hash(i, j, elements_);
  if (position >= 0) {
    elements_[position].value = value;
    setStringInTriple(elements_[position], false);
  } else {
    int newColumn = 0;
    if (j >= maximumColumns_) {
      newColumn = j + 1;
    }
    int newRow = 0;
    if (i >= maximumRows_) {
      newRow = i + 1;
    }
    CoinBigIndex newElement = 0;
    if (numberElements_ == maximumElements_) {
      newElement = (3 * numberElements_ / 2) + 1000;
    }
    if (newRow || newColumn || newElement) {
      if (newColumn)
        newColumn = (3 * newColumn) / 2 + 100;
      if (newRow)
        newRow = (3 * newRow) / 2 + 100;
      resize(newRow, newColumn, newElement);
    }
    // If columns extended - take care of that
    fillColumns(j, false);
    // If rows extended - take care of that
    fillRows(i, false);
    // treat as addRow unless only columnList_ exists
    if ((links_ & 1) != 0) {
      CoinBigIndex first = rowList_.addEasy(i, 1, &j, &value, elements_, hashElements_);
      if (links_ == 3)
        columnList_.addHard(first, elements_, rowList_.firstFree(), rowList_.lastFree(),
          rowList_.next());
      numberElements_ = std::max(numberElements_, rowList_.numberElements());
      if (links_ == 3)
        assert(columnList_.numberElements() == rowList_.numberElements());
    } else if (links_ == 2) {
      columnList_.addHard(i, 1, &j, &value, elements_, hashElements_);
      numberElements_ = std::max(numberElements_, columnList_.numberElements());
    }
    numberRows_ = std::max(numberRows_, i + 1);
    ;
    numberColumns_ = std::max(numberColumns_, j + 1);
    ;
  }
}
// Sets quadratic value for column i and j
void CoinModel::setQuadraticElement(int, int, double)
{
  printf("not written yet\n");
  abort();
  return;
}
// Sets value for row i and column j as string
void CoinModel::setElement(int i, int j, const char *value)
{
  double dummyValue = 1.0;
  if (type_ == -1) {
    // initial
    type_ = 0;
    resize(100, 100, 1000);
    createList(2);
  } else if (type_ == 3) {
    badType();
  } else if (!links_) {
    if (type_ == 0 || type_ == 2) {
      createList(1);
    } else if (type_ == 1) {
      createList(2);
    }
  }
  if (!hashElements_.maximumItems()) {
    // set up number of items
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  CoinBigIndex position = hashElements_.hash(i, j, elements_);
  if (position >= 0) {
    int iValue = addString(value);
    elements_[position].value = iValue;
    setStringInTriple(elements_[position], true);
  } else {
    int newColumn = 0;
    if (j >= maximumColumns_) {
      newColumn = j + 1;
    }
    int newRow = 0;
    if (i >= maximumRows_) {
      newRow = i + 1;
    }
    CoinBigIndex newElement = 0;
    if (numberElements_ == maximumElements_) {
      newElement = (3 * numberElements_ / 2) + 1000;
    }
    if (newRow || newColumn || newElement) {
      if (newColumn)
        newColumn = (3 * newColumn) / 2 + 100;
      if (newRow)
        newRow = (3 * newRow) / 2 + 100;
      resize(newRow, newColumn, newElement);
    }
    // If columns extended - take care of that
    fillColumns(j, false);
    // If rows extended - take care of that
    fillRows(i, false);
    // treat as addRow unless only columnList_ exists
    if ((links_ & 1) != 0) {
      CoinBigIndex first = rowList_.addEasy(i, 1, &j, &dummyValue, elements_, hashElements_);
      if (links_ == 3)
        columnList_.addHard(first, elements_, rowList_.firstFree(), rowList_.lastFree(),
          rowList_.next());
      numberElements_ = std::max(numberElements_, rowList_.numberElements());
      if (links_ == 3)
        assert(columnList_.numberElements() == rowList_.numberElements());
    } else if (links_ == 2) {
      columnList_.addHard(i, 1, &j, &dummyValue, elements_, hashElements_);
      numberElements_ = std::max(numberElements_, columnList_.numberElements());
    }
    numberRows_ = std::max(numberRows_, i + 1);
    ;
    numberColumns_ = std::max(numberColumns_, j + 1);
    ;
    CoinBigIndex position = hashElements_.hash(i, j, elements_);
    assert(position >= 0);
    int iValue = addString(value);
    elements_[position].value = iValue;
    setStringInTriple(elements_[position], true);
  }
}
// Associates a string with a value.  Returns string id (or -1 if does not exist)
int CoinModel::associateElement(const char *stringValue, double value)
{
  int position = string_.hash(stringValue);
  if (position < 0) {
    // not there -add
    position = addString(stringValue);
    assert(position == string_.numberItems() - 1);
  }
  if (sizeAssociated_ <= position) {
    int newSize = (3 * position) / 2 + 100;
    double *temp = new double[newSize];
    CoinMemcpyN(associated_, sizeAssociated_, temp);
    CoinFillN(temp + sizeAssociated_, newSize - sizeAssociated_, unsetValue());
    delete[] associated_;
    associated_ = temp;
    sizeAssociated_ = newSize;
  }
  associated_[position] = value;
  return position;
}
/* Sets rowLower (if row does not exist then
   all rows up to this are defined with default values and no elements)
*/
void CoinModel::setRowLower(int whichRow, double rowLower)
{
  assert(whichRow >= 0);
  // make sure enough room and fill
  fillRows(whichRow, true);
  rowLower_[whichRow] = rowLower;
  rowType_[whichRow] &= ~1;
}
/* Sets rowUpper (if row does not exist then
   all rows up to this are defined with default values and no elements)
*/
void CoinModel::setRowUpper(int whichRow, double rowUpper)
{
  assert(whichRow >= 0);
  // make sure enough room and fill
  fillRows(whichRow, true);
  rowUpper_[whichRow] = rowUpper;
  rowType_[whichRow] &= ~2;
}
/* Sets rowLower and rowUpper (if row does not exist then
   all rows up to this are defined with default values and no elements)
*/
void CoinModel::setRowBounds(int whichRow, double rowLower, double rowUpper)
{
  assert(whichRow >= 0);
  // make sure enough room and fill
  fillRows(whichRow, true);
  rowLower_[whichRow] = rowLower;
  rowUpper_[whichRow] = rowUpper;
  rowType_[whichRow] &= ~3;
}
/* Sets name (if row does not exist then
   all rows up to this are defined with default values and no elements)
*/
void CoinModel::setRowName(int whichRow, const char *rowName)
{
  assert(whichRow >= 0);
  // make sure enough room and fill
  fillRows(whichRow, true);
  assert(!noNames_);
  const char *oldName = rowName_.name(whichRow);
  if (oldName)
    rowName_.deleteHash(whichRow);
  if (rowName)
    rowName_.addHash(whichRow, rowName);
}
/* Sets columnLower (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnLower(int whichColumn, double columnLower)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  columnLower_[whichColumn] = columnLower;
  columnType_[whichColumn] &= ~1;
}
/* Sets columnUpper (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnUpper(int whichColumn, double columnUpper)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  columnUpper_[whichColumn] = columnUpper;
  columnType_[whichColumn] &= ~2;
}
/* Sets columnLower and columnUpper (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnBounds(int whichColumn, double columnLower, double columnUpper)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  columnLower_[whichColumn] = columnLower;
  columnUpper_[whichColumn] = columnUpper;
  columnType_[whichColumn] &= ~3;
}
/* Sets columnObjective (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnObjective(int whichColumn, double columnObjective)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  objective_[whichColumn] = columnObjective;
  columnType_[whichColumn] &= ~4;
}
/* Sets name (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnName(int whichColumn, const char *columnName)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  const char *oldName = columnName_.name(whichColumn);
  assert(!noNames_);
  if (oldName)
    columnName_.deleteHash(whichColumn);
  if (columnName)
    columnName_.addHash(whichColumn, columnName);
}
/* Sets integer (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnIsInteger(int whichColumn, bool columnIsInteger)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  integerType_[whichColumn] = (columnIsInteger) ? 1 : 0;
  columnType_[whichColumn] &= ~8;
}
// Adds one string, returns index
int CoinModel::addString(const char *string)
{
  int position = string_.hash(string);
  if (position < 0) {
    position = string_.numberItems();
    string_.addHash(position, string);
  }
  return position;
}
/* Sets rowLower (if row does not exist then
   all rows up to this are defined with default values and no elements)
*/
void CoinModel::setRowLower(int whichRow, const char *rowLower)
{
  assert(whichRow >= 0);
  // make sure enough room and fill
  fillRows(whichRow, true);
  if (rowLower) {
    int value = addString(rowLower);
    rowLower_[whichRow] = value;
    rowType_[whichRow] |= 1;
  } else {
    rowLower_[whichRow] = -COIN_DBL_MAX;
  }
}
/* Sets rowUpper (if row does not exist then
   all rows up to this are defined with default values and no elements)
*/
void CoinModel::setRowUpper(int whichRow, const char *rowUpper)
{
  assert(whichRow >= 0);
  // make sure enough room and fill
  fillRows(whichRow, true);
  if (rowUpper) {
    int value = addString(rowUpper);
    rowUpper_[whichRow] = value;
    rowType_[whichRow] |= 2;
  } else {
    rowUpper_[whichRow] = COIN_DBL_MAX;
  }
}
/* Sets columnLower (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnLower(int whichColumn, const char *columnLower)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  if (columnLower) {
    int value = addString(columnLower);
    columnLower_[whichColumn] = value;
    columnType_[whichColumn] |= 1;
  } else {
    columnLower_[whichColumn] = 0.0;
  }
}
/* Sets columnUpper (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnUpper(int whichColumn, const char *columnUpper)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  if (columnUpper) {
    int value = addString(columnUpper);
    columnUpper_[whichColumn] = value;
    columnType_[whichColumn] |= 2;
  } else {
    columnUpper_[whichColumn] = COIN_DBL_MAX;
  }
}
/* Sets columnObjective (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnObjective(int whichColumn, const char *columnObjective)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  if (columnObjective) {
    int value = addString(columnObjective);
    objective_[whichColumn] = value;
    columnType_[whichColumn] |= 4;
  } else {
    objective_[whichColumn] = 0.0;
  }
}
/* Sets integer (if column does not exist then
   all columns up to this are defined with default values and no elements)
*/
void CoinModel::setColumnIsInteger(int whichColumn, const char *columnIsInteger)
{
  assert(whichColumn >= 0);
  // make sure enough room and fill
  fillColumns(whichColumn, true);
  if (columnIsInteger) {
    int value = addString(columnIsInteger);
    integerType_[whichColumn] = value;
    columnType_[whichColumn] |= 8;
  } else {
    integerType_[whichColumn] = 0;
  }
}
//static const char * minusInfinity="-infinity";
//static const char * plusInfinity="+infinity";
//static const char * zero="0.0";
static const char *numeric = "Numeric";
/* Gets rowLower (if row does not exist then -COIN_DBL_MAX)
 */
const char *
CoinModel::getRowLowerAsString(int whichRow) const
{
  assert(whichRow >= 0);
  if (whichRow < numberRows_ && rowLower_) {
    if ((rowType_[whichRow] & 1) != 0) {
      int position = static_cast< int >(rowLower_[whichRow]);
      return string_.name(position);
    } else {
      return numeric;
    }
  } else {
    return numeric;
  }
}
/* Gets rowUpper (if row does not exist then +COIN_DBL_MAX)
 */
const char *
CoinModel::getRowUpperAsString(int whichRow) const
{
  assert(whichRow >= 0);
  if (whichRow < numberRows_ && rowUpper_) {
    if ((rowType_[whichRow] & 2) != 0) {
      int position = static_cast< int >(rowUpper_[whichRow]);
      return string_.name(position);
    } else {
      return numeric;
    }
  } else {
    return numeric;
  }
}
/* Gets columnLower (if column does not exist then 0.0)
 */
const char *
CoinModel::getColumnLowerAsString(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && columnLower_) {
    if ((columnType_[whichColumn] & 1) != 0) {
      int position = static_cast< int >(columnLower_[whichColumn]);
      return string_.name(position);
    } else {
      return numeric;
    }
  } else {
    return numeric;
  }
}
/* Gets columnUpper (if column does not exist then COIN_DBL_MAX)
 */
const char *
CoinModel::getColumnUpperAsString(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && columnUpper_) {
    if ((columnType_[whichColumn] & 2) != 0) {
      int position = static_cast< int >(columnUpper_[whichColumn]);
      return string_.name(position);
    } else {
      return numeric;
    }
  } else {
    return numeric;
  }
}
/* Gets columnObjective (if column does not exist then 0.0)
 */
const char *
CoinModel::getColumnObjectiveAsString(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && objective_) {
    if ((columnType_[whichColumn] & 4) != 0) {
      int position = static_cast< int >(objective_[whichColumn]);
      return string_.name(position);
    } else {
      return numeric;
    }
  } else {
    return numeric;
  }
}
/* Gets if integer (if column does not exist then false)
 */
const char *
CoinModel::getColumnIsIntegerAsString(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && integerType_) {
    if ((columnType_[whichColumn] & 8) != 0) {
      int position = integerType_[whichColumn];
      return string_.name(position);
    } else {
      return numeric;
    }
  } else {
    return numeric;
  }
}
/* Deletes all entries in row and bounds.*/
void CoinModel::deleteRow(int whichRow)
{
  assert(whichRow >= 0);
  if (whichRow < numberRows_) {
    if (rowLower_) {
      rowLower_[whichRow] = -COIN_DBL_MAX;
      rowUpper_[whichRow] = COIN_DBL_MAX;
      rowType_[whichRow] = 0;
      if (!noNames_)
        rowName_.deleteHash(whichRow);
    }
    // need lists
    if (type_ == 0) {
      assert(start_);
      assert(!hashElements_.numberItems());
      delete[] start_;
      start_ = NULL;
    }
    if ((links_ & 1) == 0) {
      createList(1);
    }
    assert(links_);
    // row links guaranteed to exist
    rowList_.deleteSame(whichRow, elements_, hashElements_, (links_ != 3));
    // Just need to set first and last and take out
    if (links_ == 3)
      columnList_.updateDeleted(whichRow, elements_, rowList_);
  }
}
/* Deletes all entries in column and bounds.*/
void CoinModel::deleteColumn(int whichColumn)
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_) {
    if (columnLower_) {
      columnLower_[whichColumn] = 0.0;
      columnUpper_[whichColumn] = COIN_DBL_MAX;
      objective_[whichColumn] = 0.0;
      integerType_[whichColumn] = 0;
      columnType_[whichColumn] = 0;
      if (!noNames_)
        columnName_.deleteHash(whichColumn);
    }
    // need lists
    if (type_ == 0) {
      assert(start_);
      assert(!hashElements_.numberItems());
      delete[] start_;
      start_ = NULL;
    } else if (type_ == 3) {
      badType();
    }
    if ((links_ & 2) == 0) {
      createList(2);
    }
    assert(links_);
    // column links guaranteed to exist
    columnList_.deleteSame(whichColumn, elements_, hashElements_, (links_ != 3));
    // Just need to set first and last and take out
    if (links_ == 3)
      rowList_.updateDeleted(whichColumn, elements_, columnList_);
  }
}
// Takes element out of matrix
CoinBigIndex
CoinModel::deleteElement(int row, int column)
{
  CoinBigIndex iPos = position(row, column);
  if (iPos >= 0)
    deleteThisElement(row, column, iPos);
  return iPos;
}
// Takes element out of matrix when position known
void
#ifndef NDEBUG
CoinModel::deleteThisElement(int row, int column, CoinBigIndex position)
#else
CoinModel::deleteThisElement(int, int, CoinBigIndex position)
#endif
{
  assert(row < numberRows_ && column < numberColumns_);
  assert(row == rowInTriple(elements_[position]) && column == static_cast< int >(elements_[position].column));
  if ((links_ & 1) == 0) {
    createList(1);
  }
  assert(links_);
  // row links guaranteed to exist
  rowList_.deleteRowOne(position, elements_, hashElements_);
  // Just need to set first and last and take out
  if (links_ == 3)
    columnList_.updateDeletedOne(position, elements_);
  elements_[position].column = -1;
  elements_[position].value = 0.0;
}
/* Packs down all rows i.e. removes empty rows permanently.  Empty rows
   have no elements and feasible bounds. returns number of rows deleted. */
int CoinModel::packRows()
{
  if (type_ == 3)
    badType();
  int *newRow = new int[numberRows_];
  memset(newRow, 0, numberRows_ * sizeof(int));
  int iRow;
  int n = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    if (rowLower_[iRow] != -COIN_DBL_MAX)
      newRow[iRow]++;
    if (rowUpper_[iRow] != COIN_DBL_MAX)
      newRow[iRow]++;
    if (!noNames_ && rowName_.name(iRow))
      newRow[iRow]++;
  }
  int i;
  for (i = 0; i < numberElements_; i++) {
    if (elements_[i].column >= 0) {
      iRow = rowInTriple(elements_[i]);
      assert(iRow >= 0 && iRow < numberRows_);
      newRow[iRow]++;
    }
  }
  bool doRowNames = (rowName_.numberItems() != 0);
  for (iRow = 0; iRow < numberRows_; iRow++) {
    if (newRow[iRow]) {
      rowLower_[n] = rowLower_[iRow];
      rowUpper_[n] = rowUpper_[iRow];
      rowType_[n] = rowType_[iRow];
      if (doRowNames)
        rowName_.setName(n, rowName_.getName(iRow));
      newRow[iRow] = n++;
    } else {
      newRow[iRow] = -1;
    }
  }
  int numberDeleted = numberRows_ - n;
  if (numberDeleted) {
    numberRows_ = n;
    n = 0;
    for (i = 0; i < numberElements_; i++) {
      if (elements_[i].column >= 0) {
        elements_[n] = elements_[i];
        setRowInTriple(elements_[n], newRow[rowInTriple(elements_[i])]);
        n++;
      }
    }
    numberElements_ = n;
    // now redo
    if (doRowNames) {
      rowName_.setNumberItems(numberRows_);
      rowName_.resize(rowName_.maximumItems(), true);
    }
    if (hashElements_.numberItems()) {
      hashElements_.setNumberItems(numberElements_);
      hashElements_.resize(hashElements_.maximumItems(), elements_, true);
    }
    if (start_) {
      int last = -1;
      if (type_ == 0) {
        for (i = 0; i < numberElements_; i++) {
          int now = rowInTriple(elements_[i]);
          assert(now >= last);
          if (now > last) {
            start_[last + 1] = numberElements_;
            for (int j = last + 1; j < now; j++)
              start_[j + 1] = numberElements_;
            last = now;
          }
        }
        for (int j = last + 1; j < numberRows_; j++)
          start_[j + 1] = numberElements_;
      } else {
        assert(type_ == 1);
        for (i = 0; i < numberElements_; i++) {
          int now = elements_[i].column;
          assert(now >= last);
          if (now > last) {
            start_[last + 1] = numberElements_;
            for (int j = last + 1; j < now; j++)
              start_[j + 1] = numberElements_;
            last = now;
          }
        }
        for (int j = last + 1; j < numberColumns_; j++)
          start_[j + 1] = numberElements_;
      }
    }
    if ((links_ & 1) != 0) {
      rowList_ = CoinModelLinkedList();
      links_ &= ~1;
      createList(1);
    }
    if ((links_ & 2) != 0) {
      columnList_ = CoinModelLinkedList();
      links_ &= ~2;
      createList(2);
    }
  }
  delete[] newRow;
  return numberDeleted;
}
/* Packs down all columns i.e. removes empty columns permanently.  Empty columns
   have no elements and no objective. returns number of columns deleted. */
int CoinModel::packColumns()
{
  if (type_ == 3)
    badType();
  int *newColumn = new int[numberColumns_];
  memset(newColumn, 0, numberColumns_ * sizeof(int));
  int iColumn;
  int n = 0;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if (columnLower_[iColumn] != 0.0)
      newColumn[iColumn]++;
    if (columnUpper_[iColumn] != COIN_DBL_MAX)
      newColumn[iColumn]++;
    if (objective_[iColumn] != 0.0)
      newColumn[iColumn]++;
    if (!noNames_ && columnName_.name(iColumn))
      newColumn[iColumn]++;
  }
  int i;
  for (i = 0; i < numberElements_; i++) {
    if (elements_[i].column >= 0) {
      iColumn = static_cast< int >(elements_[i].column);
      assert(iColumn >= 0 && iColumn < numberColumns_);
      newColumn[iColumn]++;
    }
  }
  bool doColumnNames = (columnName_.numberItems() != 0);
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if (newColumn[iColumn]) {
      columnLower_[n] = columnLower_[iColumn];
      columnUpper_[n] = columnUpper_[iColumn];
      objective_[n] = objective_[iColumn];
      integerType_[n] = integerType_[iColumn];
      columnType_[n] = columnType_[iColumn];
      if (doColumnNames)
        columnName_.setName(n, columnName_.getName(iColumn));
      newColumn[iColumn] = n++;
    } else {
      newColumn[iColumn] = -1;
    }
  }
  int numberDeleted = numberColumns_ - n;
  if (numberDeleted) {
    numberColumns_ = n;
    n = 0;
    for (i = 0; i < numberElements_; i++) {
      if (elements_[i].column >= 0) {
        elements_[n] = elements_[i];
        elements_[n].column = newColumn[elements_[i].column];
        n++;
      }
    }
    numberElements_ = n;
    // now redo
    if (doColumnNames) {
      columnName_.setNumberItems(numberColumns_);
      columnName_.resize(columnName_.maximumItems(), true);
    }
    if (hashElements_.numberItems()) {
      hashElements_.setNumberItems(numberElements_);
      hashElements_.resize(hashElements_.maximumItems(), elements_, true);
    }
    if (start_) {
      int last = -1;
      if (type_ == 0) {
        for (i = 0; i < numberElements_; i++) {
          int now = rowInTriple(elements_[i]);
          assert(now >= last);
          if (now > last) {
            start_[last + 1] = numberElements_;
            for (int j = last + 1; j < now; j++)
              start_[j + 1] = numberElements_;
            last = now;
          }
        }
        for (int j = last + 1; j < numberRows_; j++)
          start_[j + 1] = numberElements_;
      } else {
        assert(type_ == 1);
        for (i = 0; i < numberElements_; i++) {
          int now = elements_[i].column;
          assert(now >= last);
          if (now > last) {
            start_[last + 1] = numberElements_;
            for (int j = last + 1; j < now; j++)
              start_[j + 1] = numberElements_;
            last = now;
          }
        }
        for (int j = last + 1; j < numberColumns_; j++)
          start_[j + 1] = numberElements_;
      }
    }
    if ((links_ & 1) != 0) {
      rowList_ = CoinModelLinkedList();
      links_ &= ~1;
      createList(1);
    }
    if ((links_ & 2) != 0) {
      columnList_ = CoinModelLinkedList();
      links_ &= ~2;
      createList(2);
    }
  }
  delete[] newColumn;
  return numberDeleted;
}
/* Packs down all rows and columns.  i.e. removes empty rows and columns permanently.
   Empty rows have no elements and feasible bounds.
   Empty columns have no elements and no objective.
   returns number of rows+columns deleted. */
int CoinModel::pack()
{
  // For now do slowly (obvious overheads)
  return packRows() + packColumns();
}
// Creates a packed matrix - return sumber of errors
int CoinModel::createPackedMatrix(CoinPackedMatrix &matrix,
  const double *associated)
{
  if (type_ == 3)
    return 0; // badType();
  // Set to say all parts
  type_ = 2;
  resize(numberRows_, numberColumns_, numberElements_);
  // Do counts for CoinPackedMatrix
  int *length = new int[numberColumns_];
  CoinZeroN(length, numberColumns_);
  int i;
  int numberElements = 0;
  for (i = 0; i < numberElements_; i++) {
    int column = elements_[i].column;
    if (column >= 0) {
      length[column]++;
      numberElements++;
    }
  }
  int numberErrors = 0;
  CoinBigIndex *start = new CoinBigIndex[numberColumns_ + 1];
  int *row = new int[numberElements];
  double *element = new double[numberElements];
  start[0] = 0;
  for (i = 0; i < numberColumns_; i++) {
    start[i + 1] = start[i] + length[i];
    length[i] = 0;
  }
  numberElements = 0;
  for (i = 0; i < numberElements_; i++) {
    int column = elements_[i].column;
    if (column >= 0) {
      double value = elements_[i].value;
      if (stringInTriple(elements_[i])) {
        int position = static_cast< int >(value);
        assert(position < sizeAssociated_);
        value = associated[position];
        if (value == unsetValue()) {
          numberErrors++;
          value = 0.0;
        }
      }
      if (value) {
        numberElements++;
        CoinBigIndex put = start[column] + length[column];
        row[put] = rowInTriple(elements_[i]);
        element[put] = value;
        length[column]++;
      }
    }
  }
  for (i = 0; i < numberColumns_; i++) {
    CoinBigIndex put = start[i];
    CoinSort_2(row + put, row + put + length[i], element + put);
  }
  matrix = CoinPackedMatrix(true, numberRows_, numberColumns_, numberElements,
    element, row, start, length, 0.0, 0.0);
  delete[] start;
  delete[] length;
  delete[] row;
  delete[] element;
  return numberErrors;
}
/* Fills in startPositive and startNegative with counts for +-1 matrix.
   If not +-1 then startPositive[0]==-1 otherwise counts and
   startPositive[numberColumns]== size
      - return number of errors
*/
int CoinModel::countPlusMinusOne(CoinBigIndex *startPositive, CoinBigIndex *startNegative,
  const double *associated)
{
  if (type_ == 3)
    badType();
  memset(startPositive, 0, numberColumns_ * sizeof(int));
  memset(startNegative, 0, numberColumns_ * sizeof(int));
  // Set to say all parts
  type_ = 2;
  resize(numberRows_, numberColumns_, numberElements_);
  int numberErrors = 0;
  CoinBigIndex numberElements = 0;
  for (CoinBigIndex i = 0; i < numberElements_; i++) {
    int column = elements_[i].column;
    if (column >= 0) {
      double value = elements_[i].value;
      if (stringInTriple(elements_[i])) {
        int position = static_cast< int >(value);
        assert(position < sizeAssociated_);
        value = associated[position];
        if (value == unsetValue()) {
          numberErrors++;
          value = 0.0;
          startPositive[0] = -1;
          break;
        }
      }
      if (value) {
        numberElements++;
        if (value == 1.0) {
          startPositive[column]++;
        } else if (value == -1.0) {
          startNegative[column]++;
        } else {
          startPositive[0] = -1;
          break;
        }
      }
    }
  }
  if (startPositive[0] >= 0)
    startPositive[numberColumns_] = numberElements;
  return numberErrors;
}
/* Creates +-1 matrix given startPositive and startNegative counts for +-1 matrix.
 */
void CoinModel::createPlusMinusOne(CoinBigIndex *startPositive, CoinBigIndex *startNegative,
  int *indices,
  const double *associated)
{
  if (type_ == 3)
    badType();
  CoinBigIndex size = 0;
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex n = startPositive[iColumn];
    startPositive[iColumn] = size;
    size += n;
    n = startNegative[iColumn];
    startNegative[iColumn] = size;
    size += n;
  }
  startPositive[numberColumns_] = size;
  for (CoinBigIndex i = 0; i < numberElements_; i++) {
    int column = elements_[i].column;
    if (column >= 0) {
      double value = elements_[i].value;
      if (stringInTriple(elements_[i])) {
        int position = static_cast< int >(value);
        assert(position < sizeAssociated_);
        value = associated[position];
      }
      int iRow = rowInTriple(elements_[i]);
      if (value == 1.0) {
        CoinBigIndex position = startPositive[column];
        indices[position] = iRow;
        startPositive[column]++;
      } else if (value == -1.0) {
        CoinBigIndex position = startNegative[column];
        indices[position] = iRow;
        startNegative[column]++;
      }
    }
  }
  // and now redo starts
  for (iColumn = numberColumns_ - 1; iColumn >= 0; iColumn--) {
    startPositive[iColumn + 1] = startNegative[iColumn];
    startNegative[iColumn] = startPositive[iColumn];
  }
  startPositive[0] = 0;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex start = startPositive[iColumn];
    CoinBigIndex end = startNegative[iColumn];
    std::sort(indices + start, indices + end);
    start = startNegative[iColumn];
    end = startPositive[iColumn + 1];
    std::sort(indices + start, indices + end);
  }
}
// Fills in all associated - returning number of errors
int CoinModel::computeAssociated(double *associated)
{
  CoinYacc info;
  info.length = 0;
  int numberErrors = 0;
  for (int i = 0; i < string_.numberItems(); i++) {
    if (string_.name(i) && associated[i] == unsetValue()) {
      associated[i] = getDoubleFromString(info, string_.name(i));
      if (associated[i] == unsetValue())
        numberErrors++;
    }
  }
  return numberErrors;
}
// Creates copies of various arrays - return number of errors
int CoinModel::createArrays(double *&rowLower, double *&rowUpper,
  double *&columnLower, double *&columnUpper,
  double *&objective, int *&integerType,
  double *&associated)
{
  if (sizeAssociated_ < string_.numberItems()) {
    int newSize = string_.numberItems();
    double *temp = new double[newSize];
    CoinMemcpyN(associated_, sizeAssociated_, temp);
    CoinFillN(temp + sizeAssociated_, newSize - sizeAssociated_, unsetValue());
    delete[] associated_;
    associated_ = temp;
    sizeAssociated_ = newSize;
  }
  associated = CoinCopyOfArray(associated_, sizeAssociated_);
  int numberErrors = computeAssociated(associated);
  // Fill in as much as possible
  rowLower = CoinCopyOfArray(rowLower_, numberRows_);
  rowUpper = CoinCopyOfArray(rowUpper_, numberRows_);
  for (int iRow = 0; iRow < numberRows_; iRow++) {
    if ((rowType_[iRow] & 1) != 0) {
      int position = static_cast< int >(rowLower[iRow]);
      assert(position < sizeAssociated_);
      double value = associated[position];
      if (value != unsetValue()) {
        rowLower[iRow] = value;
      }
    }
    if ((rowType_[iRow] & 2) != 0) {
      int position = static_cast< int >(rowUpper[iRow]);
      assert(position < sizeAssociated_);
      double value = associated[position];
      if (value != unsetValue()) {
        rowUpper[iRow] = value;
      }
    }
  }
  columnLower = CoinCopyOfArray(columnLower_, numberColumns_);
  columnUpper = CoinCopyOfArray(columnUpper_, numberColumns_);
  objective = CoinCopyOfArray(objective_, numberColumns_);
  integerType = CoinCopyOfArray(integerType_, numberColumns_);
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if ((columnType_[iColumn] & 1) != 0) {
      int position = static_cast< int >(columnLower[iColumn]);
      assert(position < sizeAssociated_);
      double value = associated[position];
      if (value != unsetValue()) {
        columnLower[iColumn] = value;
      }
    }
    if ((columnType_[iColumn] & 2) != 0) {
      int position = static_cast< int >(columnUpper[iColumn]);
      assert(position < sizeAssociated_);
      double value = associated[position];
      if (value != unsetValue()) {
        columnUpper[iColumn] = value;
      }
    }
    if ((columnType_[iColumn] & 4) != 0) {
      int position = static_cast< int >(objective[iColumn]);
      assert(position < sizeAssociated_);
      double value = associated[position];
      if (value != unsetValue()) {
        objective[iColumn] = value;
      }
    }
    if ((columnType_[iColumn] & 8) != 0) {
      int position = integerType[iColumn];
      assert(position < sizeAssociated_);
      double value = associated[position];
      if (value != unsetValue()) {
        integerType[iColumn] = static_cast< int >(value);
      }
    }
  }
  return numberErrors;
}

/* Write the problem in MPS format to a file with the given filename.
 */
int CoinModel::writeMps(const char *filename, int compression,
  int formatType, int numberAcross, bool keepStrings)
{
  int numberErrors = 0;
  // Set arrays for normal use
  double *rowLower = rowLower_;
  double *rowUpper = rowUpper_;
  double *columnLower = columnLower_;
  double *columnUpper = columnUpper_;
  double *objective = objective_;
  int *integerType = integerType_;
  double *associated = associated_;
  // If strings then do copies
  if (string_.numberItems()) {
    numberErrors = createArrays(rowLower, rowUpper, columnLower, columnUpper,
      objective, integerType, associated);
  }
  CoinPackedMatrix matrix;
  if (type_ != 3) {
    createPackedMatrix(matrix, associated);
  } else {
    matrix = *packedMatrix_;
  }
  char *integrality = new char[numberColumns_];
  bool hasInteger = false;
  for (int i = 0; i < numberColumns_; i++) {
    if (integerType[i]) {
      integrality[i] = 1;
      hasInteger = true;
    } else {
      integrality[i] = 0;
    }
  }

  CoinMpsIO writer;
  writer.setInfinity(COIN_DBL_MAX);
  const char *const *rowNames = NULL;
  if (rowName_.numberItems())
    rowNames = rowName_.names();
  const char *const *columnNames = NULL;
  if (columnName_.numberItems())
    columnNames = columnName_.names();
  writer.setMpsData(matrix, COIN_DBL_MAX,
    columnLower, columnUpper,
    objective, hasInteger ? integrality : 0,
    rowLower, rowUpper,
    columnNames, rowNames);
  delete[] integrality;
  if (rowLower != rowLower_) {
    delete[] rowLower;
    delete[] rowUpper;
    delete[] columnLower;
    delete[] columnUpper;
    delete[] objective;
    delete[] integerType;
    delete[] associated;
    if (numberErrors && logLevel_ > 0 && !keepStrings)
      printf("%d string elements had no values associated with them\n", numberErrors);
  }
  writer.setObjectiveOffset(objectiveOffset_);
  writer.setProblemName(problemName_.c_str());
  if (keepStrings && string_.numberItems()) {
    // load up strings - sorted by column and row
    writer.copyStringElements(this);
  }
  return writer.writeMps(filename, compression, formatType, numberAcross);
}
/* Check two models against each other.  Return nonzero if different.
   Ignore names if that set.
   May modify both models by cleaning up
*/
int CoinModel::differentModel(CoinModel &other, bool ignoreNames)
{
  int numberErrors = 0;
  int numberErrors2 = 0;
  int returnCode = 0;
  if (numberRows_ != other.numberRows_ || numberColumns_ != other.numberColumns_) {
    if (logLevel_ > 0)
      printf("** Mismatch on size, this has %d rows, %d columns - other has %d rows, %d columns\n",
        numberRows_, numberColumns_, other.numberRows_, other.numberColumns_);
    returnCode = 1000;
  }
  // Set arrays for normal use
  double *rowLower = rowLower_;
  double *rowUpper = rowUpper_;
  double *columnLower = columnLower_;
  double *columnUpper = columnUpper_;
  double *objective = objective_;
  int *integerType = integerType_;
  double *associated = associated_;
  // If strings then do copies
  if (string_.numberItems()) {
    numberErrors += createArrays(rowLower, rowUpper, columnLower, columnUpper,
      objective, integerType, associated);
  }
  // Set arrays for normal use
  double *rowLower2 = other.rowLower_;
  double *rowUpper2 = other.rowUpper_;
  double *columnLower2 = other.columnLower_;
  double *columnUpper2 = other.columnUpper_;
  double *objective2 = other.objective_;
  int *integerType2 = other.integerType_;
  double *associated2 = other.associated_;
  // If strings then do copies
  if (other.string_.numberItems()) {
    numberErrors2 += other.createArrays(rowLower2, rowUpper2, columnLower2, columnUpper2,
      objective2, integerType2, associated2);
  }
  CoinPackedMatrix matrix;
  createPackedMatrix(matrix, associated);
  CoinPackedMatrix matrix2;
  other.createPackedMatrix(matrix2, associated2);
  if (numberErrors || numberErrors2)
    if (logLevel_ > 0)
      printf("** Errors when converting strings, %d on this, %d on other\n",
        numberErrors, numberErrors2);
  CoinRelFltEq tolerance;
  if (numberRows_ == other.numberRows_) {
    bool checkNames = !ignoreNames;
    if (!rowName_.numberItems() || !other.rowName_.numberItems())
      checkNames = false;
    int numberDifferentL = 0;
    int numberDifferentU = 0;
    int numberDifferentN = 0;
    for (int i = 0; i < numberRows_; i++) {
      if (!tolerance(rowLower[i], rowLower2[i]))
        numberDifferentL++;
      if (!tolerance(rowUpper[i], rowUpper2[i]))
        numberDifferentU++;
      if (checkNames && rowName_.name(i) && other.rowName_.name(i)) {
        if (strcmp(rowName_.name(i), other.rowName_.name(i)))
          numberDifferentN++;
      }
    }
    int n = numberDifferentL + numberDifferentU + numberDifferentN;
    returnCode += n;
    if (n && logLevel_ > 0)
      printf("Row differences , %d lower, %d upper and %d names\n",
        numberDifferentL, numberDifferentU, numberDifferentN);
  }
  if (numberColumns_ == other.numberColumns_) {
    int numberDifferentL = 0;
    int numberDifferentU = 0;
    int numberDifferentN = 0;
    int numberDifferentO = 0;
    int numberDifferentI = 0;
    bool checkNames = !ignoreNames;
    if (!columnName_.numberItems() || !other.columnName_.numberItems())
      checkNames = false;
    for (int i = 0; i < numberColumns_; i++) {
      if (!tolerance(columnLower[i], columnLower2[i]))
        numberDifferentL++;
      if (!tolerance(columnUpper[i], columnUpper2[i]))
        numberDifferentU++;
      if (!tolerance(objective[i], objective2[i]))
        numberDifferentO++;
      int i1 = (integerType) ? integerType[i] : 0;
      int i2 = (integerType2) ? integerType2[i] : 0;
      if (i1 != i2)
        numberDifferentI++;
      if (checkNames && columnName_.name(i) && other.columnName_.name(i)) {
        if (strcmp(columnName_.name(i), other.columnName_.name(i)))
          numberDifferentN++;
      }
    }
    int n = numberDifferentL + numberDifferentU + numberDifferentN;
    n += numberDifferentO + numberDifferentI;
    returnCode += n;
    if (n && logLevel_ > 0)
      printf("Column differences , %d lower, %d upper, %d objective, %d integer and %d names\n",
        numberDifferentL, numberDifferentU, numberDifferentO,
        numberDifferentI, numberDifferentN);
  }
  if (numberRows_ == other.numberRows_ && numberColumns_ == other.numberColumns_ && numberElements_ == other.numberElements_) {
    if (!matrix.isEquivalent(matrix2, tolerance)) {
      returnCode += 100;
      if (returnCode && logLevel_ > 0)
        printf("Two matrices are not same\n");
    }
  }

  if (rowLower != rowLower_) {
    delete[] rowLower;
    delete[] rowUpper;
    delete[] columnLower;
    delete[] columnUpper;
    delete[] objective;
    delete[] integerType;
    delete[] associated;
  }
  if (rowLower2 != other.rowLower_) {
    delete[] rowLower2;
    delete[] rowUpper2;
    delete[] columnLower2;
    delete[] columnUpper2;
    delete[] objective2;
    delete[] integerType2;
    delete[] associated2;
  }
  return returnCode;
}
// Returns value for row i and column j
double
CoinModel::getElement(int i, int j) const
{
  if (!hashElements_.numberItems()) {
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  CoinBigIndex position = hashElements_.hash(i, j, elements_);
  if (position >= 0) {
    return elements_[position].value;
  } else {
    return 0.0;
  }
}
// Returns value for row rowName and column columnName
double
CoinModel::getElement(const char *rowName, const char *columnName) const
{
  if (!hashElements_.numberItems()) {
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  assert(!noNames_);
  int i = rowName_.hash(rowName);
  int j = columnName_.hash(columnName);
  CoinBigIndex position;
  if (i >= 0 && j >= 0)
    position = hashElements_.hash(i, j, elements_);
  else
    position = -1;
  if (position >= 0) {
    return elements_[position].value;
  } else {
    return 0.0;
  }
}
// Returns quadratic value for columns i and j
double
CoinModel::getQuadraticElement(int, int) const
{
  printf("not written yet\n");
  abort();
  return 0.0;
}
// Returns value for row i and column j as string
const char *
CoinModel::getElementAsString(int i, int j) const
{
  if (!hashElements_.numberItems()) {
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  CoinBigIndex position = hashElements_.hash(i, j, elements_);
  if (position >= 0) {
    if (stringInTriple(elements_[position])) {
      int iString = static_cast< int >(elements_[position].value);
      assert(iString >= 0 && iString < string_.numberItems());
      return string_.name(iString);
    } else {
      return numeric;
    }
  } else {
    return NULL;
  }
}
/* Returns position of element for row i column j.
   Only valid until next modification. 
   -1 if element does not exist */
CoinBigIndex
CoinModel::position(int i, int j) const
{
  if (!hashElements_.numberItems()) {
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_, true);
  }
  return hashElements_.hash(i, j, elements_);
}

/* Returns pointer to element for row i column j.
   Only valid until next modification. 
   NULL if element does not exist */
double *
CoinModel::pointer(int i, int j) const
{
  if (!hashElements_.numberItems()) {
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  CoinBigIndex position = hashElements_.hash(i, j, elements_);
  if (position >= 0) {
    return &(elements_[position].value);
  } else {
    return NULL;
  }
}

/* Returns first element in given row - index is -1 if none.
   Index is given by .index and value by .value
*/
CoinModelLink
CoinModel::firstInRow(int whichRow) const
{
  CoinModelLink link;
  if (whichRow >= 0 && whichRow < numberRows_) {
    link.setOnRow(true);
    if (type_ == 0) {
      assert(start_);
      CoinBigIndex position = start_[whichRow];
      if (position < start_[whichRow + 1]) {
        link.setRow(whichRow);
        link.setPosition(position);
        link.setColumn(elements_[position].column);
        assert(whichRow == rowInTriple(elements_[position]));
        link.setValue(elements_[position].value);
      }
    } else {
      fillList(whichRow, rowList_, 1);
      CoinBigIndex position = rowList_.first(whichRow);
      if (position >= 0) {
        link.setRow(whichRow);
        link.setPosition(position);
        link.setColumn(elements_[position].column);
        assert(whichRow == rowInTriple(elements_[position]));
        link.setValue(elements_[position].value);
      }
    }
  }
  return link;
}
/* Returns last element in given row - index is -1 if none.
   Index is given by .index and value by .value
  */
CoinModelLink
CoinModel::lastInRow(int whichRow) const
{
  CoinModelLink link;
  if (whichRow >= 0 && whichRow < numberRows_) {
    link.setOnRow(true);
    if (type_ == 0) {
      assert(start_);
      CoinBigIndex position = start_[whichRow + 1] - 1;
      if (position >= start_[whichRow]) {
        link.setRow(whichRow);
        link.setPosition(position);
        link.setColumn(elements_[position].column);
        assert(whichRow == rowInTriple(elements_[position]));
        link.setValue(elements_[position].value);
      }
    } else {
      fillList(whichRow, rowList_, 1);
      CoinBigIndex position = rowList_.last(whichRow);
      if (position >= 0) {
        link.setRow(whichRow);
        link.setPosition(position);
        link.setColumn(elements_[position].column);
        assert(whichRow == rowInTriple(elements_[position]));
        link.setValue(elements_[position].value);
      }
    }
  }
  return link;
}
/* Returns first element in given column - index is -1 if none.
   Index is given by .index and value by .value
*/
CoinModelLink
CoinModel::firstInColumn(int whichColumn) const
{
  CoinModelLink link;
  if (whichColumn >= 0 && whichColumn < numberColumns_) {
    link.setOnRow(false);
    if (type_ == 1) {
      assert(start_);
      CoinBigIndex position = start_[whichColumn];
      if (position < start_[whichColumn + 1]) {
        link.setColumn(whichColumn);
        link.setPosition(position);
        link.setRow(rowInTriple(elements_[position]));
        assert(whichColumn == static_cast< int >(elements_[position].column));
        link.setValue(elements_[position].value);
      }
    } else {
      fillList(whichColumn, columnList_, 2);
      if ((links_ & 2) == 0) {
        // Create list
        assert(!columnList_.numberMajor());
        createList(2);
      }
      CoinBigIndex position = columnList_.first(whichColumn);
      if (position >= 0) {
        link.setColumn(whichColumn);
        link.setPosition(position);
        link.setRow(rowInTriple(elements_[position]));
        assert(whichColumn == static_cast< int >(elements_[position].column));
        link.setValue(elements_[position].value);
      }
    }
  }
  return link;
}
/* Returns last element in given column - index is -1 if none.
   Index is given by .index and value by .value
*/
CoinModelLink
CoinModel::lastInColumn(int whichColumn) const
{
  CoinModelLink link;
  if (whichColumn >= 0 && whichColumn < numberColumns_) {
    link.setOnRow(false);
    if (type_ == 1) {
      assert(start_);
      CoinBigIndex position = start_[whichColumn + 1] - 1;
      if (position >= start_[whichColumn]) {
        link.setColumn(whichColumn);
        link.setPosition(position);
        link.setRow(rowInTriple(elements_[position]));
        assert(whichColumn == static_cast< int >(elements_[position].column));
        link.setValue(elements_[position].value);
      }
    } else {
      fillList(whichColumn, columnList_, 2);
      CoinBigIndex position = columnList_.last(whichColumn);
      if (position >= 0) {
        link.setColumn(whichColumn);
        link.setPosition(position);
        link.setRow(rowInTriple(elements_[position]));
        assert(whichColumn == static_cast< int >(elements_[position].column));
        link.setValue(elements_[position].value);
      }
    }
  }
  return link;
}
/* Returns next element in current row or column - index is -1 if none.
   Index is given by .index and value by .value.
   User could also tell because input.next would be NULL
*/
CoinModelLink
CoinModel::next(CoinModelLink &current) const
{
  CoinModelLink link = current;
  CoinBigIndex position = current.position();
  if (position >= 0) {
    if (current.onRow()) {
      // Doing by row
      int whichRow = current.row();
      if (type_ == 0) {
        assert(start_);
        position++;
        if (position < start_[whichRow + 1]) {
          link.setPosition(position);
          link.setColumn(elements_[position].column);
          assert(whichRow == rowInTriple(elements_[position]));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      } else {
        assert((links_ & 1) != 0);
        position = rowList_.next()[position];
        if (position >= 0) {
          link.setPosition(position);
          link.setColumn(elements_[position].column);
          assert(whichRow == rowInTriple(elements_[position]));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      }
    } else {
      // Doing by column
      int whichColumn = current.column();
      if (type_ == 1) {
        assert(start_);
        position++;
        if (position < start_[whichColumn + 1]) {
          link.setPosition(position);
          link.setRow(rowInTriple(elements_[position]));
          assert(whichColumn == static_cast< int >(elements_[position].column));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      } else {
        assert((links_ & 2) != 0);
        position = columnList_.next()[position];
        if (position >= 0) {
          link.setPosition(position);
          link.setRow(rowInTriple(elements_[position]));
          assert(whichColumn == static_cast< int >(elements_[position].column));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      }
    }
  }
  return link;
}
/* Returns previous element in current row or column - index is -1 if none.
   Index is given by .index and value by .value.
   User could also tell because input.previous would be NULL
*/
CoinModelLink
CoinModel::previous(CoinModelLink &current) const
{
  CoinModelLink link = current;
  CoinBigIndex position = current.position();
  if (position >= 0) {
    if (current.onRow()) {
      // Doing by row
      int whichRow = current.row();
      if (type_ == 0) {
        assert(start_);
        position--;
        if (position >= start_[whichRow]) {
          link.setPosition(position);
          link.setColumn(elements_[position].column);
          assert(whichRow == rowInTriple(elements_[position]));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      } else {
        assert((links_ & 1) != 0);
        position = rowList_.previous()[position];
        if (position >= 0) {
          link.setPosition(position);
          link.setColumn(elements_[position].column);
          assert(whichRow == rowInTriple(elements_[position]));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      }
    } else {
      // Doing by column
      int whichColumn = current.column();
      if (type_ == 1) {
        assert(start_);
        position--;
        if (position >= start_[whichColumn]) {
          link.setPosition(position);
          link.setRow(rowInTriple(elements_[position]));
          assert(whichColumn == static_cast< int >(elements_[position].column));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      } else {
        assert((links_ & 2) != 0);
        position = columnList_.previous()[position];
        if (position >= 0) {
          link.setPosition(position);
          link.setRow(rowInTriple(elements_[position]));
          assert(whichColumn == static_cast< int >(elements_[position].column));
          link.setValue(elements_[position].value);
        } else {
          // signal end
          link.setPosition(-1);
          link.setColumn(-1);
          link.setRow(-1);
          link.setValue(0.0);
        }
      }
    }
  }
  return link;
}
/* Returns first element in given quadratic column - index is -1 if none.
   Index is given by .index and value by .value
*/
CoinModelLink
CoinModel::firstInQuadraticColumn(int) const
{
  printf("not written yet\n");
  abort();
  CoinModelLink x;
  return x;
}
/* Returns last element in given quadratic column - index is -1 if none.
   Index is given by .index and value by .value
*/
CoinModelLink
CoinModel::lastInQuadraticColumn(int) const
{
  printf("not written yet\n");
  abort();
  CoinModelLink x;
  return x;
}
/* Gets rowLower (if row does not exist then -COIN_DBL_MAX)
 */
double
CoinModel::getRowLower(int whichRow) const
{
  assert(whichRow >= 0);
  if (whichRow < numberRows_ && rowLower_)
    return rowLower_[whichRow];
  else
    return -COIN_DBL_MAX;
}
/* Gets rowUpper (if row does not exist then +COIN_DBL_MAX)
 */
double
CoinModel::getRowUpper(int whichRow) const
{
  assert(whichRow >= 0);
  if (whichRow < numberRows_ && rowUpper_)
    return rowUpper_[whichRow];
  else
    return COIN_DBL_MAX;
}
/* Gets name (if row does not exist then NULL)
 */
const char *
CoinModel::getRowName(int whichRow) const
{
  assert(whichRow >= 0);
  if (whichRow < rowName_.numberItems())
    return rowName_.name(whichRow);
  else
    return NULL;
}
/* Gets columnLower (if column does not exist then 0.0)
 */
double
CoinModel::getColumnLower(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && columnLower_)
    return columnLower_[whichColumn];
  else
    return 0.0;
}
/* Gets columnUpper (if column does not exist then COIN_DBL_MAX)
 */
double
CoinModel::getColumnUpper(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && columnUpper_)
    return columnUpper_[whichColumn];
  else
    return COIN_DBL_MAX;
}
/* Gets columnObjective (if column does not exist then 0.0)
 */
double
CoinModel::getColumnObjective(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && objective_)
    return objective_[whichColumn];
  else
    return 0.0;
}
/* Gets name (if column does not exist then NULL)
 */
const char *
CoinModel::getColumnName(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < columnName_.numberItems())
    return columnName_.name(whichColumn);
  else
    return NULL;
}
/* Gets if integer (if column does not exist then false)
 */
bool CoinModel::getColumnIsInteger(int whichColumn) const
{
  assert(whichColumn >= 0);
  if (whichColumn < numberColumns_ && integerType_)
    return integerType_[whichColumn] != 0;
  else
    return false;
}
// Row index from row name (-1 if no names or no match)
int CoinModel::row(const char *rowName) const
{
  assert(!noNames_);
  return static_cast< int >(rowName_.hash(rowName));
}
// Column index from column name (-1 if no names or no match)
int CoinModel::column(const char *columnName) const
{
  assert(!noNames_);
  return static_cast< int >(columnName_.hash(columnName));
}
// Resize
void CoinModel::resize(int maximumRows, int maximumColumns, CoinBigIndex maximumElements)
{
  maximumElements = std::max(maximumElements, maximumElements_);
  if (type_ == 0 || type_ == 2) {
    // need to redo row stuff
    maximumRows = std::max(maximumRows, numberRows_);
    if (maximumRows > maximumRows_) {
      bool needFill = rowLower_ == NULL;
      double *tempArray;
      tempArray = new double[maximumRows];
      CoinMemcpyN(rowLower_, numberRows_, tempArray);
#ifdef ZEROFAULT
      memset(tempArray + numberRows_, 0, (maximumRows - numberRows_) * sizeof(double));
#endif
      delete[] rowLower_;
      rowLower_ = tempArray;
      tempArray = new double[maximumRows];
      CoinMemcpyN(rowUpper_, numberRows_, tempArray);
#ifdef ZEROFAULT
      memset(tempArray + numberRows_, 0, (maximumRows - numberRows_) * sizeof(double));
#endif
      delete[] rowUpper_;
      rowUpper_ = tempArray;
      int *tempArray2;
      tempArray2 = new int[maximumRows];
      CoinMemcpyN(rowType_, numberRows_, tempArray2);
#ifdef ZEROFAULT
      memset(tempArray2 + numberRows_, 0, (maximumRows - numberRows_) * sizeof(int));
#endif
      delete[] rowType_;
      rowType_ = tempArray2;
      // resize hash
      if (!noNames_)
        rowName_.resize(maximumRows);
      // If we have links we need to resize
      if ((links_ & 1) != 0) {
        rowList_.resize(maximumRows, maximumElements);
      }
      // If we have start then we need to resize that
      if (type_ == 0) {
        CoinBigIndex *tempArray2;
        tempArray2 = new CoinBigIndex[maximumRows + 1];
#ifdef ZEROFAULT
        memset(tempArray2, 0, (maximumRows + 1) * sizeof(CoinBigIndex));
#endif
        if (start_) {
          CoinMemcpyN(start_, (numberRows_ + 1), tempArray2);
          delete[] start_;
        } else {
          tempArray2[0] = 0;
        }
        start_ = tempArray2;
      }
      maximumRows_ = maximumRows;
      // Fill
      if (needFill) {
        int save = numberRows_ - 1;
        numberRows_ = 0;
        fillRows(save, true);
      }
    }
  } else if (type_ == 3) {
    badType();
  }
  if (type_ == 1 || type_ == 2) {
    // need to redo column stuff
    maximumColumns = std::max(maximumColumns, numberColumns_);
    if (maximumColumns > maximumColumns_) {
      bool needFill = columnLower_ == NULL;
      double *tempArray;
      tempArray = new double[maximumColumns];
      CoinMemcpyN(columnLower_, numberColumns_, tempArray);
#ifdef ZEROFAULT
      memset(tempArray + numberColumns_, 0,
        (maximumColumns - numberColumns_) * sizeof(double));
#endif
      delete[] columnLower_;
      columnLower_ = tempArray;
      tempArray = new double[maximumColumns];
      CoinMemcpyN(columnUpper_, numberColumns_, tempArray);
#ifdef ZEROFAULT
      memset(tempArray + numberColumns_, 0,
        (maximumColumns - numberColumns_) * sizeof(double));
#endif
      delete[] columnUpper_;
      columnUpper_ = tempArray;
      tempArray = new double[maximumColumns];
      CoinMemcpyN(objective_, numberColumns_, tempArray);
#ifdef ZEROFAULT
      memset(tempArray + numberColumns_, 0,
        (maximumColumns - numberColumns_) * sizeof(double));
#endif
      delete[] objective_;
      objective_ = tempArray;
      int *tempArray2;
      tempArray2 = new int[maximumColumns];
      CoinMemcpyN(columnType_, numberColumns_, tempArray2);
#ifdef ZEROFAULT
      memset(tempArray2 + numberColumns_, 0,
        (maximumColumns - numberColumns_) * sizeof(int));
#endif
      delete[] columnType_;
      columnType_ = tempArray2;
      tempArray2 = new int[maximumColumns];
      CoinMemcpyN(integerType_, numberColumns_, tempArray2);
#ifdef ZEROFAULT
      memset(tempArray2 + numberColumns_, 0,
        (maximumColumns - numberColumns_) * sizeof(int));
#endif
      delete[] integerType_;
      integerType_ = tempArray2;
      // resize hash
      if (!noNames_)
        columnName_.resize(maximumColumns);
      // If we have links we need to resize
      if ((links_ & 2) != 0) {
        columnList_.resize(maximumColumns, maximumElements);
      }
      // If we have start then we need to resize that
      if (type_ == 1) {
        CoinBigIndex *tempArray2;
        tempArray2 = new CoinBigIndex[maximumColumns + 1];
#ifdef ZEROFAULT
        memset(tempArray2, 0, (maximumColumns + 1) * sizeof(CoinBigIndex));
#endif
        if (start_) {
          CoinMemcpyN(start_, (numberColumns_ + 1), tempArray2);
          delete[] start_;
        } else {
          tempArray2[0] = 0;
        }
        start_ = tempArray2;
      }
      maximumColumns_ = maximumColumns;
      // Fill
      if (needFill) {
        int save = numberColumns_ - 1;
        numberColumns_ = 0;
        fillColumns(save, true);
      }
    }
  }
  if (type_ == 3)
    badType();
  if (maximumElements > maximumElements_) {
    CoinModelTriple *tempArray = new CoinModelTriple[maximumElements];
    CoinMemcpyN(elements_, numberElements_, tempArray);
#ifdef ZEROFAULT
    memset(tempArray + numberElements_, 0,
      (maximumElements - numberElements_) * sizeof(CoinModelTriple));
#endif
    delete[] elements_;
    elements_ = tempArray;
    if (hashElements_.numberItems())
      hashElements_.resize(maximumElements, elements_);
    maximumElements_ = maximumElements;
    // If we have links we need to resize
    if ((links_ & 1) != 0) {
      rowList_.resize(maximumRows_, maximumElements_);
    }
    if ((links_ & 2) != 0) {
      columnList_.resize(maximumColumns_, maximumElements_);
    }
  }
}
void CoinModel::fillRows(int whichRow, bool forceCreation, bool fromAddRow)
{
  if (forceCreation || fromAddRow) {
    if (type_ == -1) {
      // initial
      type_ = 0;
      resize(std::max(100, whichRow + 1), 0, 1000);
    } else if (type_ == 1) {
      type_ = 2;
    }
    if (!rowLower_) {
      // need to set all
      whichRow = numberRows_ - 1;
      numberRows_ = 0;
      if (type_ != 3)
        resize(std::max(100, whichRow + 1), 0, 0);
      else
        resize(std::max(1, whichRow + 1), 0, 0);
    }
    if (whichRow >= maximumRows_) {
      if (type_ != 3)
        resize(std::max((3 * maximumRows_) / 2, whichRow + 1), 0, 0);
      else
        resize(std::max(1, whichRow + 1), 0, 0);
    }
  }
  if (whichRow >= numberRows_ && rowLower_) {
    // Need to fill
    int i;
    for (i = numberRows_; i <= whichRow; i++) {
      rowLower_[i] = -COIN_DBL_MAX;
      rowUpper_[i] = COIN_DBL_MAX;
      rowType_[i] = 0;
    }
  }
  if (!fromAddRow) {
    numberRows_ = std::max(whichRow + 1, numberRows_);
    // If simple minded then delete start
    if (start_) {
      delete[] start_;
      start_ = NULL;
      assert(!links_);
      // mixed - do linked lists for rows
      createList(1);
    }
  }
}
void CoinModel::fillColumns(int whichColumn, bool forceCreation, bool fromAddColumn)
{
  if (forceCreation || fromAddColumn) {
    if (type_ == -1) {
      // initial
      type_ = 1;
      resize(0, std::max(100, whichColumn + 1), 1000);
    } else if (type_ == 0) {
      type_ = 2;
    }
    if (!objective_) {
      // need to set all
      whichColumn = numberColumns_ - 1;
      numberColumns_ = 0;
      if (type_ != 3)
        resize(0, std::max(100, whichColumn + 1), 0);
      else
        resize(0, std::max(1, whichColumn + 1), 0);
    }
    if (whichColumn >= maximumColumns_) {
      if (type_ != 3)
        resize(0, std::max((3 * maximumColumns_) / 2, whichColumn + 1), 0);
      else
        resize(0, std::max(1, whichColumn + 1), 0);
    }
  }
  if (whichColumn >= numberColumns_ && objective_) {
    // Need to fill
    int i;
    for (i = numberColumns_; i <= whichColumn; i++) {
      columnLower_[i] = 0.0;
      columnUpper_[i] = COIN_DBL_MAX;
      objective_[i] = 0.0;
      integerType_[i] = 0;
      columnType_[i] = 0;
    }
  }
  if (!fromAddColumn) {
    numberColumns_ = std::max(whichColumn + 1, numberColumns_);
    // If simple minded then delete start
    if (start_) {
      delete[] start_;
      start_ = NULL;
      assert(!links_);
      // mixed - do linked lists for columns
      createList(2);
    }
  }
}
// Fill in default linked list information
void CoinModel::fillList(int which, CoinModelLinkedList &list, int type) const
{
  if ((links_ & type) == 0) {
    // Create list
    assert(!list.numberMajor());
    if (type == 1) {
      list.create(maximumRows_, maximumElements_, numberRows_, numberColumns_, 0,
        numberElements_, elements_);
    } else {
      list.create(maximumColumns_, maximumElements_, numberColumns_, numberRows_, 1,
        numberElements_, elements_);
    }
    if (links_ == 1 && type == 2) {
      columnList_.synchronize(rowList_);
    } else if (links_ == 2 && type == 1) {
      rowList_.synchronize(columnList_);
    }
    links_ |= type;
  }
  int number = list.numberMajor();
  if (which >= number) {
    // may still need to extend list or fill it in
    if (which >= list.maximumMajor()) {
      list.resize((which * 3) / 2 + 100, list.maximumElements());
    }
    list.fill(number, which + 1);
  }
}
/* Gets sorted row - user must provide enough space 
   (easiest is allocate number of columns).
   Returns number of elements
*/
int CoinModel::getRow(int whichRow, int *column, double *element)
{
  if (!hashElements_.maximumItems()) {
    // set up number of items
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  assert(whichRow >= 0);
  int n = 0;
  if (whichRow < numberRows_) {
    CoinModelLink triple = firstInRow(whichRow);
    bool sorted = true;
    int last = -1;
    while (triple.column() >= 0) {
      int iColumn = triple.column();
      assert(whichRow == triple.row());
      if (iColumn < last)
        sorted = false;
      last = iColumn;
      if (column)
        column[n] = iColumn;
      if (element)
        element[n] = triple.value();
      n++;
      triple = next(triple);
    }
    if (!sorted) {
      CoinSort_2(column, column + n, element);
    }
  }
  return n;
}
/* Gets sorted column - user must provide enough space 
   (easiest is allocate number of rows).
   Returns number of elements
*/
int CoinModel::getColumn(int whichColumn, int *row, double *element)
{
  if (!hashElements_.maximumItems()) {
    // set up number of items
    hashElements_.setNumberItems(numberElements_);
    hashElements_.resize(maximumElements_, elements_);
  }
  assert(whichColumn >= 0);
  int n = 0;
  if (whichColumn < numberColumns_) {
    CoinModelLink triple = firstInColumn(whichColumn);
    bool sorted = true;
    int last = -1;
    while (triple.column() >= 0) {
      int iRow = triple.row();
      assert(whichColumn == triple.column());
      if (iRow < last)
        sorted = false;
      last = iRow;
      if (row)
        row[n] = iRow;
      if (element)
        element[n] = triple.value();
      n++;
      triple = next(triple);
    }
    if (!sorted) {
      CoinSort_2(row, row + n, element);
    }
  }
  return n;
}
/* Create a linked list and synchronize free 
   type 1 for row 2 for column
   Marked as const as list is mutable */
void CoinModel::createList(int type) const
{
  type_ = 2;
  if (type == 1) {
    assert((links_ & 1) == 0);
    rowList_.create(maximumRows_, maximumElements_,
      numberRows_, numberColumns_, 0,
      numberElements_, elements_);
    if (links_ == 2) {
      // synchronize free list
      rowList_.synchronize(columnList_);
    }
    links_ |= 1;
  } else {
    assert((links_ & 2) == 0);
    columnList_.create(maximumColumns_, maximumElements_,
      numberColumns_, numberRows_, 1,
      numberElements_, elements_);
    if (links_ == 1) {
      // synchronize free list
      columnList_.synchronize(rowList_);
    }
    links_ |= 2;
  }
}
// Checks that links are consistent
void CoinModel::validateLinks() const
{
  if ((links_ & 1)) {
    // validate row links
    rowList_.validateLinks(elements_);
  }
  if ((links_ & 2)) {
    // validate column links
    columnList_.validateLinks(elements_);
  }
}
// returns jColumn (-2 if linear term, -1 if unknown) and coefficient
int CoinModel::decodeBit(char *phrase, char *&nextPhrase, double &coefficient, bool ifFirst) const
{
  char *pos = phrase;
  // may be leading - (or +)
  char *pos2 = pos;
  double value = 1.0;
  if (*pos2 == '-' || *pos2 == '+')
    pos2++;
  // next terminator * or + or -
  while (*pos2) {
    if (*pos2 == '*') {
      break;
    } else if (*pos2 == '-' || *pos2 == '+') {
      if (pos2 == pos || *(pos2 - 1) != 'e')
        break;
    }
    pos2++;
  }
  // if * must be number otherwise must be name
  if (*pos2 == '*') {
    char *pos3 = pos;
    while (pos3 != pos2) {
#ifndef NDEBUG
      char x = *pos3;
#endif
      pos3++;
      assert((x >= '0' && x <= '9') || x == '.' || x == '+' || x == '-' || x == 'e');
    }
    char saved = *pos2;
    *pos2 = '\0';
    value = atof(pos);
    *pos2 = saved;
    // and down to next
    pos2++;
    pos = pos2;
    while (*pos2) {
      if (*pos2 == '-' || *pos2 == '+')
        break;
      pos2++;
    }
  }
  char saved = *pos2;
  *pos2 = '\0';
  // now name
  // might have + or -
  if (*pos == '+') {
    pos++;
  } else if (*pos == '-') {
    pos++;
    assert(value == 1.0);
    value = -value;
  }
  int jColumn = column(pos);
  // must be column unless first when may be linear term
  if (jColumn < 0) {
    if (ifFirst) {
      char *pos3 = pos;
      while (pos3 != pos2) {
#ifndef NDEBUG
        char x = *pos3;
#endif
        pos3++;
        assert((x >= '0' && x <= '9') || x == '.' || x == '+' || x == '-' || x == 'e');
      }
      assert(*pos2 == '\0');
      // keep possible -
      value = value * atof(pos);
      jColumn = -2;
    } else {
      // bad
      *pos2 = saved;
      printf("bad nonlinear term %s\n", phrase);
      // maybe return -1
      abort();
    }
  }
  *pos2 = saved;
  pos = pos2;
  coefficient = value;
  nextPhrase = pos;
  return jColumn;
}
/* Gets correct form for a quadratic row - user to delete
   If row is not quadratic then returns which other variables are involved
   with tiny elements and count of total number of variables which could not
   be put in quadratic form
*/
CoinPackedMatrix *
CoinModel::quadraticRow(int rowNumber, double *linearRow,
  int &numberBad) const
{
  numberBad = 0;
  CoinZeroN(linearRow, numberColumns_);
  int numberElements = 0;
  assert(rowNumber >= -1 && rowNumber < numberRows_);
  if (rowNumber != -1) {
    // not objective
    CoinModelLink triple = firstInRow(rowNumber);
    while (triple.column() >= 0) {
      int iColumn = triple.column();
      const char *expr = getElementAsString(rowNumber, iColumn);
      if (strcmp(expr, "Numeric")) {
        // try and see which columns
        assert(strlen(expr) < 20000);
        char temp[20000];
        strcpy(temp, expr);
        char *pos = temp;
        bool ifFirst = true;
        while (*pos) {
          double value;
          int jColumn = decodeBit(pos, pos, value, ifFirst);
          // must be column unless first when may be linear term
          if (jColumn >= 0) {
            numberElements++;
          } else if (jColumn == -2) {
            linearRow[iColumn] = value;
          } else if (jColumn == -1) {
            // nonlinear term - we will just be marking
            numberElements++;
          } else {
            printf("bad nonlinear term %s\n", temp);
            abort();
          }
          ifFirst = false;
        }
      } else {
        linearRow[iColumn] = getElement(rowNumber, iColumn);
      }
      triple = next(triple);
    }
    if (!numberElements) {
      return NULL;
    } else {
      int *column = new int[numberElements];
      int *column2 = new int[numberElements];
      double *element = new double[numberElements];
      numberElements = 0;
      CoinModelLink triple = firstInRow(rowNumber);
      while (triple.column() >= 0) {
        int iColumn = triple.column();
        const char *expr = getElementAsString(rowNumber, iColumn);
        if (strcmp(expr, "Numeric")) {
          // try and see which columns
          assert(strlen(expr) < 20000);
          char temp[20000];
          strcpy(temp, expr);
          char *pos = temp;
          bool ifFirst = true;
          while (*pos) {
            double value;
            int jColumn = decodeBit(pos, pos, value, ifFirst);
            // must be column unless first when may be linear term
            if (jColumn >= 0) {
              column[numberElements] = iColumn;
              column2[numberElements] = jColumn;
              element[numberElements++] = value;
            } else if (jColumn == -1) {
              // nonlinear term - we will just be marking
              assert(jColumn >= 0);
              column[numberElements] = iColumn;
              column2[numberElements] = jColumn;
              element[numberElements++] = 1.0e-100;
              numberBad++;
            } else if (jColumn != -2) {
              printf("bad nonlinear term %s\n", temp);
              abort();
            }
            ifFirst = false;
          }
        }
        triple = next(triple);
      }
      CoinPackedMatrix *newMatrix = new CoinPackedMatrix(true, column2, column, element, numberElements);
      delete[] column;
      delete[] column2;
      delete[] element;
      return newMatrix;
    }
  } else {
    // objective
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      const char *expr = getColumnObjectiveAsString(iColumn);
      if (strcmp(expr, "Numeric")) {
        // try and see which columns
        assert(strlen(expr) < 20000);
        char temp[20000];
        strcpy(temp, expr);
        char *pos = temp;
        bool ifFirst = true;
        while (*pos) {
          double value;
          int jColumn = decodeBit(pos, pos, value, ifFirst);
          // must be column unless first when may be linear term
          if (jColumn >= 0) {
            numberElements++;
          } else if (jColumn == -2) {
            linearRow[iColumn] = value;
          } else if (jColumn == -1) {
            // nonlinear term - we will just be marking
            numberElements++;
          } else {
            printf("bad nonlinear term %s\n", temp);
            abort();
          }
          ifFirst = false;
        }
      } else {
        linearRow[iColumn] = getElement(rowNumber, iColumn);
      }
    }
    if (!numberElements) {
      return NULL;
    } else {
      int *column = new int[numberElements];
      int *column2 = new int[numberElements];
      double *element = new double[numberElements];
      numberElements = 0;
      for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
        const char *expr = getColumnObjectiveAsString(iColumn);
        if (strcmp(expr, "Numeric")) {
          // try and see which columns
          assert(strlen(expr) < 20000);
          char temp[20000];
          strcpy(temp, expr);
          char *pos = temp;
          bool ifFirst = true;
          while (*pos) {
            double value;
            int jColumn = decodeBit(pos, pos, value, ifFirst);
            // must be column unless first when may be linear term
            if (jColumn >= 0) {
              column[numberElements] = iColumn;
              column2[numberElements] = jColumn;
              element[numberElements++] = value;
            } else if (jColumn == -1) {
              // nonlinear term - we will just be marking
              assert(jColumn >= 0);
              column[numberElements] = iColumn;
              column2[numberElements] = jColumn;
              element[numberElements++] = 1.0e-100;
              numberBad++;
            } else if (jColumn != -2) {
              printf("bad nonlinear term %s\n", temp);
              abort();
            }
            ifFirst = false;
          }
        }
      }
      return new CoinPackedMatrix(true, column2, column, element, numberElements);
    }
  }
}
// Replaces a quadratic row
void CoinModel::replaceQuadraticRow(int rowNumber, const double *linearRow, const CoinPackedMatrix *quadraticPart)
{
  assert(rowNumber >= -1 && rowNumber < numberRows_);
  if (rowNumber >= 0) {
    CoinModelLink triple = firstInRow(rowNumber);
    while (triple.column() >= 0) {
      int iColumn = triple.column();
      deleteElement(rowNumber, iColumn);
      // triple stale - so start over
      triple = firstInRow(rowNumber);
    }
    const double *element = quadraticPart->getElements();
    const int *column = quadraticPart->getIndices();
    const CoinBigIndex *columnStart = quadraticPart->getVectorStarts();
    const int *columnLength = quadraticPart->getVectorLengths();
    int numberLook = quadraticPart->getNumCols();
    int i;
    for (i = 0; i < numberLook; i++) {
      if (!columnLength[i]) {
        // just linear part
        if (linearRow[i])
          setElement(rowNumber, i, linearRow[i]);
      } else {
        char temp[10000];
        int put = 0;
        char temp2[30];
        bool first = true;
        if (linearRow[i]) {
          sprintf(temp, "%g", linearRow[i]);
          first = false;
          /* temp is at most 10000 long, so static_cast is safe */
          put = static_cast< int >(strlen(temp));
        }
        for (CoinBigIndex j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
          int jColumn = column[j];
          double value = element[j];
          if (value < 0.0 || first)
            sprintf(temp2, "%g*c%7.7d", value, jColumn);
          else
            sprintf(temp2, "+%g*c%7.7d", value, jColumn);
          int nextPut = put + static_cast< int >(strlen(temp2));
          assert(nextPut < 10000);
          strcpy(temp + put, temp2);
          put = nextPut;
        }
        setElement(rowNumber, i, temp);
      }
    }
    // rest of linear
    for (; i < numberColumns_; i++) {
      if (linearRow[i])
        setElement(rowNumber, i, linearRow[i]);
    }
  } else {
    // objective
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      setColumnObjective(iColumn, 0.0);
    }
    const double *element = quadraticPart->getElements();
    const int *column = quadraticPart->getIndices();
    const CoinBigIndex *columnStart = quadraticPart->getVectorStarts();
    const int *columnLength = quadraticPart->getVectorLengths();
    int numberLook = quadraticPart->getNumCols();
    int i;
    for (i = 0; i < numberLook; i++) {
      if (!columnLength[i]) {
        // just linear part
        if (linearRow[i])
          setColumnObjective(i, linearRow[i]);
      } else {
        char temp[10000];
        int put = 0;
        char temp2[30];
        bool first = true;
        if (linearRow[i]) {
          sprintf(temp, "%g", linearRow[i]);
          first = false;
          /* temp is at most 10000 long, so static_cast is safe */
          put = static_cast< int >(strlen(temp));
        }
        for (CoinBigIndex j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
          int jColumn = column[j];
          double value = element[j];
          if (value < 0.0 || first)
            sprintf(temp2, "%g*c%7.7d", value, jColumn);
          else
            sprintf(temp2, "+%g*c%7.7d", value, jColumn);
          int nextPut = put + static_cast< int >(strlen(temp2));
          assert(nextPut < 10000);
          strcpy(temp + put, temp2);
          put = nextPut;
        }
        setColumnObjective(i, temp);
      }
    }
    // rest of linear
    for (; i < numberColumns_; i++) {
      if (linearRow[i])
        setColumnObjective(i, linearRow[i]);
    }
  }
}
/* If possible return a model where if all variables marked nonzero are fixed
      the problem will be linear.  At present may only work if quadratic.
      Returns NULL if not possible
*/
CoinModel *
CoinModel::reorder(const char *mark) const
{
  // redo array so 2 high priority nonlinear, 1 nonlinear, 0 linear
  char *highPriority = new char[numberColumns_];
  double *linear = new double[numberColumns_];
  CoinModel *newModel = new CoinModel(*this);
  int iRow;
  for (iRow = -1; iRow < numberRows_; iRow++) {
    int numberBad;
    CoinPackedMatrix *row = quadraticRow(iRow, linear, numberBad);
    assert(!numberBad); // fix later
    if (row) {
      // see if valid
      //const double * element = row->getElements();
      const int *column = row->getIndices();
      const CoinBigIndex *columnStart = row->getVectorStarts();
      const int *columnLength = row->getVectorLengths();
      int numberLook = row->getNumCols();
      for (int i = 0; i < numberLook; i++) {
        if (mark[i])
          highPriority[i] = 2;
        else
          highPriority[i] = 1;
        for (CoinBigIndex j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
          int iColumn = column[j];
          if (mark[iColumn])
            highPriority[iColumn] = 2;
          else
            highPriority[iColumn] = 1;
        }
      }
      delete row;
    }
  }
  for (iRow = -1; iRow < numberRows_; iRow++) {
    int numberBad;
    CoinPackedMatrix *row = quadraticRow(iRow, linear, numberBad);
    if (row) {
      // see if valid
      const double *element = row->getElements();
      const int *columnLow = row->getIndices();
      const CoinBigIndex *columnHigh = row->getVectorStarts();
      const int *columnLength = row->getVectorLengths();
      int numberLook = row->getNumCols();
      int canSwap = 0;
      for (int i = 0; i < numberLook; i++) {
        // this one needs to be available
        int iPriority = highPriority[i];
        for (CoinBigIndex j = columnHigh[i]; j < columnHigh[i] + columnLength[i]; j++) {
          int iColumn = columnLow[j];
          if (highPriority[iColumn] <= 1) {
            assert(highPriority[iColumn] == 1);
            if (iPriority == 1) {
              canSwap = -1; // no good
              break;
            } else {
              canSwap = 1;
            }
          }
        }
      }
      if (canSwap) {
        if (canSwap > 0) {
          // rewrite row
          /* get triples
	     then swap ones needed
	     then create packedmatrix
	     then replace row
	  */
          CoinBigIndex numberElements = columnHigh[numberLook];
          int *columnHigh2 = new int[numberElements];
          int *columnLow2 = new int[numberElements];
          double *element2 = new double[numberElements];
          for (int i = 0; i < numberLook; i++) {
            // this one needs to be available
            int iPriority = highPriority[i];
            if (iPriority == 2) {
              for (CoinBigIndex j = columnHigh[i]; j < columnHigh[i] + columnLength[i]; j++) {
                columnHigh2[j] = i;
                columnLow2[j] = columnLow[j];
                element2[j] = element[j];
              }
            } else {
              for (CoinBigIndex j = columnHigh[i]; j < columnHigh[i] + columnLength[i]; j++) {
                columnLow2[j] = i;
                columnHigh2[j] = columnLow[j];
                element2[j] = element[j];
              }
            }
          }
          delete row;
          row = new CoinPackedMatrix(true, columnHigh2, columnLow2, element2, numberElements);
          delete[] columnHigh2;
          delete[] columnLow2;
          delete[] element2;
          // Now replace row
          newModel->replaceQuadraticRow(iRow, linear, row);
          delete row;
        } else {
          delete row;
          delete newModel;
          newModel = NULL;
          printf("Unable to use priority - row %d\n", iRow);
          break;
        }
      }
    }
  }
  delete[] highPriority;
  delete[] linear;
  return newModel;
}
// Sets cut marker array
void CoinModel::setCutMarker(int size, const int *marker)
{
  delete[] cut_;
  cut_ = new int[maximumRows_];
  CoinZeroN(cut_, maximumRows_);
  CoinMemcpyN(marker, size, cut_);
}
// Sets priority array
void CoinModel::setPriorities(int size, const int *priorities)
{
  delete[] priority_;
  priority_ = new int[maximumColumns_];
  CoinZeroN(priority_, maximumColumns_);
  CoinMemcpyN(priorities, size, priority_);
}
/* Sets columnObjective array
 */
void CoinModel::setObjective(int numberColumns, const double *objective)
{
  fillColumns(numberColumns, true, true);
  for (int i = 0; i < numberColumns; i++) {
    objective_[i] = objective[i];
    columnType_[i] &= ~4;
  }
}
/* Sets columnLower array
 */
void CoinModel::setColumnLower(int numberColumns, const double *columnLower)
{
  fillColumns(numberColumns, true, true);
  for (int i = 0; i < numberColumns; i++) {
    columnLower_[i] = columnLower[i];
    columnType_[i] &= ~1;
  }
}
/* Sets columnUpper array
 */
void CoinModel::setColumnUpper(int numberColumns, const double *columnUpper)
{
  fillColumns(numberColumns, true, true);
  for (int i = 0; i < numberColumns; i++) {
    columnUpper_[i] = columnUpper[i];
    columnType_[i] &= ~2;
  }
}
/* Sets rowLower array
 */
void CoinModel::setRowLower(int numberRows, const double *rowLower)
{
  fillColumns(numberRows, true, true);
  for (int i = 0; i < numberRows; i++) {
    rowLower_[i] = rowLower[i];
    rowType_[i] &= ~1;
  }
}
/* Sets rowUpper array
 */
void CoinModel::setRowUpper(int numberRows, const double *rowUpper)
{
  fillColumns(numberRows, true, true);
  for (int i = 0; i < numberRows; i++) {
    rowUpper_[i] = rowUpper[i];
    rowType_[i] &= ~2;
  }
}
// Pass in CoinPackedMatrix (and switch off element updates)
void CoinModel::passInMatrix(const CoinPackedMatrix &matrix)
{
  type_ = 3;
  packedMatrix_ = new CoinPackedMatrix(matrix);
}
// Convert elements to CoinPackedMatrix (and switch off element updates)
int CoinModel::convertMatrix()
{
  int numberErrors = 0;
  if (type_ != 3) {
    // If strings then do copies
    if (string_.numberItems()) {
      numberErrors = createArrays(rowLower_, rowUpper_,
        columnLower_, columnUpper_,
        objective_, integerType_, associated_);
    }
    CoinPackedMatrix matrix;
    createPackedMatrix(matrix, associated_);
    packedMatrix_ = new CoinPackedMatrix(matrix);
    type_ = 3;
  }
  return numberErrors;
}
// Aborts with message about packedMatrix
void CoinModel::badType() const
{
  fprintf(stderr, "******** operation not allowed when in block mode ****\n");
  abort();
}

//#############################################################################
// Methods to input a problem
//#############################################################################
/** A function to convert from the lb/ub style of constraint
    definition to the sense/rhs/range style */
void convertBoundToSense(const double lower, const double upper,
  char &sense, double &right,
  double &range)
{
  double inf = 1.0e30;
  range = 0.0;
  if (lower > -inf) {
    if (upper < inf) {
      right = upper;
      if (upper == lower) {
        sense = 'E';
      } else {
        sense = 'R';
        range = upper - lower;
      }
    } else {
      sense = 'G';
      right = lower;
    }
  } else {
    if (upper < inf) {
      sense = 'L';
      right = upper;
    } else {
      sense = 'N';
      right = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
/** A function to convert from the sense/rhs/range style of
    constraint definition to the lb/ub style */
void convertSenseToBound(const char sense, const double right,
  const double range,
  double &lower, double &upper)
{
  double inf = COIN_DBL_MAX;
  switch (sense) {
  case 'E':
    lower = upper = right;
    break;
  case 'L':
    lower = -inf;
    upper = right;
    break;
  case 'G':
    lower = right;
    upper = inf;
    break;
  case 'R':
    lower = right - range;
    upper = right;
    break;
  case 'N':
    lower = -inf;
    upper = inf;
    break;
  }
}

void CoinModel::loadBlock(const CoinPackedMatrix &matrix,
  const double *collb, const double *colub,
  const double *obj,
  const double *rowlb, const double *rowub)
{
  passInMatrix(matrix);
  int numberRows = matrix.getNumRows();
  int numberColumns = matrix.getNumCols();
  setObjective(numberColumns, obj);
  setRowLower(numberRows, rowlb);
  setRowUpper(numberRows, rowub);
  setColumnLower(numberColumns, collb);
  setColumnUpper(numberColumns, colub);
}

//-----------------------------------------------------------------------------

void CoinModel::loadBlock(const CoinPackedMatrix &matrix,
  const double *collb, const double *colub,
  const double *obj,
  const char *rowsen, const double *rowrhs,
  const double *rowrng)
{
  // If any of Rhs NULLs then create arrays
  int numrows = matrix.getNumRows();
  const char *rowsenUse = rowsen;
  if (!rowsen) {
    char *rowsen = new char[numrows];
    for (int i = 0; i < numrows; i++)
      rowsen[i] = 'G';
    rowsenUse = rowsen;
  }
  const double *rowrhsUse = rowrhs;
  if (!rowrhs) {
    double *rowrhs = new double[numrows];
    for (int i = 0; i < numrows; i++)
      rowrhs[i] = 0.0;
    rowrhsUse = rowrhs;
  }
  const double *rowrngUse = rowrng;
  if (!rowrng) {
    double *rowrng = new double[numrows];
    for (int i = 0; i < numrows; i++)
      rowrng[i] = 0.0;
    rowrngUse = rowrng;
  }
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];
  for (int i = numrows - 1; i >= 0; --i) {
    convertSenseToBound(rowsenUse[i], rowrhsUse[i], rowrngUse[i], rowlb[i], rowub[i]);
  }
  if (rowsen != rowsenUse)
    delete[] rowsenUse;
  if (rowrhs != rowrhsUse)
    delete[] rowrhsUse;
  if (rowrng != rowrngUse)
    delete[] rowrngUse;
  loadBlock(matrix, collb, colub, obj, rowlb, rowub);
  delete[] rowlb;
  delete[] rowub;
}

//-----------------------------------------------------------------------------

void CoinModel::loadBlock(const int numcols, const int numrows,
  const CoinBigIndex *start, const int *index,
  const double *value,
  const double *collb, const double *colub,
  const double *obj,
  const double *rowlb, const double *rowub)
{
  CoinBigIndex numberElements = start[numcols];
  int *length = new int[numcols];
  for (int i = 0; i < numcols; i++)
    length[i] = static_cast< int >(start[i + 1] - start[i]);
  CoinPackedMatrix matrix(true, numrows, numcols, numberElements, value,
    index, start, length, 0.0, 0.0);
  loadBlock(matrix, collb, colub, obj, rowlb, rowub);
  delete[] length;
}
//-----------------------------------------------------------------------------

void CoinModel::loadBlock(const int numcols, const int numrows,
  const CoinBigIndex *start, const int *index,
  const double *value,
  const double *collb, const double *colub,
  const double *obj,
  const char *rowsen, const double *rowrhs,
  const double *rowrng)
{
  // If any of Rhs NULLs then create arrays
  const char *rowsenUse = rowsen;
  if (!rowsen) {
    char *rowsen = new char[numrows];
    for (int i = 0; i < numrows; i++)
      rowsen[i] = 'G';
    rowsenUse = rowsen;
  }
  const double *rowrhsUse = rowrhs;
  if (!rowrhs) {
    double *rowrhs = new double[numrows];
    for (int i = 0; i < numrows; i++)
      rowrhs[i] = 0.0;
    rowrhsUse = rowrhs;
  }
  const double *rowrngUse = rowrng;
  if (!rowrng) {
    double *rowrng = new double[numrows];
    for (int i = 0; i < numrows; i++)
      rowrng[i] = 0.0;
    rowrngUse = rowrng;
  }
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];
  for (int i = numrows - 1; i >= 0; --i) {
    convertSenseToBound(rowsenUse[i], rowrhsUse[i], rowrngUse[i], rowlb[i], rowub[i]);
  }
  if (rowsen != rowsenUse)
    delete[] rowsenUse;
  if (rowrhs != rowrhsUse)
    delete[] rowrhsUse;
  if (rowrng != rowrngUse)
    delete[] rowrngUse;
  CoinBigIndex numberElements = start[numcols];
  int *length = new int[numcols];
  for (int i = 0; i < numcols; i++)
    length[i] = static_cast< int >(start[i + 1] - start[i]);
  CoinPackedMatrix matrix(true, numrows, numcols, numberElements, value,
    index, start, length, 0.0, 0.0);
  loadBlock(matrix, collb, colub, obj, rowlb, rowub);
  delete[] length;
  delete[] rowlb;
  delete[] rowub;
}
/* Returns which parts of model are set
   1 - matrix
   2 - rhs
   4 - row names
   8 - column bounds and/or objective
   16 - column names
   32 - integer types
*/
int CoinModel::whatIsSet() const
{
  int type = (numberElements_) ? 1 : 0;
  bool defaultValues = true;
  if (rowLower_) {
    for (int i = 0; i < numberRows_; i++) {
      if (rowLower_[i] != -COIN_DBL_MAX) {
        defaultValues = false;
        break;
      }
      if (rowUpper_[i] != COIN_DBL_MAX) {
        defaultValues = false;
        break;
      }
    }
  }
  if (!defaultValues)
    type |= 2;
  if (rowName_.numberItems())
    type |= 4;
  defaultValues = true;
  if (columnLower_) {
    for (int i = 0; i < numberColumns_; i++) {
      if (objective_[i] != 0.0) {
        defaultValues = false;
        break;
      }
      if (columnLower_[i] != 0.0) {
        defaultValues = false;
        break;
      }
      if (columnUpper_[i] != COIN_DBL_MAX) {
        defaultValues = false;
        break;
      }
    }
  }
  if (!defaultValues)
    type |= 8;
  if (columnName_.numberItems())
    type |= 16;
  defaultValues = true;
  if (integerType_) {
    for (int i = 0; i < numberColumns_; i++) {
      if (integerType_[i]) {
        defaultValues = false;
        break;
      }
    }
  }
  if (!defaultValues)
    type |= 32;
  return type;
}
// For decomposition set original row and column indices
void CoinModel::setOriginalIndices(const int *row, const int *column)
{
  if (!rowType_)
    rowType_ = new int[numberRows_];
  memcpy(rowType_, row, numberRows_ * sizeof(int));
  if (!columnType_)
    columnType_ = new int[numberColumns_];
  memcpy(columnType_, column, numberColumns_ * sizeof(int));
}


/* code below was Clp_ampl.cpp in Clp and Cbc_ampl.cpp in Cbc before,
 * but as it implements a CoinModel method, it should be in CoinUtils
 */

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

/* Interface routines for AMPL.
*/

#ifdef COINUTILS_HAS_ASL

#ifdef HAVE_UNISTD_H
#include "unistd.h"
#endif
extern "C" {
#include "getstub.h"
#include "asl_pfgh.h"
}

#include <string>
#include <cassert>
/* so decodePhrase and clpCheck can access */
static ampl_info *saveInfo = NULL;
// Set to 1 if algorithm found
static char algFound[20] = "";
static char *
checkPhrase(Option_Info *oi, keyword *kw, char *v)
{
  if (strlen(v))
    printf("string %s\n", v);
  // Say algorithm found
  strcpy(algFound, kw->desc);
  return v;
}
static char *
checkPhrase2(Option_Info *oi, keyword *kw, char *v)
{
  if (strlen(v))
    printf("string %s\n", v);
  // put out keyword
  saveInfo->arguments = (char **)realloc(saveInfo->arguments, (saveInfo->numberArguments + 1) * sizeof(char *));
  saveInfo->arguments[saveInfo->numberArguments++] = strdup(kw->desc);
  return v;
}
static fint
decodePhrase(char *phrase, ftnlen length)
{
  char *blank = strchr(phrase, ' ');
  if (blank) {
    /* split arguments */
    *blank = '\0';
    saveInfo->arguments = (char **)realloc(saveInfo->arguments, (saveInfo->numberArguments + 2) * sizeof(char *));
    saveInfo->arguments[saveInfo->numberArguments++] = strdup(phrase);
    *blank = ' ';
    phrase = blank + 1; /* move on */
    if (strlen(phrase))
      saveInfo->arguments[saveInfo->numberArguments++] = strdup(phrase);
  } else if (strlen(phrase)) {
    saveInfo->arguments = (char **)realloc(saveInfo->arguments, (saveInfo->numberArguments + 1) * sizeof(char *));
    saveInfo->arguments[saveInfo->numberArguments++] = strdup(phrase);
  }
  return 0;
}
static void
sos_kludge(int nsos, int *sosbeg, double *sosref, int *sosind)
{
  // Adjust sosref if necessary to make monotonic increasing
  int i, j, k;
  // first sort
  for (i = 0; i < nsos; i++) {
    k = sosbeg[i];
    int end = sosbeg[i + 1];
    CoinSort_2(sosref + k, sosref + end, sosind + k);
  }
  double t, t1;
  for (i = j = 0; i++ < nsos;) {
    k = sosbeg[i];
    t = sosref[j];
    while (++j < k) {
      t1 = sosref[j];
      t += 1e-10;
      if (t1 <= t)
        sosref[j] = t1 = t + 1e-10;
      t = t1;
    }
  }
}
static char xxxxxx[20];
#define VP (char *)
static keyword keywds[] = { /* must be sorted */
  { const_cast< char * >("barrier"), checkPhrase, (char *)xxxxxx,
    const_cast< char * >("-barrier") },
  { const_cast< char * >("dual"), checkPhrase, (char *)xxxxxx,
    const_cast< char * >("-dualsimplex") },
  { const_cast< char * >("help"), checkPhrase2, (char *)xxxxxx,
    const_cast< char * >("-?") },
  { const_cast< char * >("initial"), checkPhrase, (char *)xxxxxx,
    const_cast< char * >("-initialsolve") },
  { const_cast< char * >("max"), checkPhrase2, (char *)xxxxxx,
    const_cast< char * >("-maximize") },
  { const_cast< char * >("maximize"), checkPhrase2, (char *)xxxxxx,
    const_cast< char * >("-maximize") },
  { const_cast< char * >("primal"), checkPhrase, (char *)xxxxxx,
    const_cast< char * >("-primalsimplex") },
  { const_cast< char * >("quit"), checkPhrase2, (char *)xxxxxx,
    const_cast< char * >("-quit") },
  { const_cast< char * >("wantsol"), WS_val, NULL,
    const_cast< char * >("write .sol file (without -AMPL)") }
};
static Option_Info Oinfo = {
  NULL,
  NULL,
  NULL,
  keywds,
  nkeywds,
  0,
  0,
  0,
  decodePhrase,
  0,
  0,
  0,
  20130502
};
// strdup used to avoid g++ compiler warning
static SufDecl suftab[] = {
#ifdef JJF_ZERO
  { const_cast< char * >("current"), 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
  { const_cast< char * >("current"), 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
  { const_cast< char * >("direction"), 0, ASL_Sufkind_var },
  { const_cast< char * >("down"), 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
  { const_cast< char * >("down"), 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
  { const_cast< char * >("priority"), 0, ASL_Sufkind_var },
#endif
  { const_cast< char * >("cut"), 0, ASL_Sufkind_con },
  { const_cast< char * >("direction"), 0, ASL_Sufkind_var },
  { const_cast< char * >("downPseudocost"), 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { const_cast< char * >("priority"), 0, ASL_Sufkind_var },
  { const_cast< char * >("ref"), 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { const_cast< char * >("sos"), 0, ASL_Sufkind_var },
  { const_cast< char * >("sos"), 0, ASL_Sufkind_con },
  { const_cast< char * >("sosno"), 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { const_cast< char * >("sosref"), 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { const_cast< char * >("special"), 0, ASL_Sufkind_var },
  { const_cast< char * >("special"), 0, ASL_Sufkind_con },
  /*{ const_cast<char*>("special"), 0, ASL_Sufkind_con },*/
  { const_cast< char * >("sstatus"), 0, ASL_Sufkind_var, 0 },
  { const_cast< char * >("sstatus"), 0, ASL_Sufkind_con, 0 },
  { const_cast< char * >("upPseudocost"), 0, ASL_Sufkind_var | ASL_Sufkind_real }
#ifdef JJF_ZERO
  { const_cast< char * >("unbdd"), 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
  { const_cast< char * >("up"), 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
  { const_cast< char * >("up"), 0, ASL_Sufkind_var | ASL_Sufkind_outonly }
#endif
};
#include "float.h"
#include "limits.h"
static ASL *asl = NULL;
static FILE *nl = NULL;

static void
mip_stuff(void)
{
  int i;
  double *pseudoUp, *pseudoDown;
  int *priority, *direction;
  // To label cuts (there will be other uses for special)
  int *cut;
  // To label special variables - at present 1= must be >= 1 or <= -1
  int *special;
  SufDesc *dpup, *dpdown, *dpri, *ddir, *dcut, *dspecial;

  ddir = suf_get("direction", ASL_Sufkind_var);
  direction = ddir->u.i;
  dpri = suf_get("priority", ASL_Sufkind_var);
  priority = dpri->u.i;
  dspecial = suf_get("special", ASL_Sufkind_con);
  dcut = suf_get("cut", ASL_Sufkind_con);
  cut = dcut->u.i;
  if (!cut) {
    // try special
    dcut = suf_get("special", ASL_Sufkind_con);
    cut = dcut->u.i;
  }
  dspecial = suf_get("special", ASL_Sufkind_var);
  special = dspecial->u.i;
  dpdown = suf_get("downPseudocost", ASL_Sufkind_var);
  pseudoDown = dpdown->u.r;
  dpup = suf_get("upPseudocost", ASL_Sufkind_var);
  pseudoUp = dpup->u.r;
  assert(saveInfo);
  int numberColumns = saveInfo->numberColumns;
  if (direction) {
    int baddir = 0;
    saveInfo->branchDirection = (int *)malloc(numberColumns * sizeof(int));
    for (i = 0; i < numberColumns; i++) {
      int value = direction[i];
      if (value < -1 || value > 1) {
        baddir++;
        value = 0;
      }
      saveInfo->branchDirection[i] = value;
    }
    if (baddir)
      fprintf(stderr,
        "Treating %d .direction values outside [-1, 1] as 0.\n",
        baddir);
  }
  if (priority) {
    int badpri = 0;
    saveInfo->priorities = (int *)malloc(numberColumns * sizeof(int));
    for (i = 0; i < numberColumns; i++) {
      int value = priority[i];
      if (value < 0) {
        badpri++;
        value = 0;
      }
      saveInfo->priorities[i] = value;
    }
    if (badpri)
      fprintf(stderr,
        "Treating %d negative .priority values as 0\n",
        badpri);
  }
  if (special) {
    int badspecial = 0;
    saveInfo->special = (int *)malloc(numberColumns * sizeof(int));
    for (i = 0; i < numberColumns; i++) {
      int value = special[i];
      if (value < 0) {
        badspecial++;
        value = 0;
      }
      saveInfo->special[i] = value;
    }
    if (badspecial)
      fprintf(stderr,
        "Treating %d negative special values as 0\n",
        badspecial);
  }
  int numberRows = saveInfo->numberRows;
  if (cut) {
    int badcut = 0;
    saveInfo->cut = (int *)malloc(numberRows * sizeof(int));
    for (i = 0; i < numberRows; i++) {
      int value = cut[i];
      if (value < 0) {
        badcut++;
        value = 0;
      }
      saveInfo->cut[i] = value;
    }
    if (badcut)
      fprintf(stderr,
        "Treating %d negative cut values as 0\n",
        badcut);
  }
  if (pseudoDown || pseudoUp) {
    int badpseudo = 0;
    if (!pseudoDown || !pseudoUp)
      fprintf(stderr,
        "Only one set of pseudocosts - assumed same\n");
    saveInfo->pseudoDown = (double *)malloc(numberColumns * sizeof(double));
    saveInfo->pseudoUp = (double *)malloc(numberColumns * sizeof(double));
    for (i = 0; i < numberColumns; i++) {
      double valueD = 0.0, valueU = 0.0;
      if (pseudoDown) {
        valueD = pseudoDown[i];
        if (valueD < 0) {
          badpseudo++;
          valueD = 0.0;
        }
      }
      if (pseudoUp) {
        valueU = pseudoUp[i];
        if (valueU < 0) {
          badpseudo++;
          valueU = 0.0;
        }
      }
      if (!valueD)
        valueD = valueU;
      if (!valueU)
        valueU = valueD;
      saveInfo->pseudoDown[i] = valueD;
      saveInfo->pseudoUp[i] = valueU;
    }
    if (badpseudo)
      fprintf(stderr,
        "Treating %d negative pseudoCosts as 0.0\n", badpseudo);
  }
}
static void
stat_map(int *stat, int n, int *map, int mx, const char *what)
{
  int bad, i, i1 = 0, j, j1 = 0;
  static char badfmt[] = "Coin driver: %s[%d] = %d\n";

  for (i = bad = 0; i < n; i++) {
    if ((j = stat[i]) >= 0 && j <= mx)
      stat[i] = map[j];
    else {
      stat[i] = 0;
      i1 = i;
      j1 = j;
      if (!bad++)
        fprintf(stderr, badfmt, what, i, j);
    }
  }
  if (bad > 1) {
    if (bad == 2)
      fprintf(stderr, badfmt, what, i1, j1);
    else
      fprintf(stderr,
        "Coin driver: %d messages about bad %s values suppressed.\n",
        bad - 1, what);
  }
}

int readAmpl(ampl_info *info, int argc, char **argv, void **coinModel, const char* solvername)
{
  char *stub;
  ograd *og;
  int i;
  SufDesc *csd;
  SufDesc *rsd;
  /*bool *basis, *lower;*/
  /*double *LU, *c, lb, objadj, *rshift, *shift, t, ub, *x, *x0, *x1;*/
  char tempBuffer[200];
  sprintf(tempBuffer, "%s_options", solvername);
  char *environment = getenv(tempBuffer);
  double *obj;
  double *columnLower;
  double *columnUpper;
  double *rowLower;
  double *rowUpper;
  char **saveArgv = argv;
  char fileName[1000];
  if (argc > 1)
    strcpy(fileName, argv[1]);
  else
    fileName[0] = '\0';
  int nonLinearType = -1;
  // testosi parameter - if >= 10 then go in through coinModel
  for (i = 1; i < argc; i++) {
    if (!strncmp(argv[i], "testosi", 7)) {
      char *equals = strchr(argv[i], '=');
      if (equals && atoi(equals + 1) >= 10 && atoi(equals + 1) <= 20) {
        nonLinearType = atoi(equals + 1);
        break;
      }
    }
  }
  int saveArgc = argc;
  if (info->numberRows != -1234567)
    memset(info, 0, sizeof(ampl_info)); // overwrite unless magic number set
  /* save so can be accessed by decodePhrase */
  saveInfo = info;
  info->numberArguments = 0;
  info->arguments = (char **)malloc(2 * sizeof(char *));
  info->arguments[info->numberArguments++] = strdup("ampl");
  info->arguments[info->numberArguments++] = strdup(solvername);
  asl = ASL_alloc(ASL_read_f);
  stub = getstub(&argv, &Oinfo);
  if (!stub)
    usage_ASL(&Oinfo, 1);
  Oinfo.sname = strdup(solvername); /* invocation name of solver */
  Oinfo.bsname = strdup(solvername); /* solver name in startup "banner" */
  Oinfo.opname = strdup(tempBuffer); /* name of solver_options environment var */
  nl = jac0dim(stub, 0);
  suf_declare(suftab, sizeof(suftab) / sizeof(SufDecl));

  /* set A_vals to get the constraints column-wise (malloc so can be freed) */
  A_vals = (double *)malloc(nzc * sizeof(double));
  if (!A_vals) {
    printf("no memory\n");
    return 1;
  }
  /* say we want primal solution */
  want_xpi0 = 1;
  /* for basis info */
  info->columnStatus = (int *)malloc(n_var * sizeof(int));
  for (int i = 0; i < n_var; i++)
    info->columnStatus[i] = 3;
  info->rowStatus = (int *)malloc(n_con * sizeof(int));
  for (int i = 0; i < n_con; i++)
    info->rowStatus[i] = 1;
  csd = suf_iput("sstatus", ASL_Sufkind_var, info->columnStatus);
  rsd = suf_iput("sstatus", ASL_Sufkind_con, info->rowStatus);
  if (!(nlvc + nlvo) && nonLinearType < 10) {
    /* read linear model*/
    if( sizeof(CoinBigIndex) == sizeof(int) )
      f_read(nl, 0);
    else
      f_read(nl, ASL_allow_Z | ASL_use_Z);
    // see if any sos
    if (true) {
      char *sostype;
      int nsosnz, *sosbeg, *sosind, *sospri;
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
        info->numberSos = nsos;
        info->sosType = (char *)malloc(nsos);
        info->sosPriority = (int *)malloc(nsos * sizeof(int));
        info->sosStart = (int *)malloc((nsos + 1) * sizeof(int));
        info->sosIndices = (int *)malloc(nsosnz * sizeof(int));
        info->sosReference = (double *)malloc(nsosnz * sizeof(double));
        sos_kludge(nsos, sosbeg, sosref, sosind);
        for (int i = 0; i < nsos; i++) {
          char ichar = sostype[i];
          assert(ichar == '1' || ichar == '2');
          info->sosType[i] = static_cast< char >(ichar - '0');
        }
        memcpy(info->sosPriority, sospri, nsos * sizeof(int));
        memcpy(info->sosStart, sosbeg, (nsos + 1) * sizeof(int));
        memcpy(info->sosIndices, sosind, nsosnz * sizeof(int));
        memcpy(info->sosReference, sosref, nsosnz * sizeof(double));
      }
    }

    /*sos_finish(&specialOrderedInfo, 0, &j, 0, 0, 0, 0, 0);*/
    Oinfo.uinfo = tempBuffer;
    if (getopts(argv, &Oinfo))
      return 1;
    /* objective*/
    obj = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++)
      obj[i] = 0.0;
    if (n_obj) {
      for (og = Ograd[0]; og; og = og->next)
        obj[og->varno] = og->coef;
    }
    if (objtype[0])
      info->direction = -1.0;
    else
      info->direction = 1.0;
    info->offset = objconst(0);
    /* Column bounds*/
    columnLower = (double *)malloc(n_var * sizeof(double));
    columnUpper = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++) {
      columnLower[i] = LUv[2 * i];
      if (columnLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        columnLower[i] = -COIN_DBL_MAX;
      columnUpper[i] = LUv[2 * i + 1];
      if (columnUpper[i] >= COIN_DBL_MAX /* Infinity */)
        columnUpper[i] = COIN_DBL_MAX;
    }
    /* Row bounds*/
    rowLower = (double *)malloc(n_con * sizeof(double));
    rowUpper = (double *)malloc(n_con * sizeof(double));
    for (i = 0; i < n_con; i++) {
      rowLower[i] = LUrhs[2 * i];
      if (rowLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        rowLower[i] = -COIN_DBL_MAX;
      rowUpper[i] = LUrhs[2 * i + 1];
      if (rowUpper[i] >= COIN_DBL_MAX /* Infinity */)
        rowUpper[i] = COIN_DBL_MAX;
    }
    info->numberRows = n_con;
    info->numberColumns = n_var;
    info->numberElements = nzc;
    info->numberBinary = nbv;
    info->numberIntegers = niv + nbv;
    info->objective = obj;
    info->rowLower = rowLower;
    info->rowUpper = rowUpper;
    info->columnLower = columnLower;
    info->columnUpper = columnUpper;
    if( sizeof(CoinBigIndex) == sizeof(int) )
      info->starts = reinterpret_cast<unsigned COINUTILS_BIGINDEX_T*>(A_colstarts);
    else
      info->starts = reinterpret_cast<unsigned COINUTILS_BIGINDEX_T*>(A_colstartsZ);  // FIXME doesn't this assume COINUTILS_BIGINDEX_T = long?
    /*A_colstarts=NULL;*/
    info->rows = A_rownos;
    /*A_rownos=NULL;*/
    info->elements = A_vals;
    /*A_vals=NULL;*/
    info->primalSolution = NULL;
    /* put in primalSolution if exists */
    if (X0) {
      info->primalSolution = (double *)malloc(n_var * sizeof(double));
      memcpy(info->primalSolution, X0, n_var * sizeof(double));
    }
    info->dualSolution = NULL;
    if (niv + nbv > 0)
      mip_stuff(); // get any extra info
    if ((!(niv + nbv) && (csd->kind & ASL_Sufkind_input))
      || (rsd->kind & ASL_Sufkind_input)) {
      /* convert status - need info on map */
      static int map[] = { 1, 3, 1, 1, 2, 1, 1 };
      stat_map(info->columnStatus, n_var, map, 6, "incoming columnStatus");
      stat_map(info->rowStatus, n_con, map, 6, "incoming rowStatus");
    } else {
      /* all slack basis */
      // leave status for output */
#ifdef JJF_ZERO
      free(info->rowStatus);
      info->rowStatus = NULL;
      free(info->columnStatus);
      info->columnStatus = NULL;
#endif
    }
  } else {
    // QP
    // Add .nl if not there
    if (!strstr(fileName, ".nl"))
      strcat(fileName, ".nl");
    CoinModel *model = new CoinModel((nonLinearType > 10) ? 2 : 1, fileName, info);
    if (model->numberRows() > 0 || model->numberColumns() > 0)
      *coinModel = (void *)model;
    Oinfo.uinfo = tempBuffer;
    if (getopts(argv, &Oinfo))
      return 1;
    Oinfo.wantsol = 1;
    if (objtype[0])
      info->direction = -1.0;
    else
      info->direction = 1.0;
    model->setOptimizationDirection(info->direction);
    info->offset = objconst(0);
    info->numberRows = n_con;
    info->numberColumns = n_var;
    info->numberElements = nzc;
    info->numberBinary = nbv;
    int numberIntegers = niv + nlvci + nlvoi + nbv;
    if (nlvci + nlvoi + nlvc + nlvo) {
      // Non linear
      // No idea if there are overlaps so compute
      int numberIntegers = 0;
      for (i = 0; i < n_var; i++) {
        if (model->columnIsInteger(i))
          numberIntegers++;
      }
    }
    info->numberIntegers = numberIntegers;
    // Say nonlinear if it is
    info->nonLinear = nlvc + nlvo;
    if (numberIntegers > 0) {
      mip_stuff(); // get any extra info
      if (info->cut)
        model->setCutMarker(info->numberRows, info->cut);
      if (info->priorities)
        model->setPriorities(info->numberColumns, info->priorities);
    }
  }
  /* add -solve - unless something there already
     - also check for sleep=yes */
  {
    int found = 0;
    int foundLog = 0;
    int foundSleep = 0;
    const char *something[] = { "solve", "branch", "duals", "primals", "user" };
    for (i = 0; i < info->numberArguments; i++) {
      unsigned int j;
      const char *argument = info->arguments[i];
      for (j = 0; j < sizeof(something) / sizeof(char *); j++) {
        const char *check = something[j];
        if (!strncmp(argument, check, strlen(check))) {
          found = (int)(j + 1);
        } else if (!strncmp(argument, "log", 3)) {
          foundLog = 1;
        } else if (!strncmp(argument, "sleep", 5)) {
          foundSleep = 1;
        }
      }
    }
    if (foundLog) {
      /* print options etc */
      for (i = 0; i < saveArgc; i++)
        printf("%s ", saveArgv[i]);
      printf("\n");
      if (environment)
        printf("env %s\n", environment);
      /*printf("%d rows %d columns %d elements\n",n_con,n_var,nzc);*/
    }
    if (!found) {
      if (!strlen(algFound)) {
        info->arguments = (char **)realloc(info->arguments, (info->numberArguments + 1) * sizeof(char *));
        info->arguments[info->numberArguments++] = strdup("-solve");
      } else {
        // use algorithm from keyword
        info->arguments = (char **)realloc(info->arguments, (info->numberArguments + 1) * sizeof(char *));
        info->arguments[info->numberArguments++] = strdup(algFound);
      }
    }
    if (foundSleep) {
      /* let user copy .nl file */
      fprintf(stderr, "You can copy .nl file %s for debug purposes or attach debugger\n", saveArgv[1]);
      fprintf(stderr, "Type q to quit, anything else to continue\n");
      int getChar = getc(stdin);
      if (getChar == 'q' || getChar == 'Q')
        exit(1);
    }
  }
  /* add -quit */
  info->arguments = (char **)realloc(info->arguments, (info->numberArguments + 1) * sizeof(char *));
  info->arguments[info->numberArguments++] = strdup("-quit");
  return 0;
}
void freeArrays1(ampl_info *info)
{
  free(info->objective);
  info->objective = NULL;
  free(info->rowLower);
  info->rowLower = NULL;
  free(info->rowUpper);
  info->rowUpper = NULL;
  free(info->columnLower);
  info->columnLower = NULL;
  free(info->columnUpper);
  info->columnUpper = NULL;
  /* this one not freed by ASL_free */
  free(info->elements);
  info->elements = NULL;
  free(info->primalSolution);
  info->primalSolution = NULL;
  free(info->dualSolution);
  info->dualSolution = NULL;
  /*free(info->rowStatus);
    info->rowStatus=NULL;
    free(info->columnStatus);
    info->columnStatus=NULL;*/
}
void freeArrays2(ampl_info *info)
{
  free(info->primalSolution);
  info->primalSolution = NULL;
  free(info->dualSolution);
  info->dualSolution = NULL;
  free(info->rowStatus);
  info->rowStatus = NULL;
  free(info->columnStatus);
  info->columnStatus = NULL;
  free(info->priorities);
  info->priorities = NULL;
  free(info->branchDirection);
  info->branchDirection = NULL;
  free(info->pseudoDown);
  info->pseudoDown = NULL;
  free(info->pseudoUp);
  info->pseudoUp = NULL;
  free(info->sosType);
  info->sosType = NULL;
  free(info->sosPriority);
  info->sosPriority = NULL;
  free(info->sosStart);
  info->sosStart = NULL;
  free(info->sosIndices);
  info->sosIndices = NULL;
  free(info->sosReference);
  info->sosReference = NULL;
  free(info->cut);
  info->cut = NULL;
  ASL_free(&asl);
  
  free(Oinfo.sname);
  Oinfo.sname = NULL;
  free(Oinfo.bsname);
  Oinfo.bsname = NULL;
  free(Oinfo.opname);
  Oinfo.opname = NULL;
}
void freeArgs(ampl_info *info)
{
  int i;
  for (i = 0; i < info->numberArguments; i++)
    free(info->arguments[i]);
  free(info->arguments);
}
int ampl_obj_prec()
{
  int precision = obj_prec();
  if (precision <= 0)
    precision = 15;
  return precision;
}
void writeAmpl(ampl_info *info)
{
  char buf[1000];
  typedef struct {
    const char *msg;
    int code;
    int wantObj;
  } Sol_info;
  static Sol_info solinfo[] = {
    { "optimal solution", 000, 1 },
    { "infeasible", 200, 1 },
    { "unbounded", 300, 0 },
    { "iteration limit etc", 400, 1 },
    { "solution limit", 401, 1 },
    { "ran out of space", 500, 0 },
    { "status unknown", 501, 1 },
    { "bug!", 502, 0 },
    { "best MIP solution so far restored", 101, 1 },
    { "failed to restore best MIP solution", 503, 1 },
    { "optimal (?) solution", 100, 1 }
  };
  /* convert status - need info on map */
  static int map[] = { 0, 3, 4, 1 };
  sprintf(buf, "%s %s", Oinfo.bsname, info->buffer);
  solve_result_num = solinfo[info->problemStatus].code;
  if (info->columnStatus) {
    stat_map(info->columnStatus, n_var, map, 4, "outgoing columnStatus");
    stat_map(info->rowStatus, n_con, map, 4, "outgoing rowStatus");
    suf_iput("sstatus", ASL_Sufkind_var, info->columnStatus);
    suf_iput("sstatus", ASL_Sufkind_con, info->rowStatus);
  }
  write_sol(buf, info->primalSolution, info->dualSolution, &Oinfo);
}
/* Read a problem from AMPL nl file
 */
CoinModel::CoinModel(int nonLinear, const char *fileName, const void *info)
  : CoinBaseModel()
  , maximumRows_(0)
  , maximumColumns_(0)
  , numberElements_(0)
  , maximumElements_(0)
  , numberQuadraticElements_(0)
  , maximumQuadraticElements_(0)
  , rowLower_(NULL)
  , rowUpper_(NULL)
  , rowType_(NULL)
  , objective_(NULL)
  , columnLower_(NULL)
  , columnUpper_(NULL)
  , integerType_(NULL)
  , columnType_(NULL)
  , start_(NULL)
  , elements_(NULL)
  , packedMatrix_(NULL)
  , quadraticElements_(NULL)
  , sortIndices_(NULL)
  , sortElements_(NULL)
  , sortSize_(0)
  , sizeAssociated_(0)
  , associated_(NULL)
  , numberSOS_(0)
  , startSOS_(NULL)
  , memberSOS_(NULL)
  , typeSOS_(NULL)
  , prioritySOS_(NULL)
  , referenceSOS_(NULL)
  , priority_(NULL)
  , cut_(NULL)
  , moreInfo_(NULL)
  , type_(-1)
  , noNames_(false)
  , links_(0)
{
  problemName_ = "";
  int status = 0;
  if (!strcmp(fileName, "-") || !strcmp(fileName, "stdin")) {
    // stdin
  } else {
    std::string name = fileName;
    bool readable = fileCoinReadable(name);
    if (!readable) {
      std::cerr << "Unable to open file "
                << fileName << std::endl;
      status = -1;
    }
  }
  if (!status) {
    gdb(nonLinear, fileName, info);
  }
}
#ifdef JJF_ZERO
static real
qterm(ASL *asl, fint *colq, fint *rowq, real *delsq)
{
  double t, t1, *x, *x0, *xe;
  fint *rq0, *rqe;

  t = 0.;
  x0 = x = X0;
  xe = x + n_var;
  rq0 = rowq;
  while (x < xe) {
    t1 = *x++;
    rqe = rq0 + *++colq;
    while (rowq < rqe)
      t += t1 * x0[*rowq++] * *delsq++;
  }
  return 0.5 * t;
}
#endif
// stolen from IPopt with changes
typedef struct {
  double obj_sign_;
  ASL_pfgh *asl_;
  double *non_const_x_;
  int *column_; // for jacobian
  int *rowStart_;
  double *gradient_;
  double *constraintValues_;
  int nz_h_full_; // number of nonzeros in hessian
  int nerror_;
  bool objval_called_with_current_x_;
  bool conval_called_with_current_x_;
  bool jacval_called_with_current_x_;
} ClpAmplInfo;

void CoinModel::gdb(int nonLinear, const char *fileName, const void *info)
{
  const ampl_info *amplInfo = (const ampl_info *)info;
  ograd *og = NULL;
  int i;
  SufDesc *csd = NULL;
  SufDesc *rsd = NULL;
  /*bool *basis, *lower;*/
  /*double *LU, *c, lb, objadj, *rshift, *shift, t, ub, *x, *x0, *x1;*/
  //char tempBuffer[20];
  double *objective = NULL;
  double *columnLower = NULL;
  double *columnUpper = NULL;
  double *rowLower = NULL;
  double *rowUpper = NULL;
  int *columnStatus = NULL;
  int *rowStatus = NULL;
  int numberRows = -1;
  int numberColumns = -1;
  int numberElements = -1;
  int numberBinary = -1;
  int numberIntegers = -1;
  int numberAllNonLinearBoth = 0;
  int numberIntegerNonLinearBoth = 0;
  int numberAllNonLinearConstraints = 0;
  int numberIntegerNonLinearConstraints = 0;
  int numberAllNonLinearObjective = 0;
  int numberIntegerNonLinearObjective = 0;
  double *primalSolution = NULL;
  double direction = 1.0;
  char *stub = strdup(fileName);
  CoinPackedMatrix matrixByRow;
  fint **colqp = NULL;
  int *z = NULL;
  if (nonLinear == 0) {
    // linear
    asl = ASL_alloc(ASL_read_f);
    nl = jac0dim(stub, 0);
    free(stub);
    suf_declare(suftab, sizeof(suftab) / sizeof(SufDecl));

    /* set A_vals to get the constraints column-wise (malloc so can be freed) */
    A_vals = (double *)malloc(nzc * sizeof(double));
    if (!A_vals) {
      printf("no memory\n");
      return;
    }
    /* say we want primal solution */
    want_xpi0 = 1;
    /* for basis info */
    columnStatus = (int *)malloc(n_var * sizeof(int));
    rowStatus = (int *)malloc(n_con * sizeof(int));
    csd = suf_iput("sstatus", ASL_Sufkind_var, columnStatus);
    rsd = suf_iput("sstatus", ASL_Sufkind_con, rowStatus);
    /* read linear model*/
    if( sizeof(CoinBigIndex) == sizeof(int) )
      f_read(nl, 0);
    else
      f_read(nl, ASL_allow_Z | ASL_use_Z);
    // see if any sos
    if (true) {
      char *sostype;
      int nsosnz, *sosbeg, *sosind, *sospri;
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
#ifdef JJF_ZERO
        info->numberSos = nsos;
        info->sosType = (char *)malloc(nsos);
        info->sosPriority = (int *)malloc(nsos * sizeof(int));
        info->sosStart = (int *)malloc((nsos + 1) * sizeof(int));
        info->sosIndices = (int *)malloc(nsosnz * sizeof(int));
        info->sosReference = (double *)malloc(nsosnz * sizeof(double));
        sos_kludge(nsos, sosbeg, sosref, sosind);
        for (int i = 0; i < nsos; i++) {
          int ichar = sostype[i];
          assert(ichar == '1' || ichar == '2');
          info->sosType[i] = ichar - '0';
        }
        memcpy(info->sosPriority, sospri, nsos * sizeof(int));
        memcpy(info->sosStart, sosbeg, (nsos + 1) * sizeof(int));
        memcpy(info->sosIndices, sosind, nsosnz * sizeof(int));
        memcpy(info->sosReference, sosref, nsosnz * sizeof(double));
#endif
      }
    }

    /*sos_finish(&specialOrderedInfo, 0, &j, 0, 0, 0, 0, 0);*/
    //Oinfo.uinfo = tempBuffer;
    //if (getopts(argv, &Oinfo))
    //return 1;
    /* objective*/
    objective = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++)
      objective[i] = 0.0;
    if (n_obj) {
      for (og = Ograd[0]; og; og = og->next)
        objective[og->varno] = og->coef;
    }
    if (objtype[0])
      direction = -1.0;
    else
      direction = 1.0;
    objectiveOffset_ = objconst(0);
    /* Column bounds*/
    columnLower = (double *)malloc(n_var * sizeof(double));
    columnUpper = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++) {
      columnLower[i] = LUv[2 * i];
      if (columnLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        columnLower[i] = -COIN_DBL_MAX;
      columnUpper[i] = LUv[2 * i + 1];
      if (columnUpper[i] >= COIN_DBL_MAX /* Infinity */)
        columnUpper[i] = COIN_DBL_MAX;
    }
    /* Row bounds*/
    rowLower = (double *)malloc(n_con * sizeof(double));
    rowUpper = (double *)malloc(n_con * sizeof(double));
    for (i = 0; i < n_con; i++) {
      rowLower[i] = LUrhs[2 * i];
      if (rowLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        rowLower[i] = -COIN_DBL_MAX;
      rowUpper[i] = LUrhs[2 * i + 1];
      if (rowUpper[i] >= COIN_DBL_MAX /* Infinity */)
        rowUpper[i] = COIN_DBL_MAX;
    }
    numberRows = n_con;
    numberColumns = n_var;
    numberElements = nzc;
    numberBinary = nbv;
    numberIntegers = niv;
    /* put in primalSolution if exists */
    if (X0) {
      primalSolution = (double *)malloc(n_var * sizeof(double));
      memcpy(primalSolution, X0, n_var * sizeof(double));
    }
    //double * dualSolution=NULL;
    if (niv + nbv > 0)
      mip_stuff(); // get any extra info
    if ((!(niv + nbv) && (csd->kind & ASL_Sufkind_input))
      || (rsd->kind & ASL_Sufkind_input)) {
      /* convert status - need info on map */
      static int map[] = { 1, 3, 1, 1, 2, 1, 1 };
      stat_map(columnStatus, n_var, map, 6, "incoming columnStatus");
      stat_map(rowStatus, n_con, map, 6, "incoming rowStatus");
    } else {
      /* all slack basis */
      // leave status for output */
#ifdef JJF_ZERO
      free(rowStatus);
      rowStatus = NULL;
      free(columnStatus);
      columnStatus = NULL;
#endif
    }
    CoinPackedMatrix columnCopy(true, numberRows, numberColumns, numberElements,
      A_vals, A_rownos,
#if COINUTILS_BIGINDEX_IS_INT
      A_colstarts,
#else
      reinterpret_cast<const CoinBigIndex*>(A_colstartsZ),   // FIXME doesn't that assume CoinBigIndex=long?
#endif
      NULL);
    matrixByRow.reverseOrderedCopyOf(columnCopy);
  } else if (nonLinear == 1) {
    // quadratic
    asl = ASL_alloc(ASL_read_fg);
    nl = jac0dim(stub, (ftnlen)strlen(stub));
    free(stub);
    suf_declare(suftab, sizeof(suftab) / sizeof(SufDecl));
    /* read  model*/
    X0 = (double *)malloc(n_var * sizeof(double));
    CoinZeroN(X0, n_var);
    qp_read(nl, 0);
    assert(n_obj == 1);
    int nz = 1 + n_con;
    colqp = (fint **)malloc(nz * (2 * sizeof(int *) + sizeof(double *)));
    fint **rowqp = colqp + nz;
    double **delsqp = (double **)(rowqp + nz);
    z = (int *)malloc(nz * sizeof(int));
    for (i = 0; i <= n_con; i++) {
      z[i] = nqpcheck(-i, rowqp + i, colqp + i, delsqp + i);
    }
    qp_opify();
    /* objective*/
    objective = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++)
      objective[i] = 0.0;
    if (n_obj) {
      for (og = Ograd[0]; og; og = og->next)
        objective[og->varno] = og->coef;
    }
    if (objtype[0])
      direction = -1.0;
    else
      direction = 1.0;
    objectiveOffset_ = objconst(0);
    /* Column bounds*/
    columnLower = (double *)malloc(n_var * sizeof(double));
    columnUpper = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++) {
      columnLower[i] = LUv[2 * i];
      if (columnLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        columnLower[i] = -COIN_DBL_MAX;
      columnUpper[i] = LUv[2 * i + 1];
      if (columnUpper[i] >= COIN_DBL_MAX /* Infinity */)
        columnUpper[i] = COIN_DBL_MAX;
    }
    // Build by row from scratch
    //matrixByRow.reserve(n_var,nzc,true);
    // say row orderded
    matrixByRow.transpose();
    /* Row bounds*/
    rowLower = (double *)malloc(n_con * sizeof(double));
    rowUpper = (double *)malloc(n_con * sizeof(double));
    CoinBigIndex *rowStart = new CoinBigIndex[n_con + 1];
    int *column = new int[nzc];
    double *element = new double[nzc];
    rowStart[0] = 0;
    numberElements = 0;
    for (i = 0; i < n_con; i++) {
      rowLower[i] = LUrhs[2 * i];
      if (rowLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        rowLower[i] = -COIN_DBL_MAX;
      rowUpper[i] = LUrhs[2 * i + 1];
      if (rowUpper[i] >= COIN_DBL_MAX /* Infinity */)
        rowUpper[i] = COIN_DBL_MAX;
      for (cgrad *cg = Cgrad[i]; cg; cg = cg->next) {
        column[numberElements] = cg->varno;
        element[numberElements++] = cg->coef;
      }
      rowStart[i + 1] = numberElements;
    }
    assert(numberElements == nzc);
    matrixByRow.appendRows(n_con, rowStart, column, element);
    delete[] rowStart;
    delete[] column;
    delete[] element;
    numberRows = n_con;
    numberColumns = n_var;
    //numberElements=nzc;
    numberBinary = nbv;
    numberIntegers = niv;
    numberAllNonLinearBoth = nlvb;
    numberIntegerNonLinearBoth = nlvbi;
    numberAllNonLinearConstraints = nlvc;
    numberIntegerNonLinearConstraints = nlvci;
    numberAllNonLinearObjective = nlvo;
    numberIntegerNonLinearObjective = nlvoi;
    /* say we want primal solution */
    want_xpi0 = 1;
    //double * dualSolution=NULL;
    // save asl
    // Fix memory leak one day
    ClpAmplInfo *info = new ClpAmplInfo;
    //amplGamsData_ = info;
    info->asl_ = NULL; // as wrong form asl;
    info->nz_h_full_ = -1; // number of nonzeros in hessian
    info->objval_called_with_current_x_ = false;
    info->nerror_ = 0;
    info->obj_sign_ = direction;
    info->conval_called_with_current_x_ = false;
    info->non_const_x_ = NULL;
    info->jacval_called_with_current_x_ = false;
    info->rowStart_ = NULL;
    info->column_ = NULL;
    info->gradient_ = NULL;
    info->constraintValues_ = NULL;
  } else if (nonLinear == 2) {
    // General nonlinear!
    //ASL_pfgh* asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
    asl = ASL_alloc(ASL_read_pfgh);
    nl = jac0dim(stub, (ftnlen)strlen(stub));
    free(stub);
    suf_declare(suftab, sizeof(suftab) / sizeof(SufDecl));
    /* read  model*/
    X0 = (double *)malloc(n_var * sizeof(double));
    CoinZeroN(X0, n_var);
    // code stolen from Ipopt
    int retcode = pfgh_read(nl, ASL_return_read_err | ASL_findgroups);

    switch (retcode) {
    case ASL_readerr_none: {
    } break;
    case ASL_readerr_nofile: {
      printf("Cannot open .nl file\n");
      exit(-1);
    } break;
    case ASL_readerr_nonlin: {
      assert(false); // this better not be an error!
      printf("model involves nonlinearities (ed0read)\n");
      exit(-1);
    } break;
    case ASL_readerr_argerr: {
      printf("user-defined function with bad args\n");
      exit(-1);
    } break;
    case ASL_readerr_unavail: {
      printf("user-defined function not available\n");
      exit(-1);
    } break;
    case ASL_readerr_corrupt: {
      printf("corrupt .nl file\n");
      exit(-1);
    } break;
    case ASL_readerr_bug: {
      printf("bug in .nl reader\n");
      exit(-1);
    } break;
    case ASL_readerr_CLP: {
      printf("ASL error message: \"solver cannot handle CLP extensions\"\n");
      exit(-1);
    } break;
    default: {
      printf("Unknown error in stub file read. retcode = %d\n", retcode);
      exit(-1);
    } break;
    }

    // see "changes" in solvers directory of ampl code...
    hesset(1, 0, 1, 0, nlc);

    assert(n_obj == 1);
    // find the nonzero structure for the hessian
    // parameters to sphsetup:
    int coeff_obj = 1; // coefficient of the objective fn ???
    int mult_supplied = 1; // multipliers will be supplied
    int uptri = 1; // only need the upper triangular part
    // save asl
    // Fix memory leak one day
    ClpAmplInfo *info = new ClpAmplInfo;
    moreInfo_ = (void *)info;
    //amplGamsData_ = info;
    info->asl_ = (ASL_pfgh *)asl;
    // This is not easy to get from ampl so save
    info->nz_h_full_ = sphsetup(-1, coeff_obj, mult_supplied, uptri);
    info->objval_called_with_current_x_ = false;
    info->nerror_ = 0;
    info->obj_sign_ = direction;
    info->conval_called_with_current_x_ = false;
    info->non_const_x_ = NULL;
    info->jacval_called_with_current_x_ = false;
    // Look at nonlinear
    if (nzc) {
      n_conjac[1] = nlc; // just nonlinear
      int *rowStart = new int[nlc + 1];
      info->rowStart_ = rowStart;
      // See how many
      int current_nz = 0;
      for (int i = 0; i < nlc; i++) {
        for (cgrad *cg = Cgrad[i]; cg; cg = cg->next) {
          current_nz++;
        }
      }
      // setup the structure
      int *column = new int[current_nz];
      info->column_ = column;
      current_nz = 0;
      rowStart[0] = 0;
      for (int i = 0; i < nlc; i++) {
        for (cgrad *cg = Cgrad[i]; cg; cg = cg->next) {
          cg->goff = current_nz;
          //iRow[cg->goff] = i ;
          //jCol[cg->goff] = cg->varno + 1;
          column[cg->goff] = cg->varno;
          current_nz++;
        }
        rowStart[i + 1] = current_nz;
      }
      info->gradient_ = new double[nzc];
      info->constraintValues_ = new double[nlc];
    }
    /* objective*/
    objective = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++)
      objective[i] = 0.0;
    if (n_obj) {
      for (og = Ograd[0]; og; og = og->next)
        objective[og->varno] = og->coef;
    }
    if (objtype[0])
      direction = -1.0;
    else
      direction = 1.0;
    objectiveOffset_ = objconst(0);
    /* Column bounds*/
    columnLower = (double *)malloc(n_var * sizeof(double));
    columnUpper = (double *)malloc(n_var * sizeof(double));
    for (i = 0; i < n_var; i++) {
      columnLower[i] = LUv[2 * i];
      if (columnLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        columnLower[i] = -COIN_DBL_MAX;
      columnUpper[i] = LUv[2 * i + 1];
      if (columnUpper[i] >= COIN_DBL_MAX /* Infinity */)
        columnUpper[i] = COIN_DBL_MAX;
    }
    // Build by row from scratch
    //matrixByRow.reserve(n_var,nzc,true);
    // say row orderded
    matrixByRow.transpose();
    CoinBigIndex *rowStart = new CoinBigIndex[n_con + 1];
    int *column = new int[nzc];
    double *element = new double[nzc];
    rowStart[0] = 0;
    numberElements = 0;
    /* Row bounds*/
    rowLower = (double *)malloc(n_con * sizeof(double));
    rowUpper = (double *)malloc(n_con * sizeof(double));
    for (i = 0; i < n_con; i++) {
      rowLower[i] = LUrhs[2 * i];
      if (rowLower[i] <= -COIN_DBL_MAX /* negInfinity */)
        rowLower[i] = -COIN_DBL_MAX;
      rowUpper[i] = LUrhs[2 * i + 1];
      if (rowUpper[i] >= COIN_DBL_MAX /* Infinity */)
        rowUpper[i] = COIN_DBL_MAX;
      for (cgrad *cg = Cgrad[i]; cg; cg = cg->next) {
        column[numberElements] = cg->varno;
        double value = cg->coef;
        if (!value)
          value = -1.2345e-29;
        element[numberElements++] = value;
      }
      rowStart[i + 1] = numberElements;
    }
    assert(numberElements == nzc);
    matrixByRow.appendRows(n_con, rowStart, column, element);
    delete[] rowStart;
    delete[] column;
    delete[] element;
    numberRows = n_con;
    numberColumns = n_var;
    numberElements = nzc;
    numberBinary = nbv;
    numberIntegers = niv;
    numberAllNonLinearBoth = nlvb;
    numberIntegerNonLinearBoth = nlvbi;
    numberAllNonLinearConstraints = nlvc;
    numberIntegerNonLinearConstraints = nlvci;
    numberAllNonLinearObjective = nlvo;
    numberIntegerNonLinearObjective = nlvoi;
    /* say we want primal solution */
    want_xpi0 = 1;
    //double * dualSolution=NULL;
  } else {
    abort();
  }
  // set problem name
  problemName_ = "???";

  // Build by row from scratch
  const double *element = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
  for (i = 0; i < numberRows; i++) {
    addRow(rowLength[i], column + rowStart[i],
      element + rowStart[i], rowLower[i], rowUpper[i]);
  }
  // Now do column part
  for (i = 0; i < numberColumns; i++) {
    setColumnBounds(i, columnLower[i], columnUpper[i]);
    setColumnObjective(i, objective[i]);
  }
  for (i = numberColumns - numberBinary - numberIntegers;
       i < numberColumns; i++) {
    setColumnIsInteger(i, true);
  }
  // and non linear
  for (i = numberAllNonLinearBoth - numberIntegerNonLinearBoth;
       i < numberAllNonLinearBoth; i++) {
    setColumnIsInteger(i, true);
  }
  for (i = numberAllNonLinearConstraints - numberIntegerNonLinearConstraints;
       i < numberAllNonLinearConstraints; i++) {
    setColumnIsInteger(i, true);
  }
  for (i = numberAllNonLinearObjective - numberIntegerNonLinearObjective;
       i < numberAllNonLinearObjective; i++) {
    setColumnIsInteger(i, true);
  }
  free(columnLower);
  free(columnUpper);
  free(rowLower);
  free(rowUpper);
  free(objective);
  // space for building a row
  char *temp = new char[30 * numberColumns_];
  // do names
  int iRow;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    char name[9];
    Coin8CharacterName('r',iRow,name);
    setRowName(iRow, name);
  }
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    char name[9];
    Coin8CharacterName('c',iColumn,name);
    setColumnName(iColumn, name);
  }
  if (colqp) {
    // add in quadratic
    int nz = 1 + n_con;
    int nOdd = 0;
    fint **rowqp = colqp + nz;
    double **delsqp = (double **)(rowqp + nz);
    for (i = 0; i <= n_con; i++) {
      int nels = z[i];
      if (nels) {
        double *element = delsqp[i];
        int *start = (int *)colqp[i];
        int *row = (int *)rowqp[i];
        if (!element) {
          // odd row - probably not quadratic
          nOdd++;
          continue;
        }
#ifdef JJF_ZERO
        printf("%d quadratic els\n", nels);
        for (int j = 0; j < n_var; j++) {
          for (int k = start[j]; k < start[j + 1]; k++)
            printf("%d %d %g\n", j, row[k], element[k]);
        }
#endif
        if (i) {
          int iRow = i - 1;
          for (int j = 0; j < n_var; j++) {
            for (int k = start[j]; k < start[j + 1]; k++) {
              int kColumn = row[k];
              double value = element[k];
              // ampl gives twice with assumed 0.5
              if (kColumn < j)
                continue;
              else if (kColumn == j)
                value *= 0.5;
              const char *expr = getElementAsString(iRow, j);
              double constant = 0.0;
              bool linear;
              if (expr && strcmp(expr, "Numeric")) {
                linear = false;
              } else {
                constant = getElement(iRow, j);
                linear = true;
              }
              char temp2[30];
              if (value == 1.0)
                sprintf(temp2, "c%7.7d", kColumn);
              else
                sprintf(temp2, "%g*c%7.7d", value, kColumn);
              if (linear) {
                if (!constant)
                  strcpy(temp, temp2);
                else if (value > 0.0)
                  sprintf(temp, "%g+%s", constant, temp2);
                else
                  sprintf(temp, "%g%s", constant, temp2);
              } else {
                if (value > 0.0)
                  sprintf(temp, "%s+%s", expr, temp2);
                else
                  sprintf(temp, "%s%s", expr, temp2);
              }
              assert(static_cast< int >(strlen(temp)) < 30 * numberColumns_);
              setElement(iRow, j, temp);
              if (amplInfo->logLevel > 1)
                printf("el for row %d column c%7.7d is %s\n", iRow, j, temp);
            }
          }
        } else {
          // objective
          for (int j = 0; j < n_var; j++) {
            for (int k = start[j]; k < start[j + 1]; k++) {
              int kColumn = row[k];
              double value = element[k];
              // ampl gives twice with assumed 0.5
              if (kColumn < j)
                continue;
              else if (kColumn == j)
                value *= 0.5;
              const char *expr = getColumnObjectiveAsString(j);
              double constant = 0.0;
              bool linear;
              if (expr && strcmp(expr, "Numeric")) {
                linear = false;
              } else {
                constant = getColumnObjective(j);
                linear = true;
              }
              char temp2[30];
              if (value == 1.0)
                sprintf(temp2, "c%7.7d", kColumn);
              else
                sprintf(temp2, "%g*c%7.7d", value, kColumn);
              if (linear) {
                if (!constant)
                  strcpy(temp, temp2);
                else if (value > 0.0)
                  sprintf(temp, "%g+%s", constant, temp2);
                else
                  sprintf(temp, "%g%s", constant, temp2);
              } else {
                if (value > 0.0)
                  sprintf(temp, "%s+%s", expr, temp2);
                else
                  sprintf(temp, "%s%s", expr, temp2);
              }
              assert(static_cast< int >(strlen(temp)) < 30 * numberColumns_);
              setObjective(j, temp);
              if (amplInfo->logLevel > 1)
                printf("el for objective column c%7.7d is %s\n", j, temp);
            }
          }
        }
      }
    }
    if (nOdd) {
      printf("%d non-linear constraints could not be converted to quadratic\n", nOdd);
      exit(77);
    }
  }
  delete[] temp;
  free(colqp);
  free(z);
  // see if any sos
  {
    char *sostype;
    int nsosnz, *sosbeg, *sosind, *sospri;
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
      numberSOS_ = nsos;
      typeSOS_ = new int[numberSOS_];
      prioritySOS_ = new int[numberSOS_];
      startSOS_ = new int[numberSOS_ + 1];
      memberSOS_ = new int[nsosnz];
      referenceSOS_ = new double[nsosnz];
      sos_kludge(nsos, sosbeg, sosref, sosind);
      for (int i = 0; i < nsos; i++) {
        int ichar = sostype[i];
        assert(ichar == '1' || ichar == '2');
        typeSOS_[i] = ichar - '0';
      }
      memcpy(prioritySOS_, sospri, nsos * sizeof(int));
      memcpy(startSOS_, sosbeg, (nsos + 1) * sizeof(int));
      memcpy(memberSOS_, sosind, nsosnz * sizeof(int));
      memcpy(referenceSOS_, sosref, nsosnz * sizeof(double));
    }
  }
}
#else
int readAmpl(ampl_info *, int, char **, void **, const char* solvername)
{
  return 0;
}
void freeArrays1(ampl_info *)
{
}
void freeArrays2(ampl_info *)
{
}
void freeArgs(ampl_info *)
{
}
int ampl_obj_prec()
{
  return 0;
}
void writeAmpl(ampl_info *)
{
}
#endif
// Create a name given a sequence number - allows >10000000
inline void CoinModel::Coin8CharacterName(char rowColumn,
					  int number, char * field)
{
  if (number<10000000) {
    sprintf(field,"%c%7.7d",rowColumn,number);
  } else {
    field[0] = rowColumn;
    int put = 8;
    // just using a to z in a fairly random order
    while (number>=26) {
      field[--put] = 'a'+(number%26);
      number /= 26;
    }
    if (number)
      field[--put] = 'a'+(number%26);
    // move up
    int n = 8-put;
    for (int i=0;i<n;i++)
      field[i+1] = field[put++];
    // pad out
    for (int i=n;i<7;i++)
      field[i+1] = '0';
    field[8]='\0';
  }
}
