#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <OsiSolverInterface.hpp>
#include "CbcMessage.hpp"
#include <CbcModel.hpp>
#include "CbcMipStartIO.hpp"
#include "CoinTime.hpp"

using namespace std;


bool isNumericStr( const char *str )
{
   const int l = strlen(str);

   for ( int i=0 ; (i<l) ; ++i )
      if (!(isdigit(str[i])||(str[i]=='.')))
         return false;

   return true;
}

int readMIPStart( CbcModel * model, const char *fileName,
                  vector< pair< string, double > > &colValues,
                  double &/*solObj*/ )
{
#define STR_SIZE 256
   FILE *f = fopen( fileName, "r" );
   if (!f)
      return 1;
   char line[STR_SIZE];

   int nLine = 0;
   char printLine[STR_SIZE];
   while (fgets( line, STR_SIZE, f ))
   {
      ++nLine;
      char col[4][STR_SIZE];
      int nread = sscanf( line, "%s %s %s %s", col[0], col[1], col[2], col[3] );
      if (!nread)
         continue;
      /* line with variable value */
      if (strlen(col[0])&&isdigit(col[0][0])&&(nread>=3))
      {
         if (!isNumericStr(col[0]))
         {
            sprintf( printLine, "Reading: %s, line %d - first column in mipstart file should be numeric, ignoring.", fileName, nLine );
	    model->messageHandler()->message(CBC_GENERAL, model->messages())
	      << printLine << CoinMessageEol;
            continue;
         }
         if (!isNumericStr(col[2]))
         {
            sprintf( printLine, "Reading: %s, line %d - Third column in mipstart file should be numeric, ignoring.", fileName, nLine  );
	    model->messageHandler()->message(CBC_GENERAL, model->messages())
	      << printLine << CoinMessageEol;
            continue;
         }

         //int idx = atoi( col[0] );
         char *name = col[1];
         double value = atof( col[2] );
         //double obj = 0.0;
//         if (nread >= 4)
//            obj = atof( col[3] );

         colValues.push_back( pair<string, double>(string(name),value) );
      }
   }

   if (colValues.size()) {
      sprintf( printLine,"mipstart values read for %d variables.", (int)colValues.size());
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< printLine << CoinMessageEol;
   } else
   {
      sprintf( printLine, "No mipstart solution read from %s", fileName );
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< printLine << CoinMessageEol;
      return 1;
   }

   fclose(f);
   return 0;
}

int computeCompleteSolution( CbcModel * model,
                             const vector< string > colNames,
                             const std::vector< std::pair< std::string, double > > &colValues,
                             double *sol, double &obj )
{
   int status = 0;
   double compObj = COIN_DBL_MAX;
   bool foundIntegerSol = false;
   OsiSolverInterface *lp = model->solver()->clone();
   map< string, int > colIdx;
   assert( ((int)colNames.size()) == lp->getNumCols() );

   /* for fast search of column names */
   for ( int i=0 ; (i<(int)colNames.size()) ; ++i )
      colIdx[colNames[i]] = i;

   char printLine[STR_SIZE];
   int fixed = 0;
   int notFound = 0;
   char colNotFound[256] = "";
   for ( int i=0 ; (i<(int)colValues.size()) ; ++i )
   {
      map< string, int >::const_iterator mIt = colIdx.find( colValues[i].first );
      if ( mIt == colIdx.end() )
      {
         if (!notFound)
            strcpy( colNotFound, colValues[i].first.c_str() );
         notFound++;
      }
      else
      {
         const int idx = mIt->second;
         double v = colValues[i].second;
         if (v<1e-8)
            v = 0.0;
         if (lp->isInteger(idx))  // just to avoid small
            v = floor( v+0.5 );   // fractional garbage
         lp->setColBounds( idx, v, v );
         ++fixed;
      }
   }

   if (!fixed)
   {
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< "Warning: MIPstart solution is not valid, ignoring it."
	<< CoinMessageEol;
      goto TERMINATE;
   }

   if ( notFound >= ( ((double)colNames.size()) * 0.5 ) ) {
      sprintf( printLine, "Warning: %d column names were not found (e.g. %s) while filling solution.", notFound, colNotFound );
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< printLine << CoinMessageEol;
   }

   lp->initialSolve();
   if (!lp->isProvenOptimal())
   {
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
      status = 1;
      goto TERMINATE;
   }

   /* some additional effort is needed to provide an integer solution */
   if ( lp->getFractionalIndices().size() > 0 )
   {
      sprintf( printLine,"MIPStart solution provided values for %d of %d integer variables, %d variables are still fractional.", fixed, lp->getNumIntegers(), (int)lp->getFractionalIndices().size() );
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< printLine << CoinMessageEol;
      double start = CoinCpuTime();
      CbcModel babModel( *lp );
      babModel.setLogLevel( 0 );
      babModel.setMaximumNodes( 500 );
      babModel.setMaximumSeconds( 60 );
      babModel.branchAndBound();
      clock_t end = clock();
      if (babModel.bestSolution())
      {
         sprintf( printLine,"Mini branch and bound defined values for remaining variables in %.2f seconds.", 
		  CoinCpuTime()-start);
	 model->messageHandler()->message(CBC_GENERAL, model->messages())
	   << printLine << CoinMessageEol;
         copy( babModel.bestSolution(), babModel.bestSolution()+babModel.getNumCols(), sol );
         foundIntegerSol = true;
         obj = compObj = babModel.getObjValue();
      }
      else
      {
	 model->messageHandler()->message(CBC_GENERAL, model->messages())
	   << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
         status = 1;
         goto TERMINATE;
      }
   }
   else
   {
      foundIntegerSol = true;
      obj = compObj = lp->getObjValue();
      copy( lp->getColSolution(), lp->getColSolution()+lp->getNumCols(), sol );
   }

   if ( foundIntegerSol )
   {
      sprintf( printLine,"mipstart provided solution with cost %g", compObj);
      model->messageHandler()->message(CBC_GENERAL, model->messages())
	<< printLine << CoinMessageEol;
      for ( int i=0 ; (i<lp->getNumCols()) ; ++i )
      {
         if (sol[i]<1e-8)
            sol[i] = 0.0;
         else
            if (lp->isInteger(i))
               sol[i] = floor( sol[i]+0.5 );
      }
   }

TERMINATE:
   delete lp;
   return status;
}
#undef STR_SIZE
