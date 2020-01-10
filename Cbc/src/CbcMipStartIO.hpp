#ifndef MIPSTARTIO_HPP_INCLUDED
#define MIPSTARTIO_HPP_INCLUDED

#include <vector>
#include <string>
#include <utility>
class CbcModel;
class OsiSolverInterface;
class CoinMessageHandler;
class CoinMessages;

class CbcMipStartIO{
public:  
/* tries to read mipstart (solution file) from
   fileName, filling colValues and obj
   returns 0 with success,
   1 otherwise */
static int read(OsiSolverInterface *solver, const char *fileName,
  std::vector< std::pair< std::string, double > > &colValues,
  double &solObj, CoinMessageHandler *messHandler, CoinMessages *pcoinmsgs);

/* from a partial list of variables tries to fill the
   remaining variable values.
   extraActions 0 -default, otherwise set integers not mentioned
   1 - to lower bound, 2 - to upper bound
   3,5 ones without costs as 1,2 - ones with costs to cheapest
   4,6 ones without costs as 1,2 - ones with costs to expensive
*/
static int computeCompleteSolution(CbcModel *model, OsiSolverInterface *solver,
  const std::vector< std::string > colNames,
  const std::vector< std::pair< std::string, double > > &colValues,
  double *sol, double &obj, int extraActions, CoinMessageHandler *messHandler, CoinMessages *pmessages);
  
};


#endif // MIPSTARTIO_HPP_INCLUDED

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
