#ifndef MIPSTARTIO_HPP_INCLUDED
#define MIPSTARTIO_HPP_INCLUDED

#include <vector>
#include <string>
#include <utility>
class CbcModel;

class OsiSolverInterface;

/* tries to read mipstart (solution file) from
   fileName, filling colValues and obj
   returns 0 with success,
   1 otherwise */
int readMIPStart(CbcModel *model, const char *fileName,
  std::vector< std::pair< std::string, double > > &colValues,
  double &solObj);

/* from a partial list of variables tries to fill the
   remaining variable values.
   extraActions 0 -default, otherwise set integers not mentioned
   1 - to lower bound, 2 - to upper bound
   3,5 ones without costs as 1,2 - ones with costs to cheapest
   4,6 ones without costs as 1,2 - ones with costs to expensive
*/
int computeCompleteSolution(CbcModel *model,
  const std::vector< std::string > colNames,
  const std::vector< std::pair< std::string, double > > &colValues,
  double *sol, double &obj, int extraActions);

#endif // MIPSTARTIO_HPP_INCLUDED

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
