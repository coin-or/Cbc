#ifndef MIPSTART_HPP_INCLUDED
#define MIPSTART_HPP_INCLUDED

#include <vector>
#include <string>
#include <utility>

#include "CbcConfig.h"

class CbcModel;
class OsiSolverInterface;
class CoinMessageHandler;
class CoinMessages;

class CBCLIB_EXPORT CbcMipStart {
public:
  /*! \brief Read a MIP start file and populate variable/value pairs.

      The file can be whitespace-delimited (default) or use `.csv`/`.psv`
      separators (`,` or `|`). Each record provides a column name and value;
      the reader validates numeric fields, converts text to doubles, and logs
      warnings via \a messHandler for malformed rows. After successful parsing,
      values are clipped to solver bounds when necessary and the resulting
      vector is stored in \a colValues. If the file does not exist or contains
      no valid assignments, an error is reported and a non-zero status is
      returned.

      \param solver      Solver whose column metadata (names/bounds) are used.
      \param fileName    Path to the mipstart file.
      \param colValues   Output vector filled with \c (name,value) pairs.
      \param solObj      Objective value read from the file (if provided).
      \param messHandler Message handler used for status/warning output.
      \param pcoinmsgs   Message catalog used to format emitted messages.

      \return 0 on success, non-zero if the file could not be opened or no
              valid solution was found.
   */
  static int read(OsiSolverInterface *solver, const char *fileName,
    std::vector< std::pair< std::string, double > > &colValues,
    double &solObj, CoinMessageHandler *messHandler, CoinMessages *pcoinmsgs);

  /*! \brief Filter MIP start assignments that respect bounds and integrality.

      Each named value is checked against the current solver model.
      Assignments outside their column bounds (considering the solver's primal
      tolerance) or those belonging to integer/binary columns that violate the
      solver's integer tolerance are discarded. Corresponding warning messages
      are emitted through \a messHandler using the CBC_MIPSTART_* facilities.

      \param solver       Active solver providing names, bounds, and tolerances.
      \param colValues    Raw MIP start assignments (name,value pairs).
      \param messHandler  Message handler used to report discarded values.
      \param pcoinmsgs    Message catalog to format emitted messages.
      \param preProcessedModel    If model is pre-processed and some columns
                                  may have been removed.


      \return Vector containing only the assignments that satisfied the
              solver's bounds and integrality requirements.
   */
  static std::vector< std::pair< std::string, double > > validateMIPStartValues(OsiSolverInterface *solver,
    const std::vector< std::pair< std::string, double > > &colValues,
    CoinMessageHandler *messHandler, CoinMessages *pcoinmsgs,
    bool preProcessedModel = false);

    /*! \brief Build a full feasible solution from partial MIP start values.

      The routine copies the solver, fixes variables to the provided MIP start
      values (respecting bounds and integrality), and optionally sets
      unspecified integer variables to their lower/upper bounds depending on
      \a extraActions. After solving the resulting LP, it may invoke heuristics
      (e.g., small B&B) to repair fractional values. The final solution vector
      and objective are stored in \a sol and \a obj respectively, while
      informational messages are routed through \a messHandler.

      \param model          Optional model used for SOS/object metadata.
      \param solver         Base solver to clone and manipulate.
      \param colNames       Column names aligned with the solver copy.
      \param colValues      Named values from the MIP start.
      \param sol            Output dense solution vector.
      \param obj            Output objective value of the constructed solution.
      \param extraActions   Strategy flag controlling default bounds for
                      unspecified integers (0 = none, 1/2 = force to
                      lower/upper, 3-6 combine bound fixing with cost
                      direction as implemented in CbcMipStart.cpp).
      \param messHandler    Message handler for warnings/status.
      \param pmessages      Message catalog used by \a messHandler.

      \return 0 if a feasible solution is produced, non-zero otherwise.
    */
  static int computeCompleteSolution(CbcModel *model, OsiSolverInterface *solver,
    const std::vector< std::string > &colNames,
    const std::vector< std::pair< std::string, double > > &colValues,
    double *sol, double &obj, int extraActions, CoinMessageHandler *messHandler, CoinMessages *pmessages);
};

#endif // MIPSTART_HPP_INCLUDED

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
