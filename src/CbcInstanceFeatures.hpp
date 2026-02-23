/**
 * Copyright COIN-OR Foundation (C) 2025
 * All Rights Reserved.
 * This file is distributed under the Eclipse Public License.
 *
 * Purpose: Extract OsiFeatures from a MIP instance and write them to a CSV
 * file, one row per instance, in a format compatible with -writeStatistics.
 **/

#ifndef CBC_INSTANCE_FEATURES_HPP
#define CBC_INSTANCE_FEATURES_HPP

#include <string>
#include "CbcParameters.hpp"

class OsiSolverInterface;

/**
 * Extracts all OsiFeature values from a given solver interface and can
 * write them to a CSV file (one row per instance, appended).
 *
 * Usage:
 *   CbcInstanceFeatures feat;
 *   feat.extract(solver);
 *   feat.writeCsv(parameters, "features.csv");
 */
class CbcInstanceFeatures {
public:
  CbcInstanceFeatures();
  ~CbcInstanceFeatures();

  /**
   * Extract all features from the given solver.
   * Must be called before writeCsv().
   */
  void extract(OsiSolverInterface *solver);

  /**
   * Append the extracted features to a CSV file.
   *
   * The first row is the header (written only when the file is new/empty).
   * Columns: Name, <feature_0>, <feature_1>, ..., <feature_N-1>
   *
   * @param parameters  CBC parameters (used to get IMPORTFILE name).
   * @param outFileName Path to the CSV file.
   * @return true on success, false if the file could not be opened.
   */
  bool writeCsv(CbcParameters &parameters,
                const std::string &outFileName) const;

private:
  double *values_;   ///< Feature values (size OsiFeatures::n after extract())
  bool extracted_;   ///< True after extract() has been called

  // Non-copyable
  CbcInstanceFeatures(const CbcInstanceFeatures &);
  CbcInstanceFeatures &operator=(const CbcInstanceFeatures &);
};

#endif /* CBC_INSTANCE_FEATURES_HPP */
