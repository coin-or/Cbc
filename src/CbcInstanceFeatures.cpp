/**
 * Copyright COIN-OR Foundation (C) 2025
 * All Rights Reserved.
 * This file is distributed under the Eclipse Public License.
 **/

#include "CbcInstanceFeatures.hpp"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "CbcParam.hpp"
#include "OsiFeatures.hpp"
#include "OsiSolverInterface.hpp"

namespace {

std::string stripExtension(const std::string &filename) {
  std::string base = filename;
  if (base.size() > 3 &&
      base.compare(base.size() - 3, 3, ".gz") == 0)
    base = base.substr(0, base.size() - 3);
  if (base.size() >= 4 &&
      base.compare(base.size() - 4, 4, ".mps") == 0)
    return base.substr(0, base.size() - 4);
  if (base.size() >= 3 &&
      base.compare(base.size() - 3, 3, ".lp") == 0)
    return base.substr(0, base.size() - 3);
  return filename;
}

std::string stripPath(const std::string &value) {
  std::string::size_type pos = value.find_last_of("/\\");
  if (pos == std::string::npos)
    return value;
  return value.substr(pos + 1);
}

/** Returns true if the file exists and its first line matches the expected
 *  feature header (i.e. the feature set has not changed). */
bool fileHasMatchingHeader(const std::string &outFileName,
                           const std::string &expectedHeader) {
  std::ifstream in(outFileName.c_str());
  if (!in.good())
    return false;
  std::string firstLine;
  if (!std::getline(in, firstLine))
    return false;
  return firstLine == expectedHeader;
}

std::string buildHeader() {
  std::ostringstream hdr;
  hdr << "Name";
  for (int i = 0; i < OsiFeatures::n; ++i)
    hdr << ',' << OsiFeatures::name(i);
  return hdr.str();
}

} // namespace

CbcInstanceFeatures::CbcInstanceFeatures()
    : values_(NULL), extracted_(false) {
  values_ = new double[OsiFeatures::n];
  memset(values_, 0, sizeof(double) * OsiFeatures::n);
}

CbcInstanceFeatures::~CbcInstanceFeatures() {
  delete[] values_;
}

void CbcInstanceFeatures::extract(OsiSolverInterface *solver) {
  if (!solver)
    return;
  OsiFeatures::compute(values_, solver);
  extracted_ = true;
}

bool CbcInstanceFeatures::writeCsv(CbcParameters &parameters,
                                   const std::string &outFileName) const {
  if (outFileName.empty() || !extracted_)
    return false;

  const std::string header = buildHeader();

  // Determine whether we need to write the header first.
  const bool needHeader = !fileHasMatchingHeader(outFileName, header);

  std::string inputFileName =
      parameters[CbcParam::IMPORTFILE]->strVal();
  const std::string problemName =
      stripExtension(stripPath(inputFileName));

  // Build the data row.
  std::ostringstream row;
  row << problemName;
  for (int i = 0; i < OsiFeatures::n; ++i) {
    row << ',';
    // Use enough precision to round-trip a double.
    std::ostringstream val;
    val << std::setprecision(15) << values_[i];
    row << val.str();
  }

  std::ofstream file(outFileName.c_str(),
                     std::ios::out | std::ios::app);
  if (!file.is_open())
    return false;

  if (needHeader)
    file << header << '\n';

  file << row.str() << '\n';
  return true;
}
