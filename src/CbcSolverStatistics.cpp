#include "CbcSolverStatistics.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string toLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  return value;
}

std::string stripExtension(const std::string &filename) {
    std::string base = filename;

    // First handle optional .gz
    if (base.size() > 3 &&
        base.compare(base.size() - 3, 3, ".gz") == 0) {
        base = base.substr(0, base.size() - 3);
    }

    // Now handle .mps or .lp
    auto endsWith = [&](const std::string &ext) {
        return base.size() >= ext.size() &&
               base.compare(base.size() - ext.size(), ext.size(), ext) == 0;
    };

    if (endsWith(".mps"))
        return base.substr(0, base.size() - 4);

    if (endsWith(".lp"))
        return base.substr(0, base.size() - 3);

    // No recognized extension → return original
    return filename;
}


std::string stripPath(const std::string &value) {
  std::string::size_type pos = value.find_last_of("/\\");
  if (pos == std::string::npos)
    return value;
  return value.substr(pos + 1);
}

std::string buildRuntimeOptions(const std::deque<std::string> &tokens) {
  std::ostringstream stream;
  bool first = true;
  for (const std::string &token : tokens) {
    if (token.empty())
      continue;
    if (token == "cbc" || token == "clp")
      continue;
    std::string lower = toLower(token);
    if (lower.find(".mps") != std::string::npos ||
        lower.find(".gz") != std::string::npos)
      continue;
    if (lower.rfind("-writestat", 0) == 0)
      break;
    if (!first)
      stream << ' ';
    stream << token;
    first = false;
  }
  return stream.str();
}

/** Split a CSV line into fields (no quoting support needed here). */
std::vector<std::string> splitCsv(const std::string &s) {
  std::vector<std::string> fields;
  std::string::size_type start = 0;
  while (true) {
    std::string::size_type pos = s.find(',', start);
    if (pos == std::string::npos) {
      fields.push_back(s.substr(start));
      break;
    }
    fields.push_back(s.substr(start, pos - start));
    start = pos + 1;
  }
  return fields;
}

/** Number of fixed columns before the per-generator columns. */
static const int FIXED_COLUMNS = 18;

/**
 * Read existing file contents: header line and all data lines.
 * Returns true if the file existed and contained a header.
 * headerGenerators will contain the generator column names extracted
 * from the header (columns between the fixed columns and "runtime_options").
 */
bool readExistingCsv(const std::string &outFileName,
                     std::vector<std::string> &headerGenerators,
                     std::vector<std::string> &dataLines) {
  headerGenerators.clear();
  dataLines.clear();
  std::ifstream in(outFileName.c_str());
  if (!in.good())
    return false;
  std::string line;
  if (!std::getline(in, line) || line.empty())
    return false;
  // Parse header to extract generator column names.
  // Header format: <16 fixed columns>,gen1,gen2,...,genN,runtime_options
  std::vector<std::string> hfields = splitCsv(line);
  // Generator columns are between FIXED_COLUMNS and the last column
  // (which is "runtime_options").
  if (static_cast<int>(hfields.size()) > FIXED_COLUMNS + 1) {
    for (int i = FIXED_COLUMNS; i < static_cast<int>(hfields.size()) - 1; ++i)
      headerGenerators.push_back(hfields[i]);
  }
  // Read remaining data lines.
  while (std::getline(in, line))
    if (!line.empty())
      dataLines.push_back(line);
  return true;
}

std::string formatDouble(double value, int precision,
                         std::ios_base::fmtflags floatField = std::ios_base::fmtflags(0)) {
  std::ostringstream out;
  if (floatField != std::ios_base::fmtflags(0))
    out.setf(floatField, std::ios_base::floatfield);
  std::streamsize oldPrecision = out.precision();
  out << std::setprecision(precision) << value;
  out.precision(oldPrecision);
  return out.str();
}
} // namespace

bool CbcSolverStatistics::writeCsv(CbcParameters &parameters,
                                   const std::string &outFileName,
                                   const std::deque<std::string> &inputQueue) const {
  if (outFileName.empty())
    return false;

  // Build a mapping: generator name -> cut count for the current run.
  std::vector<std::string> currentGenNames;
  std::vector<int> currentGenCuts;
  for (int i = 0; i < number_generators; ++i) {
    const char *name = (name_generators && name_generators[i])
                           ? name_generators[i]
                           : "cut";
    currentGenNames.push_back(name);
    currentGenCuts.push_back(number_cuts ? number_cuts[i] : 0);
  }

  // Read existing file (if any) to get the header's generator columns.
  std::vector<std::string> headerGenerators;
  std::vector<std::string> existingDataLines;
  bool hadHeader = readExistingCsv(outFileName, headerGenerators, existingDataLines);

  // Build the final (unified) set of generator column names:
  // start with the header's generators, then append any new generators
  // from the current run that weren't already present.
  std::vector<std::string> finalGenerators = headerGenerators;
  for (const std::string &gn : currentGenNames) {
    bool found = false;
    for (const std::string &hg : finalGenerators) {
      if (hg == gn) { found = true; break; }
    }
    if (!found)
      finalGenerators.push_back(gn);
  }

  // Determine if the header needs to be (re)written because the set of
  // generator columns has changed.
  bool headerChanged = !hadHeader ||
    (finalGenerators.size() != headerGenerators.size());

  // Build the header string.
  std::ostringstream headerStream;
  headerStream << "Name,result,time,sys,elapsed,objective,continuous,"
               << "lp_seconds,tightened,cut_time,"
               << "nodes,iterations,rows,columns,processed_rows,"
               << "processed_columns,cgraph_time,cgraph_density";
  for (const std::string &gn : finalGenerators)
    headerStream << ',' << gn;
  headerStream << ",runtime_options";
  const std::string headerLine = headerStream.str();

  // Build the new data line.
  std::string inputFileName = parameters[CbcParam::IMPORTFILE]->strVal();
  const std::string problemName = stripExtension(stripPath(inputFileName));
  const std::string runtimeOptions = buildRuntimeOptions(inputQueue);

  std::ostringstream dataStream;
  dataStream << problemName << ',' << result << ','
             << formatDouble(seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(sys_seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(elapsed_seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(obj, 16) << ','
             << formatDouble(continuous, 6) << ','
             << formatDouble(lp_seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(tighter, 6) << ','
             << formatDouble(cut_time, 2, std::ios_base::fixed) << ','
             << nodes << ',' << iterations << ',' << nrows << ',' << ncols
             << ',' << nprocessedrows << ',' << nprocessedcols
             << ',' << formatDouble(cgraph_time, 2, std::ios_base::fixed)
             << ',' << formatDouble(cgraph_density, 6);

  // Output generator cut counts aligned to finalGenerators.
  for (const std::string &gn : finalGenerators) {
    int cuts = 0;
    for (int i = 0; i < static_cast<int>(currentGenNames.size()); ++i) {
      if (currentGenNames[i] == gn) {
        cuts = currentGenCuts[i];
        break;
      }
    }
    dataStream << ',' << cuts;
  }
  dataStream << ',' << runtimeOptions;
  const std::string newDataLine = dataStream.str();

  if (headerChanged) {
    // Rewrite the entire file: new header, existing data lines (padded
    // with zeros for any newly-added generator columns), then new line.
    std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::trunc);
    if (!file.is_open())
      return false;

    file << headerLine << '\n';

    // Number of generator columns that were added beyond the old header.
    int addedCols = static_cast<int>(finalGenerators.size())
                    - static_cast<int>(headerGenerators.size());

    for (const std::string &dl : existingDataLines) {
      // Insert `addedCols` zero-valued columns just before the last field
      // (runtime_options).
      std::vector<std::string> fields = splitCsv(dl);
      if (addedCols > 0 && !fields.empty()) {
        // Last field is runtime_options.
        std::string rtOpts = fields.back();
        fields.pop_back();
        for (int j = 0; j < addedCols; ++j)
          fields.push_back("0");
        fields.push_back(rtOpts);
      }
      for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
        if (i > 0) file << ',';
        file << fields[i];
      }
      file << '\n';
    }

    file << newDataLine << std::endl;
  } else {
    // Header is unchanged — just append the new data line.
    std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::app);
    if (!file.is_open())
      return false;
    file << newDataLine << std::endl;
  }

  return true;
}
