#include "CbcSolverStatistics.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::string toLower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
    [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  return value;
}

std::string stripExtension(const std::string &filename)
{
  std::string base = filename;
  if (base.size() > 3 && base.compare(base.size() - 3, 3, ".gz") == 0)
    base = base.substr(0, base.size() - 3);
  auto endsWith = [&](const std::string &ext) {
    return base.size() >= ext.size() &&
      base.compare(base.size() - ext.size(), ext.size(), ext) == 0;
  };
  if (endsWith(".mps"))
    return base.substr(0, base.size() - 4);
  if (endsWith(".lp"))
    return base.substr(0, base.size() - 3);
  return filename;
}

std::string stripPath(const std::string &value)
{
  std::string::size_type pos = value.find_last_of("/\\");
  if (pos == std::string::npos)
    return value;
  return value.substr(pos + 1);
}

std::string buildRuntimeOptions(const std::deque<std::string> &tokens)
{
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

/** Split a CSV line into fields. */
std::vector<std::string> splitCsv(const std::string &s)
{
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

/**
 * Convert a generator/heuristic name to a safe column-name fragment.
 * Lowercases, replaces runs of non-alphanumeric characters with a single
 * underscore, and strips leading/trailing underscores.
 */
std::string sanitizeName(const std::string &name)
{
  std::string result;
  result.reserve(name.size());
  bool lastWasUnderscore = true; // suppress leading underscore
  for (unsigned char ch : name) {
    if (std::isalnum(ch)) {
      result += static_cast<char>(std::tolower(ch));
      lastWasUnderscore = false;
    } else {
      if (!lastWasUnderscore)
        result += '_';
      lastWasUnderscore = true;
    }
  }
  // strip trailing underscore
  if (!result.empty() && result.back() == '_')
    result.pop_back();
  return result;
}

std::string formatDouble(double value, int precision,
  std::ios_base::fmtflags floatField = std::ios_base::fmtflags(0))
{
  std::ostringstream out;
  if (floatField != std::ios_base::fmtflags(0))
    out.setf(floatField, std::ios_base::floatfield);
  std::streamsize oldPrecision = out.precision();
  out << std::setprecision(precision) << value;
  out.precision(oldPrecision);
  return out.str();
}

/** Number of fixed columns before the per-generator / per-heuristic columns. */
static const int FIXED_COLUMNS = 19;

/**
 * Read existing file contents: header line and all data lines.
 * Returns true if the file existed and contained a header.
 * allExtraColumns is populated with every column name beyond the fixed
 * columns and before "runtime_options".
 */
bool readExistingCsv(const std::string &outFileName,
  std::vector<std::string> &allExtraColumns,
  std::vector<std::string> &dataLines)
{
  allExtraColumns.clear();
  dataLines.clear();
  std::ifstream in(outFileName.c_str());
  if (!in.good())
    return false;
  std::string line;
  if (!std::getline(in, line) || line.empty())
    return false;
  std::vector<std::string> hfields = splitCsv(line);
  if (static_cast<int>(hfields.size()) > FIXED_COLUMNS + 1) {
    for (int i = FIXED_COLUMNS; i < static_cast<int>(hfields.size()) - 1; ++i)
      allExtraColumns.push_back(hfields[i]);
  }
  while (std::getline(in, line))
    if (!line.empty())
      dataLines.push_back(line);
  return true;
}

} // namespace

// ---------------------------------------------------------------------------
// Canonical lists
// ---------------------------------------------------------------------------

const std::vector<std::string> &CbcSolverStatistics::knownCutGenerators()
{
  // Names exactly as passed to CbcModel::addCutGenerator() in
  // CbcSolverCutSetup.cpp.  Add new names here when new generators are
  // introduced so that the CSV column set stays stable.
  static const std::vector<std::string> kList = {
    "Probing",
    "Gomory",
    "GomoryL1",
    "GomoryL2",
    "Gomory(2)",
    "Knapsack",
    "Reduce-and-split",
    "Reduce-and-split(2)",
    "Clique",
    "OddWheel",
    "MixedIntegerRounding2",
    "FlowCover",
    "TwoMirCuts",
    "TwoMirCutsL1",
    "TwoMirCutsL2",
    "LiftAndProject",
    "ResidualCapacity",
    "ZeroHalf",
    "Stored",
  };
  return kList;
}

const std::vector<std::string> &CbcSolverStatistics::knownHeuristics()
{
  // Names exactly as passed to setHeuristicName() in CbcSolverHeuristics.cpp
  // and CbcSolver.cpp.  Names that sanitize to the same string are treated
  // as the same heuristic (e.g. "feasibility pump" and "Feasibility pump").
  static const std::vector<std::string> kList = {
    "feasibility pump",
    "rounding",
    "combine solutions",
    "greedy cover",
    "greedy equality",
    "random rounding",
    "dynamic pass thru",
    "linked",
    "Partial solution given",
    "FeasibilityJump",
    "Dantzig-Wolfe-expansion",
    "RINS",
    "RENS",
    "RENSdj",
    "RENSub",
    "VND",
    "Naive",
    "DiveAny",
    "DiveCoefficient",
    "DiveFractional",
    "DiveGuided",
    "DiveLineSearch",
    "DivePseudoCost",
    "DiveVectorLength",
    "Multiple root solvers",
  };
  return kList;
}

// ---------------------------------------------------------------------------
// writeCsv
// ---------------------------------------------------------------------------

bool CbcSolverStatistics::writeCsv(CbcParameters &parameters,
  const std::string &outFileName,
  const std::deque<std::string> &inputQueue) const
{
  if (outFileName.empty())
    return false;

  // ------------------------------------------------------------------
  // 1. Build the unified cut-generator entity list.
  //    Start with the canonical list, then append any runtime extras
  //    (using sanitized names to detect duplicates).
  // ------------------------------------------------------------------
  std::vector<std::string> cutSanitized; // sanitized names (for column ids)
  auto addCutName = [&](const std::string &name) {
    std::string san = sanitizeName(name);
    for (const auto &s : cutSanitized)
      if (s == san)
        return; // already present
    cutSanitized.push_back(san);
  };
  for (const auto &n : knownCutGenerators())
    addCutName(n);
  for (const auto &cs : cutStats)
    addCutName(cs.name);

  // ------------------------------------------------------------------
  // 2. Build the unified heuristic entity list (same approach).
  // ------------------------------------------------------------------
  std::vector<std::string> heurSanitized;
  auto addHeurName = [&](const std::string &name) {
    std::string san = sanitizeName(name);
    for (const auto &s : heurSanitized)
      if (s == san)
        return;
    heurSanitized.push_back(san);
  };
  for (const auto &n : knownHeuristics())
    addHeurName(n);
  for (const auto &hs : heuristicStats)
    addHeurName(hs.name);

  // ------------------------------------------------------------------
  // 3. Build maps: sanitized-name -> aggregated stats (accumulate
  //    entries with the same sanitized name, e.g. "feasibility pump"
  //    and "Feasibility pump").
  // ------------------------------------------------------------------
  std::unordered_map<std::string, CutGeneratorStats> cutMap;
  for (const auto &cs : cutStats) {
    std::string san = sanitizeName(cs.name);
    auto &acc = cutMap[san];
    if (acc.name.empty())
      acc.name = cs.name; // first entry wins
    acc.nCuts += cs.nCuts;
    acc.nCalls += cs.nCalls;
    acc.time += cs.time;
    acc.nColumnCuts += cs.nColumnCuts;
    if (cs.minDepth >= 0 && (acc.minDepth < 0 || cs.minDepth < acc.minDepth))
      acc.minDepth = cs.minDepth;
    if (cs.maxDepth > acc.maxDepth)
      acc.maxDepth = cs.maxDepth;
  }

  std::unordered_map<std::string, HeuristicStats> heurMap;
  for (const auto &hs : heuristicStats) {
    std::string san = sanitizeName(hs.name);
    auto &acc = heurMap[san];
    if (acc.name.empty())
      acc.name = hs.name; // first entry wins
    acc.nExecutions += hs.nExecutions;
    acc.totalTime += hs.totalTime;
    acc.nSolutions += hs.nSolutions;
    if (hs.minDepth >= 0 && (acc.minDepth < 0 || hs.minDepth < acc.minDepth))
      acc.minDepth = hs.minDepth;
    if (hs.maxDepth > acc.maxDepth)
      acc.maxDepth = hs.maxDepth;
  }

  // ------------------------------------------------------------------
  // 4. Build the full ordered list of extra column names.
  //    For each cut:  cut_<san>_cuts, cut_<san>_calls, cut_<san>_time,
  //                   cut_<san>_minDepth, cut_<san>_maxDepth
  //    For each heur: heur_<san>_execs, heur_<san>_time, heur_<san>_sols,
  //                   heur_<san>_minDepth, heur_<san>_maxDepth
  // ------------------------------------------------------------------
  std::vector<std::string> extraColNames;
  for (const auto &san : cutSanitized) {
    extraColNames.push_back("cut_" + san + "_cuts");
    extraColNames.push_back("cut_" + san + "_calls");
    extraColNames.push_back("cut_" + san + "_time");
    extraColNames.push_back("cut_" + san + "_minDepth");
    extraColNames.push_back("cut_" + san + "_maxDepth");
  }
  for (const auto &san : heurSanitized) {
    extraColNames.push_back("heur_" + san + "_execs");
    extraColNames.push_back("heur_" + san + "_time");
    extraColNames.push_back("heur_" + san + "_sols");
    extraColNames.push_back("heur_" + san + "_minDepth");
    extraColNames.push_back("heur_" + san + "_maxDepth");
  }

  // ------------------------------------------------------------------
  // 5. Read existing file to check if header needs rewriting.
  // ------------------------------------------------------------------
  std::vector<std::string> existingExtraCols;
  std::vector<std::string> existingDataLines;
  bool hadHeader = readExistingCsv(outFileName, existingExtraCols, existingDataLines);

  bool headerChanged = !hadHeader || (existingExtraCols != extraColNames);

  // ------------------------------------------------------------------
  // 6. Build the header string.
  // ------------------------------------------------------------------
  std::ostringstream headerStream;
  headerStream << "Name,result,integer_feasible,time,sys,elapsed,objective,continuous,"
               << "lp_seconds,tightened,cut_time,"
               << "nodes,iterations,rows,columns,processed_rows,"
               << "processed_columns,cgraph_time,cgraph_density";
  for (const auto &col : extraColNames)
    headerStream << ',' << col;
  headerStream << ",runtime_options";
  const std::string headerLine = headerStream.str();

  // ------------------------------------------------------------------
  // 7. Build the new data line.
  // ------------------------------------------------------------------
  std::string inputFileName = parameters[CbcParam::IMPORTFILE]->strVal();
  const std::string problemName = stripExtension(stripPath(inputFileName));
  const std::string runtimeOptions = buildRuntimeOptions(inputQueue);

  std::ostringstream dataStream;
  dataStream << problemName << ',' << result << ','
             << (integer_feasible ? "1" : "0") << ','
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

  for (const auto &san : cutSanitized) {
    const auto it = cutMap.find(san);
    if (it != cutMap.end()) {
      const auto &cs = it->second;
      dataStream << ',' << cs.nCuts
                 << ',' << cs.nCalls
                 << ',' << formatDouble(cs.time, 4, std::ios_base::fixed)
                 << ',' << cs.minDepth
                 << ',' << cs.maxDepth;
    } else {
      dataStream << ",0,0,0.0000,-1,-1";
    }
  }
  for (const auto &san : heurSanitized) {
    const auto it = heurMap.find(san);
    if (it != heurMap.end()) {
      const auto &hs = it->second;
      dataStream << ',' << hs.nExecutions
                 << ',' << formatDouble(hs.totalTime, 4, std::ios_base::fixed)
                 << ',' << hs.nSolutions
                 << ',' << hs.minDepth
                 << ',' << hs.maxDepth;
    } else {
      dataStream << ",0,0.0000,0,-1,-1";
    }
  }
  dataStream << ',' << runtimeOptions;
  const std::string newDataLine = dataStream.str();

  // ------------------------------------------------------------------
  // 8. Write file.
  // ------------------------------------------------------------------
  if (headerChanged) {
    // Build a position map from existing header column names so we can
    // remap old rows to the new column layout.
    std::unordered_map<std::string, int> oldColIndex;
    oldColIndex.reserve(FIXED_COLUMNS + existingExtraCols.size() + 1);
    // Fixed columns get indices 0..FIXED_COLUMNS-1.
    // (We just copy them verbatim, so we only need to map extra cols.)
    for (int i = 0; i < static_cast<int>(existingExtraCols.size()); ++i)
      oldColIndex[existingExtraCols[i]] = FIXED_COLUMNS + i;
    // "runtime_options" is always the last old column.
    int oldRtOptsIdx = FIXED_COLUMNS + static_cast<int>(existingExtraCols.size());

    int newTotalCols = FIXED_COLUMNS + static_cast<int>(extraColNames.size()) + 1;

    std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::trunc);
    if (!file.is_open())
      return false;
    file << headerLine << '\n';

    for (const std::string &dl : existingDataLines) {
      std::vector<std::string> oldFields = splitCsv(dl);
      // Ensure old row has enough fields (pad with empty if truncated).
      while (static_cast<int>(oldFields.size()) < oldRtOptsIdx + 1)
        oldFields.emplace_back("0");

      std::vector<std::string> newFields(newTotalCols, "0");
      // Copy fixed columns verbatim.
      for (int i = 0; i < FIXED_COLUMNS && i < static_cast<int>(oldFields.size()); ++i)
        newFields[i] = oldFields[i];
      // Copy extra columns by name where available.
      for (int i = 0; i < static_cast<int>(extraColNames.size()); ++i) {
        auto it = oldColIndex.find(extraColNames[i]);
        if (it != oldColIndex.end() && it->second < static_cast<int>(oldFields.size()))
          newFields[FIXED_COLUMNS + i] = oldFields[it->second];
      }
      // Last column is runtime_options.
      newFields[newTotalCols - 1] = oldFields[oldRtOptsIdx];

      for (int i = 0; i < newTotalCols; ++i) {
        if (i > 0)
          file << ',';
        file << newFields[i];
      }
      file << '\n';
    }

    file << newDataLine << std::endl;
  } else {
    std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::app);
    if (!file.is_open())
      return false;
    file << newDataLine << std::endl;
  }

  return true;
}
