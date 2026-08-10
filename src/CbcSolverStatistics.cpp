#include "CbcSolverStatistics.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string toLower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
    [](unsigned char ch) { return static_cast< char >(std::tolower(ch)); });
  return value;
}

std::string stripExtension(const std::string &filename)
{
  std::string base = filename;

  // First handle optional .gz
  if (base.size() > 3 && base.compare(base.size() - 3, 3, ".gz") == 0) {
    base = base.substr(0, base.size() - 3);
  }

  // Now handle .mps or .lp
  auto endsWith = [&](const std::string &ext) {
    return base.size() >= ext.size() && base.compare(base.size() - ext.size(), ext.size(), ext) == 0;
  };

  if (endsWith(".mps"))
    return base.substr(0, base.size() - 4);

  if (endsWith(".lp"))
    return base.substr(0, base.size() - 3);

  // No recognized extension → return original
  return filename;
}

std::string stripPath(const std::string &value)
{
  std::string::size_type pos = value.find_last_of("/\\");
  if (pos == std::string::npos)
    return value;
  return value.substr(pos + 1);
}

std::string buildRuntimeOptions(const std::deque< std::string > &tokens)
{
  std::ostringstream stream;
  bool first = true;
  for (const std::string &token : tokens) {
    if (token.empty())
      continue;
    if (token == "cbc" || token == "clp")
      continue;
    std::string lower = toLower(token);
    if (lower.find(".mps") != std::string::npos || lower.find(".gz") != std::string::npos)
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

/**
 * Fixed, canonical list of cut generator names, in a stable order.
 *
 * Using a fixed superset (rather than only the generators active in a
 * given run) means the CSV header never needs to change across runs: any
 * generator not present/active for a particular instance simply gets 0
 * cuts / 0.0 seconds in that row. This replaces an earlier scheme that
 * discovered generator columns dynamically from whichever generators were
 * active and rewrote the whole file (header + all prior data rows,
 * padding in zero columns) whenever a run introduced a generator name not
 * seen before -- fragile, and unsafe when multiple cbc processes append to
 * the same file concurrently (e.g. from a parallel test harness).
 *
 * If a new generator name is ever added to CbcSolverCutSetup.cpp that
 * isn't listed here, it will be silently dropped from the CSV -- add its
 * name to this list too.
 */
const std::vector< std::string > &canonicalGeneratorNames()
{
  static const std::vector< std::string > names = {
    "Clique",
    "FlowCover",
    "Gomory",
    "Gomory(2)",
    "GomoryL1",
    "GomoryL2",
    "Knapsack",
    "LiftAndProject",
    "MixedIntegerRounding2",
    "OddWheel",
    "Probing",
    "Reduce-and-split",
    "Reduce-and-split(2)",
    "TwoMirCuts",
    "TwoMirCutsL1",
    "TwoMirCutsL2",
    "ZeroHalf"
  };
  return names;
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
} // namespace

bool CbcSolverStatistics::writeCsv(CbcParameters &parameters,
  const std::string &outFileName,
  const std::deque< std::string > &inputQueue) const
{
  if (outFileName.empty())
    return false;

  // Build a mapping: generator name -> (cuts, time) for the current run.
  std::vector< std::string > currentGenNames;
  std::vector< int > currentGenCuts;
  std::vector< double > currentGenTime;
  for (int i = 0; i < number_generators; ++i) {
    const char *name = (name_generators && name_generators[i])
      ? name_generators[i]
      : "cut";
    currentGenNames.push_back(name);
    currentGenCuts.push_back(number_cuts ? number_cuts[i] : 0);
    currentGenTime.push_back(time_generators ? time_generators[i] : 0.0);
  }

  const std::vector< std::string > &finalGenerators = canonicalGeneratorNames();

  // Build the header string. The generator columns are always the full,
  // fixed canonical set (2 columns each: cuts and time) so the file's
  // column layout never changes across runs/instances.
  std::ostringstream headerStream;
  headerStream << "Name,result,time,sys,elapsed,objective,continuous,"
               << "lp_seconds,tightened,cut_time,"
               << "nodes,iterations,rows,columns,processed_rows,"
               << "processed_columns,cgraph_time,cgraph_density,"
               << "clqstr_extended,clqstr_dominated,clqstr_time,"
               << "coefstr_changed,coefstr_rows,coefstr_time,"
               << "rowred_fixed,rowred_duplicate,rowred_parallel,rowred_time";
  for (const std::string &gn : finalGenerators)
    headerStream << ',' << gn << "_cuts," << gn << "_time";
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
             << ',' << formatDouble(cgraph_density, 6)
             << ',' << clqstr_extended << ',' << clqstr_dominated
             // Six decimals, not the two the other timings use: these
             // pre-root-LP steps are sub-millisecond by design (coefficient
             // tightening averages 0.026 ms over the mip-sanity corpus, 2.0 ms
             // at worst), so a %.2f column would read 0.00 on every instance
             // and the timing would be reported without being recorded.
             << ',' << formatDouble(clqstr_time, 6, std::ios_base::fixed)
             << ',' << coefstr_changed << ',' << coefstr_rows
             << ',' << formatDouble(coefstr_time, 6, std::ios_base::fixed)
             << ',' << rowred_fixed << ',' << rowred_duplicate
             << ',' << rowred_parallel
             // Six decimals for the same reason as the two timings above: this
             // step is sub-millisecond by design.
             << ',' << formatDouble(rowred_time, 6, std::ios_base::fixed);

  // Output generator cuts/time aligned to the canonical set; 0/0.0 for any
  // generator not active in this particular run.
  for (const std::string &gn : finalGenerators) {
    int cuts = 0;
    double time = 0.0;
    for (int i = 0; i < static_cast< int >(currentGenNames.size()); ++i) {
      if (currentGenNames[i] == gn) {
        cuts = currentGenCuts[i];
        time = currentGenTime[i];
        break;
      }
    }
    dataStream << ',' << cuts << ',' << formatDouble(time, 2, std::ios_base::fixed);
  }
  dataStream << ',' << runtimeOptions;
  const std::string newDataLine = dataStream.str();

  // Header is fixed/known in advance, so we only need to write it once
  // (when the file doesn't exist yet) and append afterwards -- safe even
  // when multiple cbc processes append to the same file concurrently
  // (each append is a single, small write).
  std::ifstream probe(outFileName.c_str());
  const bool fileExists = probe.good();
  probe.close();

  std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::app);
  if (!file.is_open())
    return false;
  if (!fileExists)
    file << headerLine << '\n';
  file << newDataLine << std::endl;

  return true;
}
