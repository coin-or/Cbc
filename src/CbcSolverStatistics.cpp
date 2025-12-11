#include "CbcSolverStatistics.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {
std::string toLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  return value;
}

#include <string>
#include <algorithm>

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

bool fileHasHeader(const std::string &outFileName) {
  std::ifstream in(outFileName.c_str());
  if (!in.good())
    return false;
  std::string line;
  if (!std::getline(in, line))
    return false;
  return !line.empty();
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

  const bool hasHeader = fileHasHeader(outFileName);
  std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::app);
  if (!file.is_open())
    return false;

  if (!hasHeader) {
    std::ostringstream header;
    header << "Name,result,time,sys,elapsed,objective,continuous,tightened,cut_time,";
    header << "nodes,iterations,rows,columns,processed_rows,processed_columns";
    for (int i = 0; i < number_generators; ++i) {
      const char *name = (name_generators && name_generators[i]) ? name_generators[i] : "cut";
      header << ',' << name;
    }
    header << ",runtime_options";
    file << header.str() << '\n';
  }

  std::string inputFileName = parameters[CbcParam::IMPORTFILE]->strVal();

  const std::string problemName = stripExtension(stripPath(inputFileName));

  const std::string runtimeOptions = buildRuntimeOptions(inputQueue);

  std::ostringstream line;
  line << problemName << ',' << result << ','
       << formatDouble(seconds, 2, std::ios_base::fixed) << ','
       << formatDouble(sys_seconds, 2, std::ios_base::fixed) << ','
       << formatDouble(elapsed_seconds, 2, std::ios_base::fixed) << ','
       << formatDouble(obj, 16) << ','
       << formatDouble(continuous, 6) << ','
       << formatDouble(tighter, 6) << ','
       << formatDouble(cut_time, 2, std::ios_base::fixed) << ','
       << nodes << ',' << iterations << ',' << nrows << ',' << ncols << ','
       << nprocessedrows << ',' << nprocessedcols;
  for (int i = 0; i < number_generators; ++i) {
    const int cuts = number_cuts ? number_cuts[i] : 0;
    line << ',' << cuts;
  }
  line << ',' << runtimeOptions;

  file << line.str() << std::endl;

  return true;
}
