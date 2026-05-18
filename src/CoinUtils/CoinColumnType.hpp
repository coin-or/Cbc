/**
 * Helper structure to map column type codes to human friendly
 * names and descriptions. Moved from Osi to CoinUtils.
 */

#ifndef CoinColumnType_H
#define CoinColumnType_H

#include "CoinUtilsConfig.h"

struct COINUTILSLIB_EXPORT CoinColumnType {
  /// Canonical codes for column types. Stored as char to match existing API.
  enum Code : char {
    Continuous = 0,
    Binary = 1,
    GeneralInteger = 2,
    SemiContinuous = 3,
    SemiInteger = 4
  };

  static const char *name(Code code)
  {
    switch (code) {
    case Continuous:
      return "continuous";
    case Binary:
      return "binary";
    case GeneralInteger:
      return "general-integer";
    case SemiContinuous:
      return "semi-continuous";
    case SemiInteger:
      return "semi-integer";
    default:
      return "unknown";
    }
  }

  static const char *description(Code code)
  {
    switch (code) {
    case Continuous:
      return "Continuous column";
    case Binary:
      return "Binary column (bounds in {0,1})";
    case GeneralInteger:
      return "General integer column";
    case SemiContinuous:
      return "Semi-continuous column";
    case SemiInteger:
      return "Semi-continuous integer column";
    default:
      return "Unknown column type";
    }
  }

  static Code fromChar(char code)
  {
    switch (code) {
    case Continuous:
    case Binary:
    case GeneralInteger:
    case SemiContinuous:
    case SemiInteger:
      return static_cast<Code>(code);
    default:
      return Continuous;
    }
  }

  static const char *nameFromChar(char code) { return name(fromChar(code)); }
  static const char *descriptionFromChar(char code) { return description(fromChar(code)); }
};

#endif // CoinColumnType_H
