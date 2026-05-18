/** Copyright COIN-OR Foundation (C) 2025
 * All Rights Reserved.
 * This file is distributed under the Eclipse Public License.
 *
 * Purpose: helper structure to map Osi column type codes to human friendly
 * names and descriptions.
 **/

#ifndef OsiColumnType_H
#define OsiColumnType_H

#include "OsiConfig.h"

/** \brief Helper for interpreting the integer codes returned by
 *  OsiSolverInterface::getColType().
 */
struct OSILIB_EXPORT OsiColumnType {
  /// Canonical codes for column types. Stored as char to match existing API.
  enum Code : char {
    Continuous = 0,
    Binary = 1,
    GeneralInteger = 2,
    SemiContinuous = 3,
    SemiInteger = 4
  };

  /// Return a short name for the supplied code.
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

  /// Return a human readable description for the supplied code.
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

  /// Convert from a raw char code to a Code value, defaulting to Continuous.
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

  /// Convenience wrapper returning the canonical short name for a raw code.
  static const char *nameFromChar(char code) { return name(fromChar(code)); }

  /// Convenience wrapper returning the description for a raw code.
  static const char *descriptionFromChar(char code)
  {
    return description(fromChar(code));
  }
};

#endif // OsiColumnType_H
