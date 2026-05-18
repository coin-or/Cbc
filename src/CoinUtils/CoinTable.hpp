// Copyright (C) 2024, COIN-OR Foundation
// All Rights Reserved. This code is published under the Eclipse Public License.

#ifndef CoinTable_H
#define CoinTable_H

/** \file CoinTable.hpp
 *  \brief Lightweight helper for generating fixed-width text tables.
 *
 *  Used by ClpOutput, CbcOutput, and future phases (B&B progress, cut
 *  summaries, presolve summary, etc.) to produce consistent, aligned,
 *  tabular output that works in both ASCII and UTF-8 terminals.
 *
 *  UTF-8 mode uses Unicode box-drawing characters:
 *    ─ (U+2500)  horizontal rule
 *    │ (U+2502)  vertical bar (column separator)
 *    ┬ (U+252C)  top joint
 *    ┼ (U+253C)  crossing joint
 *    ┴ (U+2534)  bottom joint
 *
 *  ASCII fallback:  -  |  -+-  -+-  -+-
 *
 *  Table layout:
 *    indent + col0 + ─┬─ + col1 + ─┬─ + ... + colN-1   (top sep)
 *    indent + col0 + ─┼─ + col1 + ─┼─ + ... + colN-1   (mid sep)
 *    indent + col0 + ─┴─ + col1 + ─┴─ + ... + colN-1   (bot sep)
 *
 *  A "marker row" shows joints for the first few columns then embeds a
 *  label (e.g. "── restart ") and fills the rest with dashes, giving:
 *    "  ─────────┼──────────┼── restart ──────────────────────────────────────"
 *
 *  A "section rule" is a standalone horizontal rule with embedded text:
 *    "── Root LP relaxation (dual simplex) ────────────────────────────────"
 */

#include <string>
#include <vector>

#include "CoinUtilsConfig.h"

class COINUTILSLIB_EXPORT CoinTable {
public:
  /** Separator line type. */
  enum SepType {
    Top,   ///< Uses ┬ joints (opens the table)
    Middle,///< Uses ┼ joints (divides header from data, or sub-sections)
    Bottom ///< Uses ┴ joints (closes the table)
  };

  /** Column definition. */
  struct Col {
    Col() = default;
    Col(std::string n, int w, bool la = false) : name(std::move(n)), width(w), leftAlign(la) {}
    std::string name;      ///< Column header label
    int width;             ///< Visual character width (not including sep chars)
    bool leftAlign = false;///< true = left-align header; false (default) = right-align
  };

  /**
   * @param cols     Column definitions.
   * @param utf8     true = Unicode box-drawing; false = ASCII.
   * @param indent   Number of leading spaces prepended to every generated line.
   * @param compact  true = compact style (no box borders; columns separated by
   *                 spaces; one thin rule under the header).
   */
  CoinTable(const std::vector<Col> &cols, bool utf8 = true, int indent = 2,
    bool compact = false);

  /** Generate a full-width horizontal separator line.
   *
   *  Bordered mode: joints are ─┬─ (Top), ─┼─ (Middle), ─┴─ (Bottom).
   *
   *  Compact mode: @p type is ignored.
   *    UTF-8: one continuous ─────────────────── rule spanning the full width.
   *    ASCII: individual dashes per column separated by 2 spaces (e.g. ------  -----).
   */
  std::string sepLine(SepType type = Middle) const;

  /** Generate the column header row with the column names.
   *
   *  Bordered mode: columns separated by " │ " / " | ".
   *  Compact mode:  columns separated by "  " (two spaces).
   */
  std::string headerLine() const;

  /** Column separator string used between data cells.
   *  Returns " │ " (UTF-8 bordered), " | " (ASCII bordered), or "  " (compact). */
  const char *colSep() const;

  /** Generate a "marker" row.
   *
   *  Shows full-dash segments with cross joints for the first
   *  @p jointsToShow columns, then embeds @p label and fills the
   *  remaining width with dashes.
   *
   *  @param label         The label to embed, e.g. "── restart " (UTF-8)
   *                       or "-- restart " (ASCII).  Should include any
   *                       leading dashes and trailing space desired.
   *  @param jointsToShow  Number of column joints to print before the label.
   *
   *  Example (jointsToShow=2, LP table):
   *    "  ─────────┼──────────┼── restart ─────────────────────────────────"
   */
  std::string markerRow(const std::string &label, int jointsToShow = 2) const;

  /** Total visual width of the table line (indent + all cols + all joints). */
  int totalWidth() const;

  bool utf8() const { return utf8_; }
  bool compact() const { return compact_; }
  int indent() const { return indent_; }
  int numCols() const { return static_cast<int>(cols_.size()); }

  /** Generate a standalone section-rule line with @p text embedded on the left.
   *
   *  Layout:  leaderDashes×hr + " " + text + " " + fill×hr
   *
   *  E.g. sectionRule("Root LP relaxation (dual simplex)", true) →
   *    "── Root LP relaxation (dual simplex) ────────────────────────────────"
   *
   *  If the text itself (plus leader and spaces) fills or exceeds @p totalWidth,
   *  the fill is omitted and the text is returned as-is with the leader.
   *
   *  @note text is assumed to be pure ASCII for width calculation.
   */
  static std::string sectionRule(const std::string &text, bool utf8,
    int totalWidth = 80, int leaderDashes = 2);

  /** @name Phase markers
   *  Use these for consistent, visually distinct phase start/end lines.
   *
   *  phaseStart("Root LP relaxation (dual simplex)", true)  →  "▶ Root LP..."
   *  phaseEnd  ("LP Optimal — Frac: 2/10 ...",       true)  →  "✔ LP Optimal..."
   *
   *  UTF-8 symbols:  ▶ (start)  ✔ (end)  ★ (incumbent)
   *  ASCII fallbacks: >          +         *
   *
   *  dashSep() returns the em-dash separator with spaces (" — " / " - ").
   */
  //@{
  static std::string phaseStart(const std::string &text, bool utf8);
  static std::string phaseEnd(const std::string &text, bool utf8);
  static const char *incumbentMarker(bool utf8);
  static const char *dashSep(bool utf8);
  //@}

private:
  std::vector<Col> cols_;
  bool utf8_;
  bool compact_;
  int indent_;
};

#endif /* CoinTable_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
