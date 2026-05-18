// Copyright (C) 2024, COIN-OR Foundation
// All Rights Reserved. This code is published under the Eclipse Public License.

#include "CoinTable.hpp"

#include <cassert>
#include <cstring>
#include <string>
#include <vector>

// ─── UTF-8 box-drawing sequences ──────────────────────────────────────────────
//
//   ─  U+2500  \xe2\x94\x80   horizontal rule
//   │  U+2502  \xe2\x94\x82   vertical bar
//   ┬  U+252C  \xe2\x94\xac   top joint (T-down)
//   ┼  U+253C  \xe2\x94\xbc   crossing joint
//   ┴  U+2534  \xe2\x94\xb4   bottom joint (T-up)

static const char *U8_HR   = "\xe2\x94\x80"; // ─
static const char *U8_VBAR = "\xe2\x94\x82"; // │
static const char *U8_TOP  = "\xe2\x94\xac"; // ┬
static const char *U8_MID  = "\xe2\x94\xbc"; // ┼
static const char *U8_BOT  = "\xe2\x94\xb4"; // ┴

// ─── helpers ──────────────────────────────────────────────────────────────────

/** Count the visual (on-screen) character width of a string.
 *  ASCII bytes count as 1.  UTF-8 multi-byte sequences (each grapheme cluster
 *  we emit is a single code point) count as 1 visual char per sequence. */
static int visualWidth(const std::string &s)
{
  int w = 0;
  for (int k = 0; k < (int)s.size(); ) {
    unsigned char c = (unsigned char)s[k];
    if ((c & 0x80) == 0) {
      ++k;
    } else if ((c & 0xF0) == 0xE0) { // 3-byte sequence (all box-drawing chars)
      k += 3;
    } else if ((c & 0xE0) == 0xC0) { // 2-byte sequence
      k += 2;
    } else if ((c & 0xF8) == 0xF0) { // 4-byte sequence
      k += 4;
    } else {
      ++k; // malformed — advance one byte
    }
    ++w;
  }
  return w;
}

/** Append @p n copies of the horizontal-rule string to @p s. */
static void appendHr(std::string &s, const char *hr, int n)
{
  for (int i = 0; i < n; ++i)
    s += hr;
}

// ─── CoinTable ────────────────────────────────────────────────────────────────

CoinTable::CoinTable(const std::vector<Col> &cols, bool utf8, int indent, bool compact)
  : cols_(cols)
  , utf8_(utf8)
  , compact_(compact)
  , indent_(indent)
{
}

int CoinTable::totalWidth() const
{
  int w = indent_;
  int n = (int)cols_.size();
  for (int i = 0; i < n; ++i) {
    w += cols_[i].width;
    if (i < n - 1)
      w += compact_ ? 1 : 3; // compact: " " (1); bordered: "─┬─" (3 visual)
  }
  return w;
}

std::string CoinTable::sepLine(SepType type) const
{
  const int n = (int)cols_.size();
  const char *hr = utf8_ ? U8_HR : "-";

  if (compact_) {
    // Compact: per-column rule segments separated by one space.
    // UTF-8 uses ─ characters; ASCII uses - characters.
    std::string s(indent_, ' ');
    for (int i = 0; i < n; ++i) {
      appendHr(s, hr, cols_[i].width);
      if (i < n - 1)
        s += " ";
    }
    return s;
  }

  // Bordered mode — type selects the joint character.
  // The joint string (─┬─ / ─┼─ / ─┴─) consists of one leading dash, the
  // joint character, and one trailing dash — total 3 visual chars.
  std::string jt;
  if (utf8_) {
    jt = std::string(U8_HR);
    switch (type) {
    case Top:    jt += U8_TOP; break;
    case Bottom: jt += U8_BOT; break;
    default:     jt += U8_MID; break;
    }
    jt += U8_HR;
  } else {
    jt = "-+-";
  }

  std::string s(indent_, ' ');
  for (int i = 0; i < n; ++i) {
    appendHr(s, hr, cols_[i].width);
    if (i < n - 1)
      s += jt;
  }
  return s;
}

std::string CoinTable::headerLine() const
{
  const int n = (int)cols_.size();
  const char *sep = colSep(); // " │ " / " | " / " "

  std::string s(indent_, ' ');
  for (int i = 0; i < n; ++i) {
    const std::string &nm = cols_[i].name;
    int w = cols_[i].width;
    int pad = w - (int)nm.size();
    if (cols_[i].leftAlign) {
      s += nm;
      if (pad > 0)
        s += std::string(pad, ' ');
    } else {
      if (pad > 0)
        s += std::string(pad, ' ');
      s += nm;
    }
    if (i < n - 1)
      s += sep;
  }
  return s;
}

const char *CoinTable::colSep() const
{
  if (compact_)
    return " ";
  return utf8_ ? " \xe2\x94\x82 " : " | "; // " │ " or " | "
}

std::string CoinTable::markerRow(const std::string &label, int jointsToShow) const
{
  const int n = (int)cols_.size();
  const char *hr = utf8_ ? U8_HR : "-";

  // The joint used inside a marker row is the Middle (┼) joint.
  std::string jt;
  if (utf8_) {
    jt = std::string(U8_HR) + U8_MID + U8_HR;
  } else {
    jt = "-+-";
  }

  std::string s(indent_, ' ');
  int visualPos = indent_;
  const int total = totalWidth();

  int cols_to_show = (jointsToShow < n) ? jointsToShow : n - 1;

  for (int i = 0; i <= cols_to_show; ++i) {
    if (i < cols_to_show) {
      // Print this column's dashes + joint
      appendHr(s, hr, cols_[i].width);
      visualPos += cols_[i].width;
      s += jt; // ─┼─
      visualPos += 3;
    } else {
      // Print the label, then fill remaining width with dashes
      s += label;
      visualPos += visualWidth(label);
      int fill = total - visualPos;
      if (fill > 0)
        appendHr(s, hr, fill);
      break;
    }
  }
  return s;
}

std::string CoinTable::sectionRule(const std::string &text, bool utf8,
  int totalWidth, int leaderDashes)
{
  const char *hr = utf8 ? U8_HR : "-";

  // Layout: leaderDashes×hr + " " + text + " " + fill×hr
  const int fixedVis = leaderDashes + 1 + (int)text.size() + 1;
  const int fillVis = totalWidth - fixedVis;

  std::string s;
  appendHr(s, hr, leaderDashes);
  s += ' ';
  s += text;
  if (fillVis > 0) {
    s += ' ';
    appendHr(s, hr, fillVis - 1); // -1 because we already added the space
  }
  return s;
}

// ---------------------------------------------------------------------------
// Phase markers
// ---------------------------------------------------------------------------

std::string CoinTable::phaseStart(const std::string &text, bool utf8)
{
  // UTF-8: ▶ (U+25B6)   ASCII: >
  return std::string(utf8 ? "\xe2\x96\xb6" : ">") + " " + text;
}

std::string CoinTable::phaseEnd(const std::string &text, bool utf8)
{
  // UTF-8: ✔ (U+2714)   ASCII: +
  return std::string(utf8 ? "\xe2\x9c\x94" : "+") + " " + text;
}

const char *CoinTable::incumbentMarker(bool utf8)
{
  // UTF-8: ★ (U+2605)   ASCII: *
  return utf8 ? "\xe2\x98\x85" : "*";
}

const char *CoinTable::dashSep(bool utf8)
{
  // UTF-8: em-dash with spaces   ASCII: simple dash with spaces
  return utf8 ? " \xe2\x80\x94 " : " - ";
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
