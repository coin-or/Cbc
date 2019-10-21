/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Guido Tack <guido.tack@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef __MINIZINC_UTILS_H__
#define __MINIZINC_UTILS_H__

#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <limits>
#include <cstring>
#include <iomanip>
#include <cassert>
//#include <chrono>
//#include <ratio>
#include <exception>


#ifdef MZN_HAS_LLROUND
#include <cmath>
namespace MiniZinc {
  inline
  long long int round_to_longlong(double v) {
    return ::llround(v);
  }
}
#else
namespace MiniZinc {
  inline
  long long int round_to_longlong(double v) {
    return static_cast<long long int>(v < 0 ? v-0.5 : v+0.5);
  }
}
#endif

namespace MiniZinc {
  
// #define __MZN_PRINTATONCE__
#ifdef __MZN_PRINTATONCE__
  #define __MZN_PRINT_SRCLOC(e1, e2) \
    std::cerr << '\n' << __FILE__ << ": " << __LINE__ << " (" << __func__ \
     << "): not " << e1 << ":  " << std::flush; \
     std::cerr << e2 << std::endl
#else
  #define __MZN_PRINT_SRCLOC(e1, e2)
#endif
#define MZN_ASSERT_HARD( c ) \
   do { if ( !(c) ) { __MZN_PRINT_SRCLOC( #c, "" ); throw std::runtime_error( #c ); } } while (0)
#define MZN_ASSERT_HARD_MSG( c, e ) \
   do { if ( !(c) ) { __MZN_PRINT_SRCLOC( #c, e ); \
     std::ostringstream oss; oss << "not " << #c << ":  " << e; \
     throw std::runtime_error( oss.str() ); } } while (0)

  inline bool beginswith(std::string s, std::string t) {
    return s.compare(0, t.length(), t)==0;
  }

  inline void checkIOStatus( bool fOk, std::string msg, bool fHard=1 )
  {
    if ( !fOk ) {
#ifdef _MSC_VER
      char errBuf[1024];
      strerror_s(errBuf, sizeof(errBuf), errno);
#else
      char* errBuf = strerror(errno);
#endif
      std::cerr << "\n  " << msg
        << ":   " << errBuf << "." << std::endl;
      MZN_ASSERT_HARD_MSG ( !fHard, msg << ": " << errBuf );
    }
  }
  
  template <class T> inline bool assignStr(T*, const std::string ) { return false; }
  template<> inline bool assignStr(std::string* pS, const std::string s ) {
    *pS = s;
    return true;
  }
  
  /// A simple per-cmdline option parser
  class CLOParser {
    int& i;              // current item
    std::vector<std::string>& argv;
    
  public:
    CLOParser( int& ii, std::vector<std::string>& av )
      : i(ii), argv(av) { }
    template <class Value=int>
    inline bool get(  const char* names, // space-separated option list
                      Value* pResult=nullptr, // pointer to value storage
                      bool fValueOptional=false // if pResult, for non-string values
                ) {
      return getOption( names, pResult, fValueOptional );
    }
    template <class Value=int>
    inline bool getOption(  const char* names, // space-separated option list
                            Value* pResult=nullptr, // pointer to value storage
                            bool fValueOptional=false // if pResult, for non-string values
                ) {
      assert(0 == strchr(names, ','));
      assert(0 == strchr(names, ';'));
      if( i>=argv.size() )
        return false;
      std::string arg( argv[i] );
      /// Separate keywords
      std::string keyword;
      std::istringstream iss( names );
      while ( iss >> keyword ) {
        if ( ((2<keyword.size() || 0==pResult) && arg!=keyword) ||  // exact cmp
          (0!=arg.compare( 0, keyword.size(), keyword )) )           // truncated cmp
          continue;
        /// Process it
        bool combinedArg = false; // whether arg and value are combined in one string (like -Ggecode)
        if ( keyword.size() < arg.size() ) {
          if ( 0==pResult )
            continue;
          combinedArg = true;
          arg.erase( 0, keyword.size() );
        } else {
          if ( 0==pResult )
            return true;
          i++;
          if( i>=argv.size() ) {
            --i;
            return fValueOptional;
          }
          arg = argv[i];
        }
        assert( pResult );
        if ( assignStr( pResult, arg ) ) 
          return true;
        std::istringstream iss( arg );
        Value tmp;
        if ( !( iss >> tmp ) ) {
          if (!combinedArg)
            --i;
          if ( fValueOptional ) {
            return true;
          }
          // Not print because another agent can handle this option
//           cerr << "\nBad value for " << keyword << ": " << arg << endl;
          return false;
        }
        *pResult = tmp;
        return true;
      }
      return false;
    }
  };  // class CLOParser
  
  /// This class prints a value if non-0 and adds comma if not 1st time
  class HadOne {
    bool fHadOne=false;
  public:
    template <class N>
    std::string operator()(const N& val, const char* descr=0) {
      std::ostringstream oss;
      if ( val ) {
        if ( fHadOne )
          oss << ", ";
        fHadOne=true;
        oss << val;
        if ( descr )
          oss << descr;
      }
      return oss.str();
    }
    void reset() { fHadOne=false; }
    operator bool() const { return fHadOne; }
    bool operator!() const { return !fHadOne; }
  };
  
  /// Split a string into words
  /// Add the words into the given vector
  inline void split(const std::string& str, std::vector<std::string>& words) {
    std::istringstream iss(str);
    std::string buf;
    while (iss) {
      iss >> buf;
      words.push_back(buf);
    }
  }
  
  /// Puts the strings' c_str()s into the 2nd argument.
  /// The latter is only valid as long as the former isn't changed.
  inline void vecString2vecPChar(const std::vector<std::string>& vS, std::vector<const char*>& vPC) {
    vPC.resize(vS.size());
    for ( size_t i=0; i<vS.size(); ++i ) {
      vPC[i] = vS[i].c_str();
    }
  }

}

#endif  // __MINIZINC_UTILS_H__

