// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinUtilsConfig.h"

#include <cassert>
#include <cerrno>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "CoinParam.hpp"
#include "CoinFileIO.hpp"

/* Visible functions */

namespace CoinParamUtils {

//#############################################################################
//#############################################################################

std::string printString(const std::string &input, int maxWidth)
{
   std::ostringstream buffer;
   std::istringstream iss(input);
   int currentWidth = 0;
   std::string word;
   while (iss >> word){
      if (currentWidth >= maxWidth){
         buffer << std::endl << word;
         currentWidth = 0;
      } else {
         buffer << word;
      }         
      currentWidth += word.length();
   }
   return buffer.str();
}

std::string printString(const char *input, int maxWidth)
{
   return CoinParamUtils::printString(std::string(input), maxWidth);
}

//#############################################################################
// Read into the input queue from any stream (file, string, etc.)   
//#############################################################################

void
readFromStream(std::deque<std::string> &inputQueue,
               std::istream &inputStream)
{
   std::string field;
   while (inputStream >> field){
      std::string::size_type found = field.find('=');
      if (found != std::string::npos) {
         inputQueue.push_back(field.substr(0, found));
         inputQueue.push_back(field.substr(found + 1));
      } else {
         inputQueue.push_back(field);
      }
   }
}

//#############################################################################
// Get input from the interactive command prompt and pass it to readFromStrem
//#############################################################################

//#############################################################################
// Get the next field from the input queue   
//#############################################################################

std::string
getNextField(std::deque<std::string> &inputQueue, bool /*interactiveMode*/,
             std::string /*prompt*/)
{
  if (inputQueue.empty()){
     return "";
  }else{
     std::string field = inputQueue.front();
     inputQueue.pop_front();
     const std::string::size_type found = field.find('=');
     if (found != std::string::npos) {
        // Adds a '-' at the beginning of the parameter name if not present.
        const std::string s1 = field.at(0) == '-' ? field.substr(0, found) :
                             '-' + field.substr(0, found);
        const std::string s2 = field.substr(found + 1);
        inputQueue.push_front(s2);
        inputQueue.push_front(s1);

        field = inputQueue.front();
        inputQueue.pop_front();
     }
     return field;
  }
}

//#############################################################################
//#############################################################################

/*
  Function to look up a parameter keyword (name) in the parameter vector and
  deal with the result. The keyword may end in one or more `?' characters;
  this is a query for information about matching parameters.

  If we have a single match satisfying the minimal match requirements, and
  there's no query, we simply return the index of the matching parameter in
  the parameter vector. If there are no matches, and no query, the return
  value will be -3. No matches on a query returns -1.

  A single short match, or a single match of any length with a query, will
  result in a short help message

  If present, these values are set as follows:
    * matchCntp is set to the number of parameters that matched.
    * shortCntp is set to the number of matches that failed to meet the minimum
      match requirement.
    * queryCntp is set to the number of trailing `?' characters at the end
      of name.

  Return values:
    >0:	index of the single unique match for the name
    -1: query present
    -2: no query, one or more short matches
    -3: no query, no match
    -4: multiple full matches (indicates configuration error)

  The final three parameters (matchCnt, shortCnt, queryCnt) are optional and
  default to null. Use them if you want more detail on the match.
*/

int lookupParam(std::string name, CoinParamVec &paramVec,
                int *matchCntp, int *shortCntp, int *queryCntp)

{
  int retval = -3;

  // crude fix to stop no match message on clp when in cbc
  bool printMessage = true;
  if (queryCntp && *queryCntp == 999)
    printMessage = false;

  *matchCntp = 0;
  *shortCntp = 0;
  *queryCntp = 0;

  /*
  Is there anything here at all? 
*/
  if (name.length() == 0) {
    return (retval);
  }
  /*
  Scan the parameter name to see if it ends in one or more `?' characters. If
  so, take it as a request to return a list of parameters that match name up
  to the first `?'.  The strings '?' and '???' are considered to be valid
  parameter names (short and long help, respectively) and are handled as
  special cases: If the whole string is `?'s, one and three are commands as
  is, while 2 and 4 or more are queries about `?' or `???'.
*/
  int numQuery = 0;
  {
    int length = static_cast< int >(name.length());
    int i;
    for (i = length - 1; i >= 0 && name[i] == '?'; i--) {
      numQuery++;
    }
    if (numQuery == length) {
      switch (length) {
      case 1:
      case 3: {
        numQuery = 0;
        break;
      }
      case 2: {
        numQuery -= 1;
        break;
      }
      default: {
        numQuery -= 3;
        break;
      }
      }
    }
    name = name.substr(0, length - numQuery);
    if (queryCntp != 0) {
      *queryCntp = numQuery;
    }
  }
  /*
  See if we can match the parameter name. On return, matchNdx is set to the
  last match satisfying the minimal match criteria, or -1 if there's no
  match.  matchCnt is the number of matches satisfying the minimum match
  length, and shortCnt is possible matches that were short of the minimum
  match length,
*/
  int matchNdx = -1;
  int shortCnt = 0;
  int matchCnt = CoinParamUtils::matchParam(paramVec, name, matchNdx, shortCnt);
  /*
  Set up return values before we get into further processing.
*/
  if (matchCntp != 0) {
    *matchCntp = matchCnt;
  }
  if (shortCntp != 0) {
    *shortCntp = shortCnt;
  }
  if (numQuery > 0) {
    if (/* was numQuery! */  (matchCnt + shortCnt == 0) && printMessage){
        std::cout << "Sorry, no match found!" << std::endl;
     }
    retval = -1;
  } else {
    if (matchCnt + shortCnt == 0) {
      retval = -3;
    } else if (matchCnt > 1) {
      retval = -4;
    } else {
      retval = -2;
    }
  }
  /*
  No matches? Nothing more to be done here.
*/
  if (matchCnt + shortCnt == 0) {
    return (retval);
  }
  /*
  A unique match and no `?' in the name says we have our parameter. Return
  the result.
  allow if exact match (e.g. -netlib)
*/
  if (matchCnt == 1 /* && shortCnt == 0 */ && numQuery == 0) {
    assert(matchNdx >= 0 && matchNdx < static_cast< int >(paramVec.size()));
    return (matchNdx);
  }
  /*
  A single match? There are two possibilities:
    * The string specified is shorter than the match length requested by the
      parameter. (Useful for avoiding inadvertent execution of commands that
      the client might regret.)
    * The string specified contained a `?', in which case we print the help.
      The match may or may not be short.
*/
  if (matchCnt + shortCnt == 1) {
    CoinParamUtils::shortOrHelpOne(paramVec, matchNdx, name, numQuery);
    return (retval);
  }
  /*
  The final case: multiple matches. Most commonly this will be multiple short
  matches. If we have multiple matches satisfying the minimal length
  criteria, we have a configuration problem.  The other question is whether
  the user wanted help information. Two question marks gets short help.
*/
  if (matchCnt > 1) {
    std::cout
      << "Configuration error! `" << name
      << "' was fully matched " << matchCnt << " times!"
      << std::endl;
  }
  std::cout
    << "Multiple matches for `" << name << "'; possible completions:"
    << std::endl;
  CoinParamUtils::shortOrHelpMany(paramVec, name, numQuery);

  return (retval);
}

//#############################################################################
//#############################################################################

/*
  Utility function to scan a parameter vector for matches. Sets matchNdx to
  the index of the last parameter that meets the minimal match criteria (but
  note there should be at most one such parameter if the parameter vector is
  properly configured). Sets shortCnt to the number of short matches (should
  be zero in a properly configured vector if a minimal match is found).
  Returns the number of matches satisfying the minimal match requirement
  (should be 0 or 1 in a properly configured vector).

  The routine allows for the possibility of null entries in the parameter
  vector.

  In order to handle `?' and `???', there's nothing to it but to force a
  unique match if we match `?' exactly. (This is another quirk of clp/cbc
  parameter parsing, which we need to match for historical reasons.)
*/

int matchParam(const CoinParamVec &paramVec, std::string name,
               int &matchNdx, int &shortCnt)

{
  int vecLen = static_cast< int >(paramVec.size());
  int matchCnt = 0;

  matchNdx = -1;
  shortCnt = 0;

  for (int i = 0; i < vecLen; i++) {
    int match = paramVec[i]->matches(name);
    if (match == 1) {
      matchNdx = i;
      matchCnt++;
      if (name == "?") {
        matchCnt = 1;
        break;
      }
    } else {
      shortCnt += match >> 1;
    }
  }

  return (matchCnt);
}

/*
  Now a bunch of routines that are useful in the context of generating help
  messages.
*/

/*
  Simple formatting routine for long messages. Used to print long help for
  parameters. Lines are broken at the first white space after 65 characters,
  or when an explicit return (`\n') character is scanned. Leading spaces are
  suppressed.
*/

void printIt(const char *msg)

{
  int length = static_cast< int >(strlen(msg));
  char temp[101];
  int i;
  int n = 0;
  for (i = 0; i < length; i++) {
    if (msg[i] == '\n' || (n >= 65 && (msg[i] == ' ' || msg[i] == '\t'))) {
      temp[n] = '\0';
      std::cout << temp << std::endl;
      n = 0;
    } else if (n || msg[i] != ' ') {
      temp[n++] = msg[i];
    }
  }
  if (n > 0) {
    temp[n] = '\0';
    std::cout << temp << std::endl;
  }

  return;
}

//#############################################################################
//#############################################################################

/*
  Utility function for the case where a name matches a single parameter, but
  either it's short, or the user wanted help, or both.

  The routine allows for the possibility that there are null entries in the
  parameter vector, but matchNdx should point to a valid entry if it's >= 0.
*/

void shortOrHelpOne(CoinParamVec &paramVec,
  int matchNdx, std::string name, int numQuery)

{
  int i;
  int numParams = static_cast< int >(paramVec.size());
  int lclNdx = -1;
  /*
  For a short match, we need to look up the parameter again. This should find
  a short match, given the conditions where this routine is called. But be
  prepared to find a full match.
  
  If matchNdx >= 0, just use the index we're handed.
*/
  if (matchNdx < 0) {
    int match = 0;
    for (i = 0; i < numParams; i++) {
      int match = paramVec[i]->matches(name);
      if (match != 0) {
        lclNdx = i;
        break;
      }
    }

    assert(lclNdx >= 0);

    if (match == 1) {
      std::cout
        << "Match for '" << name << "': "
        << paramVec[matchNdx]->matchName() << ".";
    } else {
      std::cout
        << "Short match for '" << name << "'; possible completion: "
        << paramVec[lclNdx]->matchName() << ".";
    }
  } else {
    assert(matchNdx >= 0 && matchNdx < static_cast< int >(paramVec.size()));
    std::cout << "Match for `" << name << "': "
              << paramVec[matchNdx]->matchName();
    lclNdx = matchNdx;
  }
  /*
  Print some help, if there was a `?' in the name. `??' gets the long help.
*/
  if (numQuery > 0) {
    std::cout << std::endl;
    if (numQuery == 1) {
      std::cout << paramVec[lclNdx]->shortHelp();
    } else {
      std::cout << paramVec[lclNdx]->printLongHelp();
    }
  }
  std::cout << std::endl;

  return;
}

//#############################################################################
//#############################################################################

/*
  Utility function for the case where a name matches multiple parameters.
  Zero or one `?' gets just the matching names, while `??' gets short help
  with each match.

  The routine allows for the possibility that there are null entries in the
  parameter vector.
*/

void shortOrHelpMany(CoinParamVec &paramVec, std::string name, int numQuery)

{
  int numParams = static_cast< int >(paramVec.size());
  /*
  Scan the parameter list. For each match, print just the name, or the name
  and short help.
*/
  int lineLen = 0;
  bool printed = false;
  for (int i = 0; i < numParams; i++) {
    int match = paramVec[i]->matches(name);
    if (match > 0 &&
        paramVec[i]->getDisplayPriority() != CoinParam::displayPriorityNone) {
      std::string nme = paramVec[i]->matchName();
      int len = static_cast< int >(nme.length());
      if (numQuery >= 2) {
        std::cout << nme << " : " << paramVec[i]->shortHelp();
        std::cout << std::endl;
      } else {
        lineLen += 2 + len;
        if (lineLen > 80) {
          std::cout << std::endl;
          lineLen = 2 + len;
        }
        std::cout << "  " << nme;
        printed = true;
      }
    }
  }

  if (printed) {
    std::cout << std::endl;
  }

  return;
}

//#############################################################################
//#############################################################################

/*
  A generic help message that explains the basic operation of parameter
  parsing.
*/

void printGenericHelp()

{
  std::cout << std::endl;
  std::cout
    << "For command line arguments, keywords have a leading `-' or '--'; "
    << std::endl;
  std::cout
    << "-stdin or just - switches to stdin with a prompt."
    << std::endl;
  std::cout
    << "When prompted, one command per line, without the leading `-'."
    << std::endl;
  std::cout
    << "abcd value sets abcd to value."
    << std::endl;
  std::cout
    << "abcd without a value (where one is expected) gives the current value."
    << std::endl;
  std::cout
    << "abcd? gives a list of possible matches; if there's only one, a short"
    << std::endl;
  std::cout
    << "help message is printed."
    << std::endl;
  std::cout
    << "abcd?? prints the short help for all matches; if there's only one"
    << std::endl;
  std::cout
    << "match, a longer help message and current value are printed."
    << std::endl;

  return;
}

//#############################################################################
//#############################################################################

/*
  Utility function for various levels of `help' command. The entries between
  paramVec[firstParam] and paramVec[lastParam], inclusive, will be printed.
  If shortHelp is true, the short help message will be printed for each
  parameter. If longHelp is true, the long help message will be printed for
  each parameter.  If hidden is true, even parameters with display = false
  will be printed. Each line is prefaced with the specified prefix.

  The routine allows for the possibility that there are null entries in the
  parameter vector.
*/

void printHelp(CoinParamVec &paramVec, int firstParam, int lastParam,
  std::string prefix, bool shortHelp, bool longHelp, bool hidden)
{
  bool noHelp = !(shortHelp || longHelp);
  int i;
  int pfxLen = static_cast< int >(prefix.length());
  bool printed = false;
  CoinParam *param; 

  if (noHelp) {
    int lineLen = 0;
    for (i = firstParam; i <= lastParam; i++) {
      param = paramVec[i]; 
      if (param == 0){
         continue;
      }
      if (param->getDisplayPriority() || hidden) {
        std::string nme = param->matchName();
        int len = static_cast< int >(nme.length());
        if (!printed) {
          std::cout << std::endl
                    << prefix;
          lineLen += pfxLen;
          printed = true;
        }
        lineLen += 2 + len;
        if (lineLen > 80) {
          std::cout << std::endl
                    << prefix;
          lineLen = pfxLen + 2 + len;
        }
        std::cout << "  " << nme;
      }
    }
    if (printed) {
      std::cout << std::endl;
    }
  } else if (shortHelp) {
    for (i = firstParam; i <= lastParam; i++) {
      param = paramVec[i]; 
      if (param == 0){
         continue;
      }
      if (param->getDisplayPriority() || hidden) {
        std::cout << std::endl
                  << prefix;
        std::cout << param->matchName();
        std::cout << ": ";
        std::cout << param->shortHelp();
      }
    }
    std::cout << std::endl;
  } else if (longHelp) {
    for (i = firstParam; i <= lastParam; i++) {
      param = paramVec[i]; 
      if (param == 0){
         continue;
      }
      if (param->getDisplayPriority() || hidden) {
        std::cout << std::endl
                  << prefix;
        std::cout << "Command: " << param->matchName();
        std::cout << std::endl
                  << prefix;
        std::cout << "---- description" << std::endl;
        printIt(param->longHelp().c_str());
        std::cout << prefix << "----" << std::endl;
      }
    }
  }

  std::cout << std::endl;

  return;
}

//#############################################################################
//#############################################################################

/* Take in file name and directory name, determine whether file name 
   is absolute and if not, prepend with directory name. */

void processFile(std::string &fileName, std::string dirName,
                 bool *fileExists){

   if (fileName == "--" || fileName == "stdin"){
      // This is for the case of stdin
      fileName = "-";
      if (fileExists != NULL){
         *fileExists = true;
      }
   } 
   if (fileName == "stdin_lp"){
      // This is for the case of stdin
      fileName = "-lp";
      if (fileExists != NULL){
         *fileExists = true;
      }
   } 

   if (fileName == "$" || fileName == ""){
      fileName = "";
      return;
   }
   if (fileName[0] != '/' && fileName[0] != '\\' &&
       !(fileName.length() >= 2 && fileName[1] == ':' &&
         ((fileName[0] >= 'a' && fileName[0] <= 'z') ||
          (fileName[0] >= 'A' && fileName[0] <= 'Z')))) {
      fileName = dirName + fileName;
   }
   if (fileExists != NULL){
      *fileExists = fileCoinReadable(fileName);
   }
}

//#############################################################################
//#############################################################################

void formInputQueue(std::deque<std::string> &inputQueue,
                    std::string commandName,
                    int argc, char **argv)
{
   std::string::size_type found;
   for (int i = 1; i < argc; i++){
      std::string tmp(argv[i]);
      if (tmp == commandName) {
         // For some reason, the command can sometimes be listed more than once
         // in argv
         continue;
      }
      found = tmp.find('=');
      if (found != std::string::npos) {
         inputQueue.push_back(tmp[0] == '-' ? tmp.substr(0, found) :
                              '-' + tmp.substr(0, found));
         inputQueue.push_back(tmp.substr(found + 1));
      } else {
         inputQueue.push_back(tmp);
      }
   }
}

} // end namespace CoinParamUtils

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
