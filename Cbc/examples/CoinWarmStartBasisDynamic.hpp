/*! \legal
  Copyright (C) 2006, International Business Machines Corporation
  and others.  All Rights Reserved.
*/

/*! \file CoinWarmStartBasisDynamic.hpp
  \brief Declaration of the dynamic simplex (basis-oriented) warm start
  class.
*/

#ifndef CoinWarmStartBasisDynamic_H
#define CoinWarmStartBasisDynamic_H

#include "CoinWarmStartBasis.hpp"

//#############################################################################

/*! \class CoinWarmStartBasisDynamic
  \brief The dynamic COIN simplex (basis-oriented) warm start class
  
  CoinWarmStartBasisDynamic provides for a warm start object which contains the
  status of each variable (structural and artificial) and a list of variables in problem.
  
*/

class CoinWarmStartBasisDynamic : public CoinWarmStartBasis {
  
public:
  
  /*! \name Methods to get and set basis information.
    
  By "common" I mean in all problems
  
  Much of this is in CoinWarmStartBasis
  */
  //@{
  /// Return the number of common variables
  inline int getNumCommonVariables() const { return numberCommonVariables_; }
  /// Set the number of common variables
  inline void setNumCommonVariables(int number) { numberCommonVariables_ = number; }
  
  /// Return the number of dynamic variables
  inline int getNumDynamicVariables() const { return numberDynamicVariables_; }
  
  /// Return the number of common rows (probably not used)
  inline int getNumCommonRows() const { return numberCommonRows_; }
  
  /** Return the list of identifiers for the dynamic structural variables
      
  There will be numDynamicVariables entries
  The encoding is not defined
  */
  inline const int * getDynamicVariables() const { return dynamicVariables_; }
  
  /**  Save list of dynamic variables */
  
  void setDynamicVariables(int numberDynamicVariables, const int * dynamicVariables);
  
  //@}
  
  /** \brief Delete a set of columns from the basis
      
  \warning
  The resulting basis is guaranteed valid only if all deleted variables
  are nonbasic.
  
  Removal of a basic variable implies that some nonbasic variable must be
  made basic. This correction is left to the client.
  */
  
  virtual void deleteColumns(int number, const int * which);
  
  //@}
  
  /*! \name Constructors, destructors, and related functions */
  
  //@{
  
  /** Default constructor
      
  Creates a warm start object representing an empty basis
  (0 rows, 0 columns).
  */
  CoinWarmStartBasisDynamic();
  
  /** Constructs a warm start object with the specified status vectors.
      
  The parameters are copied.
  Consider assignBasisStatus(int,int,char*&,char*&) if the object should
  assume ownership.
  
  \sa CoinWarmStartBasis::Status for a description of the packing used in
  the status arrays.
  */
  CoinWarmStartBasisDynamic(int ns, int na, const char* sStat, const char* aStat,
                            int numberCommon, int numberDynamicVariables, const int * dynamicVariables) ;
  
  /** Copy constructor */
  CoinWarmStartBasisDynamic(const CoinWarmStartBasisDynamic& rhs) ;
  
  /** `Virtual constructor' */
  virtual CoinWarmStart *clone() const
  {
    return new CoinWarmStartBasisDynamic(*this);
  }
  
  /** Destructor */
  virtual ~CoinWarmStartBasisDynamic();
  
  /** Assignment */
  
  virtual CoinWarmStartBasisDynamic& operator=(const CoinWarmStartBasisDynamic& rhs) ;
  
  /** Assign the status vectors to be the warm start information.
      
  In this method the CoinWarmStartBasis object assumes ownership of the
  pointers and upon return the argument pointers will be NULL.
  If copying is desirable, use the
  \link CoinWarmStartBasis(int,int,const char*,const char*)
  array constructor \endlink
  or the
  \link operator=(const CoinWarmStartBasis&)
  assignment operator \endlink.
  
  \note
  The pointers passed to this method will be
  freed using delete[], so they must be created using new[].
  */
  virtual void assignBasisStatus(int ns, int na, char*& sStat, char*& aStat,
                                 int numberCommon, int numberDynamicVariables, int *& dynamicVariables) ;
  //@}
  
  /*! \name Miscellaneous methods */
  //@{
  
  /// Prints in readable format (for debug)
  virtual void print() const;
  
  /*! \brief Generate a `diff' that can convert the warm start basis passed as
	     a parameter to the warm start basis specified by \c this.

    The capabilities are limited: the basis passed as a parameter can be no
    larger than the basis pointed to by \c this.
  */

  virtual CoinWarmStartDiff*
  generateDiff (const CoinWarmStart *const oldCWS) const ;
  //@}
  
private:
  /** \name Private data members
      
  \sa CoinWarmStartBasis::Status for a description of the packing used in
  the status arrays.
  */
  //@{
  /// The number of common variables
  int numberCommonVariables_;
  /// The number of common rows
  int numberCommonRows_;
  /// The number of dynamic variables
  int numberDynamicVariables_;
  /** The dynamic variables. */
  int * dynamicVariables_;
  //@}
};


/*! \class CoinWarmStartBasisDiffDynamic
    \brief A `diff' between two CoinWarmStartBasisDynamic objects

    See CoinWarmStartBasis(Diff) for most ideas
    This version just saves new list
*/

class CoinWarmStartBasisDiffDynamic : public CoinWarmStartBasisDiff
{ 
public:

  /*! \brief `Virtual constructor' */
  virtual CoinWarmStartDiff *clone() const
  { CoinWarmStartBasisDiffDynamic *cwsbd =  new CoinWarmStartBasisDiffDynamic(*this) ;
    return (dynamic_cast<CoinWarmStartDiff *>(cwsbd)) ; }

  /*! \brief Assignment */
  virtual
    CoinWarmStartBasisDiffDynamic &operator= (const CoinWarmStartBasisDiffDynamic &rhs) ;

  /*! \brief Destructor */
  virtual ~CoinWarmStartBasisDiffDynamic()
  { delete[] dynamic_; }

  protected:

  /*! \brief Default constructor
  
    This is protected (rather than private) so that derived classes can
    see it when they make <i>their</i> default constructor protected or
    private.
  */
  CoinWarmStartBasisDiffDynamic () : CoinWarmStartBasisDiff(),numberDynamic_(0),
                                     dynamic_(NULL) { } ;

  /*! \brief Copy constructor
  
    For convenience when copying objects containing CoinWarmStartBasisDiff
    objects. But consider whether you should be using #clone() to retain
    polymorphism.

    This is protected (rather than private) so that derived classes can
    see it when the make <i>their</i> copy constructor protected or
    private.
  */
  CoinWarmStartBasisDiffDynamic (const CoinWarmStartBasisDiffDynamic &rhs) ;
  private:

  /*! \brief Standard constructor */
  CoinWarmStartBasisDiffDynamic (int sze, const unsigned int *const diffNdxs,
                                 const unsigned int *const diffVals,
                                 int numberDynamic, const int * dynamic) ;
  friend CoinWarmStartDiff *
    CoinWarmStartBasisDynamic::generateDiff(const CoinWarmStart *const oldCWS) const ;
  /* friend void
     CoinWarmStartBasis::applyDiff(const CoinWarmStartDiffDynamic *const diff) ;*/

  
  /*! \brief Number of entries (and allocated capacity), in units of \c int. */
  int sze_ ;

  /*! \brief Number of entries  */

  int numberDynamic_ ;

  /*! \brief Array of column indices */

  int *dynamic_ ;
} ;

#endif
