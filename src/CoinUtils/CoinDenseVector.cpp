// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Resized.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <cassert>
#include "CoinDenseVector.hpp"
#include "CoinHelperFunctions.hpp"

//#############################################################################

template < typename T >
void CoinDenseVector< T >::clear()
{
  elements_.clear();
}

//#############################################################################

template < typename T >
CoinDenseVector< T > &
CoinDenseVector< T >::operator=(const CoinDenseVector< T > &rhs)
{
  if (this != &rhs) {
    setVector(rhs.getNumElements(), rhs.getElements());
  }
  return *this;
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::setVector(int size, const T *elems)
{
  elements_ = std::vector<T>(elems, elems + size);
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::setConstant(int size, T value)
{
  resize(size);
  for (int i = 0; i < size; i++)
    elements_[i] = value;
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::resize(size_t newsize, T value)
{
  if (newsize != elements_.size()) {
    assert(newsize > 0);
    elements_.resize(newsize, value);
  }
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::setElement(size_t index, T element)
{
  assert(index < elements_.size());
  elements_[index] = element;
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::append(const CoinDenseVector< T > &caboose)
{
  const int cs = caboose.getNumElements();
  const T *celem = caboose.getElements();
  for (int i = 0; i< cs; ++i)
    elements_.push_back(celem[i]);
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::operator+=(T value)
{
  for (size_t i = 0; i < elements_.size(); i++)
    elements_[i] += value;
}

//-----------------------------------------------------------------------------

template < typename T >
void CoinDenseVector< T >::operator-=(T value)
{
  for (size_t i = 0; i < elements_.size(); i++)
    elements_[i] -= value;
}

//-----------------------------------------------------------------------------

template < typename T >
void CoinDenseVector< T >::operator*=(T value)
{
  for (size_t i = 0; i < elements_.size(); i++)
    elements_[i] *= value;
}

//-----------------------------------------------------------------------------

template < typename T >
void CoinDenseVector< T >::operator/=(T value)
{
  for (size_t i = 0; i < elements_.size(); i++)
    elements_[i] /= value;
}

//#############################################################################

template < typename T >
CoinDenseVector< T >::CoinDenseVector()
  : elements_()
{
}

//#############################################################################

template < typename T >
CoinDenseVector< T >::CoinDenseVector(size_t size, const T *elems)
  : elements_(elems, elems + size)
{
}

//-----------------------------------------------------------------------------

template < typename T >
CoinDenseVector< T >::CoinDenseVector(size_t size, T value)
  : elements_(size, value)
{
}

//-----------------------------------------------------------------------------

template < typename T >
CoinDenseVector< T >::CoinDenseVector(const CoinDenseVector< T > &rhs)
  : elements_()
{
  setVector(rhs.getNumElements(), rhs.getElements());
}

//-----------------------------------------------------------------------------

template < typename T >
CoinDenseVector< T >::~CoinDenseVector()
{
}

//#############################################################################

template < typename T >
void CoinDenseVector< T >::gutsOfSetVector(size_t size, const T *elems)
{
  if (size != 0) {
    elements_= std::vector<T>(elems, elems + size);
  }
}

//-----------------------------------------------------------------------------

template < typename T >
void CoinDenseVector< T >::gutsOfSetConstant(size_t size, T value)
{
  if (size != 0) {
    elements_ = std::vector<T>(size, value);
  }
}

//#############################################################################
/** Access the i'th element of the dense vector.  */
template < typename T >
T &
  CoinDenseVector< T >::operator[](size_t index)
{
  assert(index < elements_.size());
  return elements_[index];
}

/** Access the i'th element of the dense vector.  */
template < typename T >
const T &
  CoinDenseVector< T >::operator[](size_t index) const
{
  assert(index < elements_.size());
  return elements_[index];
}
//#############################################################################

// template class CoinDenseVector<int>; This works but causes warning messages
template class COINUTILSLIB_EXPORT CoinDenseVector<float>;
template class COINUTILSLIB_EXPORT CoinDenseVector<double>;

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
