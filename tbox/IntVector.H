/*************************************************************************
 *
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * Written by Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute,
 * Amos Eaton 301, 110 8th St., Troy, NY 12180); Jeffrey Hittinger
 * hittinger1@llnl.gov, William Arrighi arrighi2@llnl.gov, Richard Berger
 * berger5@llnl.gov, Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808,
 * Livermore, CA 94551); Stephan Brunner stephan.brunner@epfl.ch (Ecole
 * Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13,
 * CH-1015 Lausanne, Switzerland).
 * CODE-744849
 *
 * All rights reserved.
 *
 * This file is part of Loki.  For details, see.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 ************************************************************************/
#ifndef included_tbox_IntVector
#define included_tbox_IntVector

#include "Dimension.H"
#include "HDF5_DataBase.h"
#include "MathUtilities.H"

#include <vector>

namespace Loki {
namespace tbox {
   
/**
 * Class IntVector implements a simple N-dimensional integer vector.  This
 * class is the base class for most of the simple indexing classes.
 *
 */
class IntVector
{
   public:

      /**
       * Creates an uninitialized vector.
       */
      IntVector(const Dimension& dim);

      /**
       * Construct an integer vector with all components equal to the argument.
       */
      IntVector(const Dimension& dim, const int i);

      /**
       * Construct a n-dimensional integer vector with the value with
       * values provided by the array.
       *
       * Dimension inferred from array size.
       */
      IntVector(const std::vector<int>& a);

      /**
       * Construct a n-dimensional integer vector with the value with
       * values provided by the array.
       *
       */
      IntVector( const Dimension& dim, const int array[]);

      /**
       * Construct an integer vector equal to the argument.
       */
      IntVector(const IntVector& rhs);

      /**
       * The assignment operator sets the integer vector equal to the argument.
       *
       * An assignment to an uninitialized IntVector is allowed but assigning
       * from an uninitialized IntVector will result in an assert.
       */
      IntVector& operator=(const IntVector& rhs)
      {
         /*
          * Allow assignment of to an uninitialized but do not allow
          * assignment from an uninitialized.
          */
         if (d_dim.isValid()) {
            TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
         } else {
            TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());
            d_dim = rhs.getDim();
         }

         for (int i = 0; i < d_dim; i++)
            d_vector[i] = rhs.d_vector[i];

//#ifdef DEBUG_INITIALIZE_UNDEFINED
         for (int i = d_dim; i < Dimension::MAXIMUM_DIMENSION; i++) {
            d_vector[i] = MathUtilities<int>::getMin();
         }
//#endif

         return *this;
      }

      /**
       * The integer vector destructor does nothing interesting.
       */
      virtual ~IntVector();

      /**
       * Return the specified component of the vector.  No bounds checking.
       */
      int& operator[](const int i)
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
         TBOX_ASSERT(i >= 0 && i < d_dim);
         return d_vector[i];
      }

      /**
       * Return the specified component of the vector as a const integer.
       * No bounds checking.
       */
      const int& operator[](const int i) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
         TBOX_ASSERT(i >= 0 && i < d_dim);
         return d_vector[i];
      }

      /**
       * Return the specified component of the vector.  No bounds checking.
       */
      int& operator()(const int i)
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
         TBOX_ASSERT(i >= 0 && i < d_dim);
         return d_vector[i];
      }

      /**
       * Return the specified component of the vector as a const integer.
       * No bounds checking.
       */
      const int& operator()(const int i) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
         TBOX_ASSERT(i >= 0 && i < d_dim);
         return d_vector[i];
      }

      /**
       * Plus-equals operator for two integer vectors.
       */
      IntVector& operator+=(const IntVector& rhs)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] += rhs.d_vector[i];

         return *this;
      }

      /**
       * Plus operator for two integer vectors.
       */
      IntVector operator+(const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         IntVector tmp(*this);
         tmp += rhs;
         return tmp;
      }

      /**
       * Plus-equals operator for an integer vector and an integer.
       */
      IntVector& operator+=(const int rhs)
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] += rhs;

         return *this;
      }

      /**
       * Plus operator for an integer vector and an integer.
       */
      IntVector operator+(const int rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         IntVector tmp = *this;
         tmp += rhs;
         return tmp;
      }

      /**
       * Minus-equals operator for two integer vectors.
       */
      IntVector& operator-=(const IntVector& rhs)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] -= rhs.d_vector[i];

         return *this;
      }

      /**
       * Minus operator for two integer vectors.
       */
      IntVector operator-(const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         IntVector tmp = *this;
         tmp -= rhs;
         return tmp;
      }

      /**
       * Minus-equals operator for an integer vector and an integer.
       */
      IntVector& operator-=(const int rhs)
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] -= rhs;

         return *this;
      }

      /**
       * Minus operator for an integer vector and an integer.
       */
      IntVector operator-(const int rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         IntVector tmp = *this;
         tmp -= rhs;
         return tmp;
      }

      /**
       * Times-equals operator for two integer vectors.
       */
      IntVector&
      operator *= (
         const IntVector& rhs)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] *= rhs.d_vector[i];

         return *this;
      }

      /**
       * Times operator for two integer vectors.
       */
      IntVector
      operator * (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         IntVector tmp = *this;
         tmp *= rhs;
         return tmp;
      }

      /**
       * Times-equals operator for an integer vector and an integer.
       */
      IntVector&
      operator *= (
         const int rhs)
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] *= rhs;

         return *this;
      }

      /**
       * Times operator for an integer vector and an integer.
       */
      IntVector
      operator * (
         const int rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         IntVector tmp = *this;
         tmp *= rhs;
         return tmp;
      }

      /**
       * Assign-quotient operator for two integer vectors.
       */
      IntVector&
      operator /= (
         const IntVector& rhs)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] /= rhs.d_vector[i];

         return *this;
      }

      /**
       * Component-wise ceiling quotient (integer divide with rounding up).
       */
      void
      ceiling(
         const IntVector& denominator)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, denominator);

         for (int i = 0; i < getDim(); i++) {
            /*
             * This is the formula for integer divide, rounding away from
             * zero.  It is meant as an extension of the ceiling quotient of
             * 2 positive integers.
             *
             * The ceiling is the integer divide plus 0, -1 or 1 representing
             * the results of rounding.
             * - Add zero if there's no remainder to round.
             * - Round remainder to 1 if numerator and denominator has same
             *   sign.
             * - Round remainder to -1 if numerator and denominator has
             *   opposite sign.
             */
            d_vector[i] = (d_vector[i] / denominator[i])
               + ((d_vector[i] % denominator[i]) ?
                  ((d_vector[i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
         }
      }

      /**
       * Component-wise ceiling quotient (integer divide with rounding up).
       */
      static IntVector
      ceiling(
         const IntVector& numerator,
         const IntVector& denominator)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(numerator, denominator);

         IntVector rval(numerator.getDim());
         for (int i = 0; i < numerator.getDim(); i++) {
            /*
             * This is the formula for integer divide, rounding away from
             * zero.  It is meant as an extension of the ceiling quotient of
             * 2 positive integers.
             *
             * The ceiling is the integer divide plus 0, -1 or 1 representing
             * the results of rounding.
             * - Add zero if there's no remainder to round.
             * - Round remainder to 1 if numerator and denominator has same
             *   sign.
             * - Round remainder to -1 if numerator and denominator has
             *   opposite sign.
             */
            rval[i] = (numerator[i] / denominator[i])
               + ((numerator[i] % denominator[i]) ?
                  ((numerator[i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
         }
         return rval;
      }

      /**
       * Quotient operator for two integer vectors.
       */
      IntVector
      operator / (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         IntVector tmp = *this;
         tmp /= rhs;
         return tmp;
      }

      /**
       * Assign-quotient operator for an integer vector and an integer.
       */
      IntVector&
      operator /= (
         const int rhs)
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         for (int i = 0; i < d_dim; i++)
            d_vector[i] /= rhs;

         return *this;
      }

      /**
       * Quotient operator for an integer vector and an integer.
       */
      IntVector
      operator / (
         const int rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         IntVector tmp = *this;
         tmp /= rhs;
         return tmp;
      }

      /**
       * Unary minus to negate an integer vector.
       */
      IntVector
      operator - () const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         IntVector tmp(d_dim);
         for (int i = 0; i < d_dim; i++)
            tmp.d_vector[i] = -d_vector[i];

         return tmp;
      }

      /**
       * Returns true if all components are equal to a given integer.
       */
      bool
      operator == (
         int rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = d_vector[i] == rhs;

         return result;
      }

      /**
       * Returns true if some components are not equal to a given integer.
       */
      bool
      operator != (
         int rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = d_vector[i] != rhs;

         return result;
      }

      /**
       * Returns true if two vector objects are equal.  All components
       * must be the same for equality.
       */
      bool
      operator == (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = result && (d_vector[i] == rhs.d_vector[i]);

         return result;
      }

      /**
       * Returns true if two vector objects are not equal.  Any of
       * the components may be different for inequality.
       */
      bool
      operator != (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         return !(*this == rhs);
      }

      /**
       * Returns true if each integer in vector is less than
       * corresponding integer in comparison vector.
       */
      bool
      operator < (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = result && (d_vector[i] < rhs.d_vector[i]);

         return result;
      }

      /**
       * Returns true if each integer in vector is less or equal to
       * corresponding integer in comparison vector.
       */
      bool
      operator <= (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = result && (d_vector[i] <= rhs.d_vector[i]);

         return result;

      }

      /**
       * Returns true if each integer in vector is greater than
       * corresponding integer in comparison vector.
       */
      bool
      operator > (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = result && (d_vector[i] > rhs.d_vector[i]);

         return result;
      }

      /**
       * Returns true if each integer in vector is greater or equal to
       * corresponding integer in comparison vector.
       */
      bool
      operator >= (
         const IntVector& rhs) const
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         bool result = true;

         for (int i = 0; result && (i < d_dim); i++)
            result = result && (d_vector[i] >= rhs.d_vector[i]);

         return result;
      }

      /**
       * Return the component-wise minimum of two integer vector objects.
       */
      void
      min(
         const IntVector& rhs)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         for (int i = 0; i < d_dim; i++)
            if (rhs.d_vector[i] < d_vector[i])
               d_vector[i] = rhs.d_vector[i];

      }

      /**
       * Return the minimum entry in an integer vector.
       */
      int
      min() const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         int min = d_vector[0];

         for (int i = 1; i < d_dim; i++)
            if (d_vector[i] < min)
               min = d_vector[i];

         return min;
      }

      /**
       * Return the component-wise maximum of two integer vector objects.
       */
      void
      max(
         const IntVector& rhs)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

         for (int i = 0; i < d_dim; i++)
            if (rhs.d_vector[i] > d_vector[i])
               d_vector[i] = rhs.d_vector[i];
      }

      /**
       * Return the maximum entry in an integer vector.
       */
      int
      max() const
      {
         TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

         int max = d_vector[0];

         for (int i = 1; i < d_dim; i++)
            if (d_vector[i] > max)
               max = d_vector[i];

         return max;
      }

      /**
       * Utility function to take the minimum of two integer vector objects.
       */
      static IntVector
      min(
         const IntVector& a,
         const IntVector& b)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(a, b);

         IntVector tmp = a;
         tmp.min(b);
         return tmp;
      }

      /**
       * Utility function to take the maximum of two integer vector objects.
       */
      static IntVector
      max(
         const IntVector& a,
         const IntVector& b)
      {
         TBOX_DIM_ASSERT_CHECK_ARGS2(a, b);

         IntVector tmp = a;
         tmp.max(b);
         return tmp;
      }

      /**
       * Return the product of the entries in the integer vector.
       */
      int
      getProduct() const
      {
         int prod = 1;
         for (int i = 0; i < d_dim; i++) {
            prod *= d_vector[i];
         }
         return prod;
      }

      /**
       * Store the object state to the specified database
       * with the provided name.
       *
       */
      void
      putToDatabase(
         HDF_DataBase& database,
         const std::string& name) const;

      /**
       * Restores the object state from the specified database
       * with the provided name.
       *
       */
      void
      getFromDatabase(
         const HDF_DataBase& database,
         const std::string& name);

      /**
       * Return the dimension of this object.
       */
      const Dimension&
      getDim() const
      {
         return d_dim;
      }

      /*!
       * @brief Return an IntVector of zeros of the specified dimension.
       *
       * Can be used to avoid object creation overheads.
       */
      static const IntVector&
      Zero(
         const Dimension& dim)
      {
         return *(s_zeros[dim - 1]);
      }

      /*!
       * @brief Return an IntVector of ones of the specified dimension.
       *
       * Can be used to avoid object creation overheads.
       */
      static const IntVector&
      Unit(
         const Dimension& dim)
      {
         return *(s_ones[dim - 1]);
      }

      /**
       * Read an integer vector from an input stream.  The format for
       * the input is (i0,...,in) for an n-dimensional vector.
       */
      friend std::istream&
      operator >> (
         std::istream& s,
         IntVector& rhs);

      /**
       * Write an integer vector into an output stream.  The format for
       * the output is (i0,...,in) for an n-dimensional vector.
       */
      friend std::ostream&
      operator << (
         std::ostream& s,
         const IntVector& rhs);

      static
      void
      shutdown()
      {
         for (int i = 0; i <Dimension::MAXIMUM_DIMENSION; ++i) {
            if (s_zeros[i]) {
               delete s_zeros[i];
            }
            if (s_ones[i]) {
               delete s_ones[i];
            }
         }
      }

      friend class std::vector<IntVector>;

   protected:
      /**
       * Default ctor for IntVector is protected to disallow normal use.
       * This is needed by the poorly designed STL container library.
       *
       *
       */
      IntVector();

   private:
      /*!
       * @brief Initialize static objects and register shutdown routine.
       */
      static bool
      initializeStatics();

      /*!
       * @brief Method registered with ShutdownRegister to cleanup statics.
       */
      static void
      freeStatics();

      Dimension d_dim;

      int d_vector[Dimension::MAXIMUM_DIMENSION];

      static bool s_registered_shutdown;

      static IntVector* s_zeros[Dimension::MAXIMUM_DIMENSION];
      static IntVector* s_ones[Dimension::MAXIMUM_DIMENSION];
   };
}
}


#define INTVECTOR1D_TO_FORT(__iv) \
   (__iv)[0]

#define INTVECTOR2D_TO_FORT(__iv) \
   (__iv)[0], (__iv)[1]

#define INTVECTOR3D_TO_FORT(__iv) \
   (__iv)[0], (__iv)[1], (__iv)[2]

#define INTVECTOR4D_TO_FORT(__iv) \
   (__iv)[0], (__iv)[1], (__iv)[2], (__iv)[3]

#endif
