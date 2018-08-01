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
#ifndef included_tbox_Dimension
#define included_tbox_Dimension

#include "Utilities.H"

#include <limits>

namespace Loki {
namespace tbox {

class IntVector;

/**
 * Class Dimension is used to represent the dimension of a Loki object.
 *
 * The maximum dimension is set at compile time using a flag to the
 * configure script.  This is used to allocate arrays in some lower
 * level classes such as IntVector.  If dynamic memory allocation is
 * used the performance impact is significant; a maximum dimension
 * allows for stack based memory allocation in performance critical
 * classes at the expense of wasting storage for objects with
 * dimension less than the maximum dimension.
 *
 * A class is used rather than a simple short or integer to provide
 * enhanced type safety.
 *
 */
class Dimension
{
public:
   /**
    * Constructor for Dimension, object is built using the specified dimension
    *
    * Note that the constructor is "explicit" thus making automatic
    * type conversions from integers impossible.  This is intentionally to avoid
    * unintended conversions.
    *
    * When dimensional assertion checking is active an assert is
    * thrown when dim < 1 or dim > MAX_DIMENSION value specified when
    * the library is configured (defaults to 3).  dim also cannot be
    * the INVALID_DIMENSION (the largest unsigned short value).
    *
    */
   explicit Dimension(
      const unsigned short& dim);

   /**
    * Construct a dimension equal to the argument.
    */
   Dimension(
      const Dimension& dimension);

   /**
    * Returns true if Dimension is valid.
    *
    * A valid Dimension != 0; != INVALID_DIMENSION, and <= MAXIMUM_DIMENSION.
    *
    */
   bool
   isValid() const
   {
      /*
       * Note we don't check for != INVALID_DIMENSION because it is
       * maximum value for a unsigned short so this condition will
       * be caught by the < MAXIMUM_DIMENSION clause.
       *
       * If for some reason INVALID_DIMENSION implementation is changed
       * may need to be fixed.
       */
      return (d_dim != 0) && (d_dim <= Dimension::MAXIMUM_DIMENSION);
   }

   /**
    * Returns true if Dimension is initialized (not set to INVALID_DIMENSION).
    *
    */
   bool
   isInitialized() const
   {
      return d_dim != INVALID_DIMENSION_VALUE;
   }

   /**
    * Equality operator.
    */
   bool
   operator == (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim == rhs.d_dim;
   }

   /**
    * Inequality operator.
    */
   bool
   operator != (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim != rhs.d_dim;
   }

   /**
    * Greater than operator.
    */
   bool
   operator > (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim > rhs.d_dim;
   }

   /**
    * Greater than or equal operator.
    */
   bool
   operator >= (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim >= rhs.d_dim;
   }

   /**
    * Less than operator.
    */
   bool
   operator < (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim < rhs.d_dim;
   }

   /**
    * Less than or equal operator.
    */
   bool
   operator <= (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim <= rhs.d_dim;
   }

   /**
    * Returns the dimension of the Dimension.
    *
    * The method is provided to allow sizing of arrays based on the
    * dimension and for iteration.  In general this should not be
    * used for comparisons, the Dimension comparison operations are
    * better suited for that purpose.
    */
   operator unsigned short () const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      return d_dim;
   }

   /**
    * Returns the maximum dimension for the currently compiled library.
    *
    * When the Loki library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This method is typically used
    * to allocate arrays.
    *
    */
   unsigned short
   getMaxDim() const
   {
      return MAXIMUM_DIMENSION;
   }

   /**
    * The maximum dimension for the currently compiled library.
    *
    * When the Loki library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This value is typically used to
    * allocate arrays for performance critical classes.  E.G.
    *
    *  double array[Loki::tbox::Dimension::MAXIMUM_DIMENSION];
    *
    * The value must be >= 1 and < numeric_limits<unsigned short>::max()
    *
    */
   static const unsigned short MAXIMUM_DIMENSION = 6;

   /**
    * An invalid dimension value as a Dimension object.
    */
   static const Dimension INVALID_DIMENSION;

   /**
    * An invalid dimension value as an unsigned short.
    *
    * Currently this value is numeric_limits<unsigned short>::max() but
    * use this symbol as it is more readable.
    *
    */
   static unsigned short INVALID_DIMENSION_VALUE;

   /*
    * Classes that are friends of dimension in order to access the
    * private ctor which builds invalid dimensions.
    *
    * This is obviously not a very good design but so far
    * a better solution has been elusive.   Allowing
    * any code to create invalid dimensions seemed too
    * error prone.
    */
   friend class IntVector;

   /**
    * Output operator for debugging and error messages.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Dimension& rhs);

private:
   /**
    * @brief Create an invalid dimension object.
    *
    * This ctor is private to prevent a default constructor call.
    * Currently Dimension objects must always created with a dimension
    * specified for normal code.  Several special classes are allowed
    * and are declared to be friends to access this ctor.
    */
   Dimension();

   /**
    * Assignment operator is private to prevent dimensions
    * from being assigned.  This was done to improve type
    * safety.
    */
   Dimension&
   operator = (
      const Dimension& rhs)
   {
      d_dim = rhs.d_dim;
      return *this;
   }

   unsigned short d_dim;
};

}
}

#endif
