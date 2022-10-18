/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _TZ_SOURCE_FACTORY_H_
#define _TZ_SOURCE_FACTORY_H_

namespace Loki {

class LokiInputParser;
class TZSourceInterface;

/**
 * Factory for creating twilight zone source objects.  Only contains static
 * create method.
 */
class TZSourceFactory
{
public:
   /**
    * @brief Creates the requested TZSourceInterface.
    *
    * @param[in] a_pp Where the twilight zone source's input parameters will be
    *                 found.
    *
    * @return The requested TZSourceInterface.
    */
   static
   TZSourceInterface*
   create(
      LokiInputParser& a_pp);

private:
   // Unimplemented default constructor.
   TZSourceFactory();

   // Unimplemented destructor.
   ~TZSourceFactory();

   // Unimplemented copy constructor.
   TZSourceFactory(const TZSourceFactory& other);

   // Unimplemented assignment operator.
   TZSourceFactory&
   operator = (
      const TZSourceFactory& rhs);
};

} // end namespace Loki

#endif
