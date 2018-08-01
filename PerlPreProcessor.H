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
#ifndef _PERL_PRE_PROCESSOR_H_
#define _PERL_PRE_PROCESSOR_H_

#include "OvertureParser.h"

#include <string>

namespace Loki {

/**
 * Pre-processes Loki input files' perl syntax producing files with all perl
 * variables evaluated and suitable for parsing by Loki.  This scheme is
 * shamelessly stolen from Overture so a Loki input file has the same format
 * as an Overture input file.  The contents are Loki specific but the way that
 * embedded Perl is used is identical to Overture.  Hence, "all" we do here is
 * to use the OvertureParser to parse each logical line and, if appropriate, to
 * write that processed line to the output file.
 */
class PerlPreProcessor
{
public:
   /**
    * @brief Constructor.
    *
    * @param[in] a_argc Command line arg count.
    * @param[in] a_argv Command line args
    */
   PerlPreProcessor(
      int a_argc,
      char** a_argv);

   /**
    * @brief Destructor.
    */
   ~PerlPreProcessor();

   /**
    * @brief Evaluate all expressions in input file and produce resulting
    * output file.
    *
    * @param[in] a_input_file_name Name of unprocessed input file.
    *
    * @return Name of processed output file.
    */
   std::string
   process(
      const std::string& a_input_file_name);

private:
   // Unimplemented default constrcutor.
   PerlPreProcessor();

   // Unimplemented copy constructor.
   PerlPreProcessor(
      const PerlPreProcessor& other);

   // Unimplemented assignment operator.
   PerlPreProcessor&
   operator = (
      const PerlPreProcessor& rhs);

   // This function takes the result of the OvertureParser and determines if
   // the output line from the parser should be output.  The OvertureParser
   // return the following values:
   // 0 if the parsed line was unchanged (most likely a comment) or if it
   //   was an expression that had its perl variables replaced by their
   //   corresponding values
   // 1 if the parsed line was a perl expression that was passed to the perl
   //   interpreter
   // 2 if the parsed line was an expression that had its perl variables
   //   replaced by their corresponding values and the result contains multiple
   //   lines
   // Given this, a line should be output if the OvertureParser returned 0 or 2.
   bool
   outputLine(
      const int& a_result)
   {
      return (a_result == 0 || a_result == 2);
   }

   // The underlying parser.
   OvertureParser* m_ovparser;
};

} // end namespace Loki

#endif
