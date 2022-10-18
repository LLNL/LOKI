/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#ifndef _PERL_PRE_PROCESSOR_H_
#define _PERL_PRE_PROCESSOR_H_

#include "LokiParser.H"

using namespace std;

namespace Loki {

/**
 * Pre-processes Loki input files' perl syntax producing files with all perl
 * variables evaluated and suitable for parsing by Loki.  This scheme is
 * shamelessly stolen from Overture so a Loki input file has the same format
 * as an Overture input file.  The contents are Loki specific but the way that
 * embedded Perl is used is identical to Overture.  Hence, "all" we do here is
 * to use the LokiParser to parse each logical line and, if appropriate, to
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
    *        output file.
    *
    * @param[in] a_input_file_name Name of unprocessed input file.
    *
    * @return Name of processed output file.
    */
   string
   process(
      const string& a_input_file_name);

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

   // This function takes the result of the LokiParser and determines if the
   // output line from the parser should be output.  The LokiParser returns the
   // following values:
   // 0 if the parsed line was unchanged (most likely a comment) or if it
   //   was an expression that had its perl variables replaced by their
   //   corresponding values
   // 1 if the parsed line was a perl expression that was passed to the perl
   //   interpreter
   // 2 if the parsed line was an expression that had its perl variables
   //   replaced by their corresponding values and the result contains multiple
   //   lines
   // Given this, a line should be output if the LokiParser returned 0 or 2.
   bool
   outputLine(
      const int& a_result)
   {
      return (a_result == 0 || a_result == 2);
   }

   // The underlying parser.
   LokiParser* m_parser;
};

} // end namespace Loki

#endif
