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
#include "PerlPreProcessor.H"

#include "Overture.h"

#include <iostream>
#include <fstream>

namespace Loki {

PerlPreProcessor::PerlPreProcessor(
   int a_argc,
   char** a_argv)
{
   // Construct the OvertureParser which is what does the actual parsing.
   m_ovparser = new OvertureParser(a_argc, a_argv);
}


PerlPreProcessor::~PerlPreProcessor()
{
   delete m_ovparser;
}


std::string
PerlPreProcessor::process(
   const std::string& a_input_file_name)
{
   // Make an input file stream for the input file.
   std::ifstream ifile(a_input_file_name.c_str());
   if (!ifile) {
      aString msg("Input file " + a_input_file_name + " could not be opened!");
      Overture::abort(msg);
   }

   // The processed output file is named the same as the input but with
   // .processed tacked on.
   std::string output_file_name(a_input_file_name + ".processed");

   int myid = max(0, Communication_Manager::My_Process_Number);

   // Only process 0 preprocesses the input file.
   if (myid == 0) {
      // Make an output file stream for the output file.
      std::ofstream ofile(output_file_name.c_str());
      if (!ofile) {
         aString msg("Output file " + output_file_name + " could not be opened!");
         Overture::abort(msg);
      }

      // Read and process each line of the input file until the end of file is
      // reached.
      while (!ifile.eof()) {
         // Get the physical line from the file.
         aString physical_line;
         getline(ifile, physical_line);

         // Set the logical line to the physical line.  If the phyical line is
         // not empty and it contained a line continuation character then the
         // logical line consists of multiple physical lines.  Read each of the
         // logical line's remaining physical lines.
         aString logical_line = physical_line;
         while (physical_line.size() > 0 &&
                physical_line.rfind('\\') == physical_line.size()-1) {
            // Get rid of the line continuation character in the logical line,
            // get the next physical line, and append it to the logical line.
            logical_line.resize(logical_line.size()-1);
            getline(ifile, physical_line);
            logical_line += physical_line;
         }

         // Parse the logical line with the OvertureParser and output it if
         // appropriate.
         int result(m_ovparser->parse(logical_line));
         if (outputLine(result)) {
            ofile << logical_line << std::endl;
         }
      }
      ofile.close();
      fflush(0);
   }

   ifile.close();

   return output_file_name;
}

} // end namespace Loki
