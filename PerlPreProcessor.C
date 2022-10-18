/*************************************************************************
 *
 * Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
 * See the top-level LICENSE file for details.
 * Produced at the Lawrence Livermore National Laboratory
 *
 * SPDX-License-Identifier: MIT
 *
 ************************************************************************/
#include "PerlPreProcessor.H"
#include "Loki_Defines.H"
#include "Loki_Utilities.H"

#include <iostream>
#include <fstream>
#include <sstream>

namespace Loki {

PerlPreProcessor::PerlPreProcessor(
   int a_argc,
   char** a_argv)
{
   // Construct the LokiParser which is what does the actual parsing.
   m_parser = new LokiParser(a_argc, a_argv);
}


PerlPreProcessor::~PerlPreProcessor()
{
   delete m_parser;
}


string
PerlPreProcessor::process(
   const string& a_input_file_name)
{
   // Make an input file stream for the input file.
   ifstream ifile(a_input_file_name.c_str());
   if (!ifile) {
      ostringstream msg;
      msg << "Input file " << a_input_file_name << " could not be opened!";
      LOKI_ABORT(msg);
   }

   // The processed output file is named the same as the input but with
   // .processed tacked on.
   string output_file_name(a_input_file_name + ".processed");

   // Only process 0 preprocesses the input file.
   if (Loki_Utilities::s_my_id == 0) {
      // Make an output file stream for the output file.
      ofstream ofile(output_file_name.c_str());
      if (!ofile) {
         ostringstream msg;
         msg << "Output file " << output_file_name << " could not be opened!";
         LOKI_ABORT(msg);
      }

      // Read and process each line of the input file until the end of file is
      // reached.
      while (!ifile.eof()) {
         // Get the physical line from the file.
         string physical_line;
         getline(ifile, physical_line);

         // Set the logical line to the physical line.  If the phyical line is
         // not empty and it contained a line continuation character then the
         // logical line consists of multiple physical lines.  Read each of the
         // logical line's remaining physical lines.
         string logical_line = physical_line;
         while (physical_line.size() > 0 &&
                physical_line.rfind('\\') == physical_line.size()-1) {
            // Get rid of the line continuation character in the logical line,
            // get the next physical line, and append it to the logical line.
            logical_line.resize(logical_line.size()-1);
            getline(ifile, physical_line);
            logical_line += physical_line;
         }

         // Parse the logical line with the LokiParser and output it if
         // appropriate.
         int result(m_parser->parse(logical_line));
         if (outputLine(result)) {
            ofile << logical_line << endl;
         }
      }
      ofile.close();
      fflush(0);
   }

   ifile.close();

   return output_file_name;
}

} // end namespace Loki
