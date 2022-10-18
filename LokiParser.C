#include "LokiParser.H"

extern "C"
{
#include <EXTERN.h>
#include <perl.h>
#include <XSUB.h>
}

namespace Loki {

int LokiParser::s_DEBUG = 1;

//==============================================================================
//  Class for parsing commands using perl
//  Commands containing a semi-colon are processed by perl; otherwise any perl
//  variables (beginning with a '$') are replaced.
//==============================================================================
LokiParser::
LokiParser(
   int argc,
   char** argv)
   : m_parserPointer(0)
{
   PerlInterpreter* const& my_perl_c = (PerlInterpreter*)m_parserPointer;

   PerlInterpreter*& my_perl = (PerlInterpreter*&)my_perl_c;
   my_perl = NULL;

   char** embedding = 0;
   int npa = 3;
   embedding = new char*[npa+argc];
   embedding[0] = (char*)"";
   // execute one line at a time starting with the argument that follows -e
   embedding[1] = (char*)"-e";
   string cmd("use Getopt::Long; use Getopt::Std;");
   embedding[2] = (char*)cmd.c_str();

   // tack on the given argument list, probably should strip out things like
   // noplot/nopause... etc ?
   for (int i = npa; i < (npa+argc); ++i) {
      embedding[i] = argv[i-npa];
   }

   // perldoc.perl.org suggest that PERL_SYS_INIT3() should be invoked before
   // the first interpreter is created and PERL_SYS_TERM() invoked after the
   // last interpreter is freed.  The NULL argument forces PERL_SYS_INIT3() to
   // take the current environment:
   PERL_SYS_INIT3(&argc, &argv, NULL);

   my_perl = perl_alloc();
   perl_construct(my_perl);

   perl_parse((PerlInterpreter*)my_perl, NULL, npa+argc, embedding, NULL);
   perl_run(my_perl);

   m_parserPointer = my_perl;
   delete [] embedding;
}


LokiParser::
~LokiParser()
{
   PerlInterpreter* const& my_perl_c = (PerlInterpreter*)m_parserPointer;

   PerlInterpreter*& my_perl = (PerlInterpreter*&)my_perl_c;
   if(my_perl != NULL) {
      perl_destruct(my_perl);
      perl_free(my_perl);
   }

   PERL_SYS_TERM();
}


// =============================================================================
//   Parse a string with perl.
//   If the string contains a semi-colon, ';', then treat the string as perl
//   commands to process. Otherwise evaluate the string to replace any perl
//   variables.
//
// answer (input/output): On input a string to parse. On output the string with
//                        perl variables replaced IF the input string contained
//                        no semi-colons. If the input string had semi-colons 
//                        then answer remains unchanged on output.
// Return value: 0: answer was not changed or was evaluated replacing any perl
//                  variables with their values.
//               1: answer contained a semi-colon and was processed by perl,
//                  answer remains unchanged. 
//               2: answer contained a newline indicating multiple commands
// =============================================================================
int
LokiParser::parse(
   string & answer)
{
   int returnValue = 0;

   const char* canswer = answer.c_str();
   if (strchr(canswer, ';') != NULL) {
      // answer has a semi-colon -- just evaluate
      eval_pv(canswer, TRUE);

      returnValue = 1;
   }
   else if (strchr(canswer, '$') != NULL) {
      string line = "$overtureParserline = \"" + answer + "\";";
      eval_pv((const char*)line.c_str(), TRUE);
      STRLEN numChars;
      const char* result = SvPV(get_sv("overtureParserline", FALSE), numChars);
      if (s_DEBUG > 0 && strlen(result) < 160) {
         printf("LokiParser::result = [%s]\n", result);
      }

      if (strchr(result, '\n') != NULL) {
         // answer contains a new line
         returnValue=2;
      }
      answer = result;
   }
   return returnValue;
}

} // end namespace Loki
