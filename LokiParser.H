#ifndef _LOKI_PARSER_H_
#define _LOKI_PARSER_H_

#include <string>

using namespace std;

namespace Loki {

class LokiParser
{
public:
   LokiParser(
      int argc = 0,
      char** argv = 0);

   ~LokiParser();

   int
   parse(
      string& answer);

   // set to a positive value to get debug info
   static int s_DEBUG;

private:
   void* m_parserPointer;
};

} // end namespace Loki

#endif
