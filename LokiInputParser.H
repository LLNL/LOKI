#ifndef _LOKI_INPUT_PARSER_H_
#define _LOKI_INPUT_PARSER_H_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>

using namespace std;

namespace Loki {

class LokiInputParser
{
   enum LIPType
   {
      lipDefn,
      lipOption,
      lipInt,
      lipFloat,
      lipDouble,
      lipString,
      lipEQ_sign,
      lipEOF
   };

private:
   class LIP_entry
   {
   private:
      friend class LokiInputParser;
      LIP_entry()
      {
      }

      LIP_entry(
         string& name,
         LokiInputParser::LIPType typ,
         list<string>& vals);

      ~LIP_entry()
      {
      }

      string defname;

      LokiInputParser::LIPType deftype;

      vector<string> val;

      void dump(
         ostream& os) const;
   };

public:
   LokiInputParser(
      int argc,
      char** argv,
      const char* prefix = 0,
      const char* parfile = 0);

   LokiInputParser(
      const char* prefix = 0);

   ~LokiInputParser();

   bool
   contains(
      const char* name);

   int
   countval(
      const char* name,
      int n = -1)
   {
      //
      // First find n'th occurance of name in table.
      //
      const LIP_entry* def = lipindex(n, name);
      return def == 0 ? 0 : def->val.size();
   }

   void
   dumpTable(
      ostream& os);

   void
   get(
      const char* name,
      int& ref,
      int ival = 0)
   {
      getval(name, lipInt, &ref, ival, -1);
   }

   int
   query(
      const char* name,
      int& ref,
      int ival = 0)
   {
      return queryval(name, lipInt, &ref, ival, -1);
   }

   void
   get(
      const char* name,
      float& ref,
      int ival = 0)
   {
      getval(name, lipFloat, &ref, ival, -1);
   }

   int
   query(
      const char* name,
      float& ref,
      int ival = 0)
   {
      return queryval(name, lipFloat, &ref, ival, -1);
   }

   void
   get(const char* name,
      double& ref,
      int ival = 0)
   {
      getval(name, lipDouble, &ref, ival, -1);
   }

   int
   query(
      const char* name,
      double& ref,
      int ival = 0)
   {
      return queryval(name, lipDouble, &ref, ival, -1);
   }

   void
   get(
      const char* name,
      string& ref,
      int ival = 0)
   {
      getval(name, lipString, &ref, ival, -1);
   }

   int
   query(
      const char* name,
      string& ref,
      int ival = 0)
   {
      return queryval(name, lipString, &ref, ival, -1);
   }

   void
   getarr(
      const char* name,
      vector<int>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      getarr(name, lipInt, &(ref[0]), start_ix, num_val, -1);
   }

   int
   queryarr(
      const char* name,
      vector<int>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      return queryarr(name, lipInt, &(ref[0]), start_ix, num_val, -1);
   }

   void
   getarr(
      const char* name,
      vector<float>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      getarr(name, lipFloat, &(ref[0]), start_ix, num_val, -1);
   }

   int
   queryarr(
      const char* name,
      vector<float>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      return queryarr(name, lipFloat, &(ref[0]), start_ix, num_val, -1);
   }

   void
   getarr(
      const char* name,
      vector<double>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      getarr(name, lipDouble, &(ref[0]), start_ix, num_val, -1);
   }

   int
   queryarr(
      const char* name,
      vector<double>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      return queryarr(name, lipDouble, &(ref[0]), start_ix, num_val, -1);
   }

   void
   getarr(
      const char* name,
      vector<string>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      getarr(name, lipString, &(ref[0]), start_ix, num_val, -1);
   }

   int
   queryarr(
      const char* name,
      vector<string>& ref,
      int start_ix,
      int num_val)
   {
      if (ref.size() < num_val) {
         ref.resize(num_val);
      }
      return queryarr(name, lipString, &(ref[0]), start_ix, num_val, -1);
   }

private:
   //
   // Table of entries common to all objects.
   //
   static list<LIP_entry*> table;

   //
   // Command line arguments.
   //
   static int xargc;
   static char** xargv;

   //
   // Keep track of number of LokiInputParser objects out there.
   //
   static int num_obj;

   //
   // Parses string and builds table.
   //
   void
   bldTable(
      const char* str,
      int lenstr,
      list<LIP_entry*>& tab);

   //
   // Add defn to table, check for file inclusion.
   //
   void
   addDefn(
      string& def,
      list<string>&  val,
      list<LIP_entry*>& tab);

   //
   // Reads file into string then parses it with call to bldTable.
   //
   void
   read_file(
      const char* fname,
      list<LIP_entry*>& tab);

   //
   // Lexical analyser called by bldTable.
   //
   LIPType
   getToken(
      const char* str,
      int& i,
      int slen,
      char* ostr);

   //
   // Prefix used in keyword search.
   //
   string thePrefix;

   //
   // Used by constructor to build table.
   //
   void
   lipinit(
      const char* parfile);

   //
   // Find n'th occurence of name in table.
   //
   const LIP_entry*
   lipindex(
      int n,
      const char* name) const;

   //
   // Get ival_th value of k_th occurence of given name.
   // If k_th occurence does not exist or ival_th value does
   // not exist or if ival_type does not match with type, an
   // error message is reported and the program halts.
   // If successful, value is stored in ptr.
   // same as above but searches for last occurence of name.
   //
   void
   getval(
      const char* name,
      const LIPType type,
      void* ptr,
      int ival,
      int k = -1);

   //
   // Get an array of values.
   //
   void
   getarr(
      const char* name,
      const LIPType type,
      void* ptr,
      int start_ix,
      int num_val,
      int k = -1);

   int
   queryval(
      const char* name,
      const LIPType type,
      void* ptr,
      int ival,
      int k = -1);

   int
   queryarr(
      const char* name,
      const LIPType type,
      void* ptr,
      int start_ix,
      int num_val,
      int k = -1);

   bool
   isInteger(
      const string& str,
      int& val)
   {
      char* endp = 0;
      val = (int) strtol(str.c_str(), &endp, 10);
      return *endp == 0;
   }

   int
   isDouble(
      const string& str,
      double& val)
   {
      char* endp = 0;
      val = strtod(str.c_str(), &endp);
      return *endp == 0;
   }
};

} // end namespace Loki

#endif
