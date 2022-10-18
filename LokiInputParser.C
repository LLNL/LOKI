#include <sstream>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "LokiInputParser.H"
#include "Loki_Defines.H"

namespace Loki {

int LokiInputParser::xargc = -1;
int LokiInputParser::num_obj = 0;
char** LokiInputParser::xargv = 0;
list<LokiInputParser::LIP_entry*> LokiInputParser::table;

static const char* tok_name[] =
{
   "NAME", "OPTION", "INTEGER", "FLOAT", "DOUBLE", "STRING", "EQ_SIGN", "EOF"
};


LokiInputParser::LIP_entry::LIP_entry(
   string& name,
   LokiInputParser::LIPType typ,
   list<string>& vals)
   : defname(name),
     deftype(typ),
     val(vals.size())
{
   list<string>::const_iterator li(vals.begin());
   for (int i = 0; li != vals.end(); ++i, ++li) {
      val[i] = *li;
   }
}


void
LokiInputParser::LIP_entry::dump(
   ostream& os) const
{
   static const char TokenInitial[] =
   {
      'N', 'O', 'I', 'F', 'D', 'S', '=', 'E'
   };

   char tmp[200];
   long nval = val.size();
   sprintf(tmp, "(%c,%1d) %15s :: ",
           TokenInitial[deftype],
           int(nval),
           defname.c_str());
   os << tmp;
   for (int i = 0; i < nval; i++) {
      os << " (" << TokenInitial[LokiInputParser::lipString] << ','
         << val[i] << ')';
   }
   os << endl;

   if (os.fail()) {
      LOKI_ABORT("Failed.");
   }
}


LokiInputParser::LokiInputParser(
   int argc,
   char** argv,
   const char* prefix,
   const char* parfile)
{
   if (!table.empty()) {
      LOKI_ABORT("Table already built");
   }
   num_obj++;
   xargc = argc;
   xargv = argv;
   if (prefix != 0) {
      thePrefix = prefix;
   }
   lipinit(parfile);
}


LokiInputParser::LokiInputParser(
   const char* prefix)
{
   num_obj++;
   if (xargc < 0) {
      LOKI_ABORT("Parser not properly initialized.");
   }
   if (prefix != 0) {
      thePrefix = prefix;
   }
}


void
LokiInputParser::dumpTable(
   ostream& os)
{
   for (list<LIP_entry*>::const_iterator li(table.begin());
        li != table.end(); ++li) {
      (*li)->dump(os);
   }
}


//
// Initialize LokiInputParser.
//
void
LokiInputParser::lipinit(
   const char* parfile)
{
   if (xargc < 0) {
      LOKI_ABORT("Not initiated in main.");
   }

   if (parfile != 0) {
      read_file(parfile, table);
   }

   if (xargc > 0) {
      string argstr;
      const char SPACE = ' ';
      for (int i = 0; i < xargc; i++) {
         argstr += xargv[i];
         argstr += SPACE;
      }
      list<LIP_entry*> arg_table;
      bldTable(argstr.c_str(), argstr.length()+1, arg_table);
      //
      // Append arg_table to end of existing table.
      //
      table.splice(table.end(), arg_table);
   }
}


LokiInputParser::~LokiInputParser()
{
   if (--num_obj == 0) {
      for (list<LIP_entry*>::iterator li(table.begin());
          li != table.end(); ++li) {
         delete *li;
      }
      table.clear();
   }
}

//
// Keyword aware string comparison.
//

//
// Return true if key==keyword || key == prefix.keyword.
//
static
bool
lipfound(
   const char* keyword,
   const string& key,
   const string& prefix)
{
   if (prefix.length() != 0) {
      string tmp(prefix);
      tmp += '.';
      tmp += keyword;
      return (key.compare(tmp) == 0);
   }
   else {
      return (key.compare(keyword) == 0);
   }
}


//
// Return true if name in table.
//
bool
LokiInputParser::contains(
   const char* name)
{
   for (list<LIP_entry*>::const_iterator li(table.begin());
        li != table.end(); ++li) {
      if (lipfound(name, (*li)->defname, thePrefix)) {
         return true;
      }
   }
   return false;
}


//
// Return the index of the n'th occurence of a parameter name,
// except if n==-1, return the index of the last occurence.
// Return 0 if the specified occurence does not exist.
//
const LokiInputParser::LIP_entry*
LokiInputParser::lipindex(
   int n,
   const char* name) const
{
   const LIP_entry* fnd = 0;
   if (n < 0) {
      //
      // Search from back of list.
      //
      for (list<LIP_entry*>::const_reverse_iterator li = table.rbegin();
           li != table.rend(); ++li) {
         if (lipfound(name, (*li)->defname, thePrefix)) {
            fnd = *li;
            break;
         }
      }
   }
   else {
      for (list<LIP_entry*>::const_iterator li(table.begin());
           li != table.end(); ++li) {
         if (lipfound(name, (*li)->defname, thePrefix)) {
            fnd = *li;
            n--;
            if (n < 0) {
               break;
            }
         }
      }
      if (n >= 0) {
         fnd = 0;
      }
   }
   return fnd;
}


void
LokiInputParser::getval(
   const char* name,
   const LIPType type,
   void* ptr,
   int ival,
   int occurence)
{
   if (queryval(name, type, ptr, ival, occurence) == 0) {
      ostringstream msg;
      if (occurence >= 0) {
         msg << "Occurence number " << occurence << " of ";
      }
      if (thePrefix.length() != 0) {
         msg << thePrefix << '.';
      }
      msg << name << " not found in table" << endl;
      dumpTable(cerr);
      LOKI_ABORT(msg.str());
   }
}


int
LokiInputParser::queryval(
   const char* name,
   const LIPType type,
   void* ptr,
   int ival,
   int occurence)
{
   //
   // Get last occurrance of name in table.
   //
   const LIP_entry *def = lipindex(occurence, name);
   if (def == 0) {
      return 0;
   }
   //
   // Does it have ival values?
   //
   if (ival >= def->val.size()) {
      ostringstream msg;
      msg << "No value number" << ival << " for ";
      if (occurence < 0) {
         msg << "last occurence of ";
      }
      else {
         msg << " occurence " << occurence << " of ";
      }
      msg << def->defname << endl;
      def->dump(msg);
      LOKI_ABORT(msg.str());
   }

   const string& valname = def->val[ival];

   int ok;
   double val_dbl;
   //
   // Retrieve value.
   //
   switch (type) {
      case lipInt:
         ok = isInteger(valname, *(int*)ptr);
         break;
      case lipFloat:
         ok = isDouble(valname, val_dbl);
         if (ok) {
            *(float*)ptr = (float) val_dbl;
         }
         break;
      case lipDouble:
         ok = isDouble(valname, val_dbl);
         *(double*)ptr = val_dbl;
         break;
      case lipString:
         ok = true;
         *(string*)ptr = valname;
         break;
      default:
         LOKI_ABORT("Invalid type.");
   }
   if (!ok) {
      ostringstream msg;
      msg << "Type mismatch on value number " << ival << " of " << endl;
      if (occurence < 0) {
         msg << " last occurence of ";
      }
      else {
         msg << " occurence number " << occurence << " of ";
      }
      msg << def->defname << endl;
      msg << " Expected " << tok_name[type]
           << " value = " << valname << endl;
      def->dump(msg);
      LOKI_ABORT(msg.str());
   }
   return 1;
}


void
LokiInputParser::getarr(
   const char* name,
   const LIPType type,
   void* ptr,
   int start_ix,
   int num_val,
   int occurence)
{
   if (queryarr(name, type, ptr, start_ix, num_val, occurence) == 0) {
      ostringstream msg;
      if (occurence >= 0) {
         msg << "Occurence number " << occurence << " of ";
      }
      if (thePrefix.length() != 0) {
         msg << thePrefix << '.';
      }
      msg << name << " not found in table" << endl;
      dumpTable(msg);
      LOKI_ABORT(msg.str());
   }
}


int
LokiInputParser::queryarr(
   const char* name,
   const LIPType type,
   void* ptr,
   int start_ix,
   int num_val,
   int occurence)
{
   //
   // Get last occurrance of name in table.
   //
   const LIP_entry *def = lipindex(occurence, name);
   if (def == 0) {
      return 0;
   }
   //
   // Does it have sufficient number of values and are they all
   // the same type?
   //
   int stop_ix = start_ix + num_val - 1;
   if (stop_ix >= def->val.size()) {
      ostringstream msg;
      msg << "Too many values requested for";
      if (occurence < 0) {
         msg << " last occurence of ";
      }
      else {
         msg << " occurence " << occurence << " of ";
      }
      msg << def->defname << endl;
      def->dump(msg);
      LOKI_ABORT(msg.str());
   }

   for (int n = start_ix; n <= stop_ix; n++) {
      const string& valname = def->val[n];
      //
      // Retrieve value.
      //
      int ok = false;
      double val_dbl;
      switch (type) {
         case lipInt:
            ok = isInteger(valname, *(int*)ptr);
            ptr = (int*)ptr+1;
            break;
         case lipFloat:
            ok = isDouble(valname, val_dbl);
            if (ok) {
               *(float*)ptr = (float) val_dbl;
            }
            ptr = (float*)ptr+1;
            break;
         case lipDouble:
            ok = isDouble(valname, *(double*)ptr);
            ptr = (double*)ptr+1;
            break;
         case lipString:
            ok = true;
            *(string*)ptr = valname;
            ptr = (string*)ptr+1;
            break;
         default:
            LOKI_ABORT("Invalid type");
      }
      if (!ok) {
         ostringstream msg;
         msg << "Type mismatch on value number " <<  n << " of ";
         if (occurence < 0) {
            msg << " last occurence of ";
         }
         else {
            msg << " occurence number " << occurence << " of ";
         }
         msg << def->defname << endl;
         msg << " Expected " << tok_name[type]
             << " value = " << valname << endl;
         def->dump(msg);
         LOKI_ABORT(msg.str());
      }
   }

   return 1;
}


void
LokiInputParser::read_file(
   const char* fname,
   list<LIP_entry*>& tab)
{
   //
   // Space for input file if it exists.
   //
   if (fname != 0 && fname[0] != 0) {
      FILE* pffd = fopen(fname, "r");
      if (pffd == 0) {
         ostringstream msg;
         msg << "Couldn't open \"" << fname  << "\"" << endl;
         LOKI_ABORT(msg.str());
      }
      //
      // Get the length.
      //
      fseek(pffd, 0, 2);
      int pflen = (int)ftell(pffd);
      rewind(pffd);
      char* str = new char[pflen+1];
      if (str == 0) {
         LOKI_ABORT("Out of memory, input file too long.");
      }
      memset(str, 0, pflen+1);
      int nread = fread(str, 1, pflen, pffd);
      if (nread != pflen) {
         ostringstream msg;
         msg << "Only read "  << nread
             << " bytes out of "
             << pflen
             << " from "
             << fname << endl;
         LOKI_ABORT(msg.str());
      }
      fclose(pffd);
      bldTable(str, pflen+1, tab);
      delete [] str;
   }
}


static
void
eat_garbage(
   const char* str,
   int& i,
   int len)
{
   for (;;) {
      if (i < len && str[i] == '#') {
         while (i < len && str[i] != '\n') {
            i++;
         }
      }
      else if (i < len && isspace(str[i])) {
         i++;
      }
      else {
         break;
      }
   }
}


//
// Simple lexical analyser.
//

enum lexState
{
    START,
    MINUS,
    SIGN,
    OPTION,
    STRING,
    QUOTED_STRING,
    INTEGER,
    START_FRACTION,
    FRACTION,
    START_EXP,
    SIGNED_EXP,
    EXP,
    PREFIX,
    SUFFIX,
    STAR
};

static const char* state_name[] =
{
   "START",
   "MINUS",
   "SIGN",
   "OPTION",
   "STRING",
   "QUOTED_STRING",
   "INTEGER",
   "START_FRACTION",
   "FRACTION",
   "START_EXP",
   "SIGNED_EXP",
   "EXP",
   "PREFIX",
   "SUFFIX",
   "STAR"
};

LokiInputParser::LIPType
LokiInputParser::getToken(
   const char* str,
   int& i,
   int slen,
   char* ostr)
{
#define ERROR_MESS \
   ostr[k++] = '\0'; \
   ostringstream msg; \
   msg << "Invalid string = " << ostr << endl; \
   msg << "STATE = " << state_name[state] \
       << ", next char = " << ch << endl; \
   msg << ", rest of input = \n" << (str+i) << endl; \
   LOKI_ABORT(msg.str())

   //
   // Eat white space and comments.
   //
   eat_garbage(str, i, slen);

   //
   // Check for end of file.
   //
   if (i >= slen || str[i] == '\0') {
      return lipEOF;
   }

   //
   // Start token scan.
   //
   lexState state = START;
   int k = 0;
   while (1) {
      char ch = str[i];
      switch (state)
      {
         case START:
            if (ch == '=') {
               ostr[k++] = ch;
               i++;
               ostr[k++] = 0;
               return lipEQ_sign;
            }
            else if (ch == '"') {
               i++;
               state = QUOTED_STRING;
            }
            else if (ch == '*') {
               ostr[k++] = ch;
               i++;
               state = STAR;
            }
            else if (isalpha(ch) || ch == '_') {
               ostr[k++] = ch;
               i++;
               state = PREFIX;
            }
            else if (ch == '-') {
               ostr[k++] = ch;
               i++;
               state = MINUS;
            }
            else {
               ostr[k++] = ch;
               i++;
               state = STRING;
            }
            break;

         case MINUS:
            if (isalpha(ch) || ch == '_') {
               k--;
               ostr[k++] = ch;
               i++;
               state = OPTION;
            }
            else {
               ostr[k++] = ch;
               i++;
               state = STRING;
            }
            break;

         case OPTION:
            if (isalnum(ch) || ch == '_' || ch == '.') {
               ostr[k++] = ch;
               i++;
            }
            else if (isspace(ch)) {
               ostr[k++] = 0;
               return lipOption;
            }
            else {
               ERROR_MESS;
            }
            break;

         case STAR:
            if (ch == '.') {
               ostr[k++] = ch;
               i++;
               state = SUFFIX;
            }
            else {
               ostr[k++] = ch;
               i++;
               state = STRING;
            }
            break;

         case PREFIX:
            if (isalnum(ch) || ch == '_') {
               ostr[k++] = ch;
               i++;
            }
            else if (ch == '.') {
               ostr[k++] = ch;
               i++;
               state = SUFFIX;
            }
            else if (isspace(ch) || ch == '=') {
               ostr[k++] = 0;
               return lipDefn;
            }
            else {
               ostr[k++] = ch;
               i++;
               state = STRING;
            }
            break;

         case SUFFIX:
            if (isalnum(ch) || ch == '_') {
               ostr[k++] = ch;
               i++;
            }
            else if (isspace(ch) || ch == '=') {
               ostr[k++] = 0;
               return lipDefn;
            }
            else {
               ostr[k++] = ch;
               i++;
               state = STRING;
            }
            break;

         case QUOTED_STRING:
            if (ch == '"') {
               i++;
               ostr[k++] = 0;
               return lipString;
            }
            else {
               ostr[k++] = ch;
               i++;
            }
            break;

         case STRING:
            if (isspace(ch) || ch == '=') {
               ostr[k++] = 0;
               return lipString;
            }
            else {
               ostr[k++] = ch;
               i++;
            }
            break;

         default:
            ERROR_MESS;
      }
   }
#undef ERROR_MESS
}


void
LokiInputParser::bldTable(
   const char* str,
   int lenstr,
   list<LIP_entry*>& tab)
{
   string cur_name;
   list<string> cur_list;
   string tmp_str;
   LIP_entry *lip;

   int i = 0;
   LIPType token;
   const int SCRATCH_STR_LEN  = 200;
   char tokname[SCRATCH_STR_LEN];

   while (true) {
      token = getToken(str, i, lenstr, tokname);

      switch (token) {
         case lipEOF:
            addDefn(cur_name, cur_list, tab);
            return;
         case lipOption:
            addDefn(cur_name, cur_list, tab);
            tmp_str = tokname;
            lip = new LIP_entry(tmp_str, lipOption, cur_list);
            if (lip == 0) {
               LOKI_ABORT("Out of memory, too many table entries.");
            }
            tab.push_back(lip);
            break;
         case lipEQ_sign:
            if (cur_name.length() == 0) {
               LOKI_ABORT("EQ with no current defn.");
            }
            if (cur_list.empty()) {
               //
               // First time we see equal sign, just ignore.
               //
               break;
            }
            //
            // Read one too far, remove last name on list.
            //
            tmp_str = cur_list.back();
            cur_list.pop_back();
            addDefn(cur_name, cur_list, tab);
            cur_name = tmp_str;
            break;
         case lipDefn:
            if (cur_name.length() == 0) {
               cur_name = tokname;
               break;
            }
         //
         // Otherwise, fall through.  This may be a string.
         //
         case lipInt:
         case lipFloat:
         case lipDouble:
         case lipString:
            if (cur_name.length() == 0) {
               LOKI_ABORT("Value with no defn.");
            }
            cur_list.push_back(tokname);
            break;
      }
   }
}


void
LokiInputParser::addDefn(
   string& def,
   list<string>& val,
   list<LIP_entry*>& tab)
{
   static const string FileKeyword("FILE");
   //
   // Check that defn exists.
   //
   if (def.length() == 0) {
      val.clear();
      return;
   }

   //
   // Check that it has values.
   //
   if (val.empty()) {
      ostringstream msg;
      msg << "No values for definition " << def << endl;
      LOKI_ABORT(msg.str());
   }

   //
   // Check if this defn is a file include directive.
   //
   if (def == FileKeyword && val.size() == 1) {
      //
      // Read file and add to this table.
      //
      const char* fname = val.front().c_str();
      read_file(fname, tab);
   }
   else {
      LIP_entry* lip = new LIP_entry(def, lipDefn, val);
      if (lip == 0) {
         LOKI_ABORT("Out of memory, too many table entries.");
      }
      tab.push_back(lip);
   }
   val.clear();
   def = string();
}

} // end namespace Loki
