#!/usr/local/bin/perl
#########################################################################
##
## Copyright (c) 2018-2022, Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
## Produced at the Lawrence Livermore National Laboratory
##
## SPDX-License-Identifier: MIT
##
## perl script to compare two files but ignore CVS comments.
##
##########################################################################
## Usage: cmp.pl <file1> <file2>
##

$ANAME = shift(@ARGV);
$BNAME = shift(@ARGV);

open(AFILE, "$ANAME") || die "Cannot open input file $ANAME...";
open(BFILE, "$BNAME") || die "Cannot open input file $BNAME...";

while (!eof(AFILE) && !eof(BFILE)) {
   $ALINE = <AFILE>;
   $BLINE = <BFILE>;
   $_ = $ALINE;

   if (!/^(\/\/|c|C|#|##| \*)[ ]*(Release:[\t ]*\$Name|Revision:[\t ]*\$LastChangedRevision|Modified:[\t ]*\$LastChangedDate):[^\$]*\$/o) {
      if ($ALINE ne $BLINE) {
         exit 1;
      }
   }
}

exit 0 if (eof(AFILE) && eof(BFILE));
exit 1;
