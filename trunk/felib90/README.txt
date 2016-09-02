=======================================================================

@@@@@@@ @@@@@@@ @         @@@   @@@@@@          @@@@      @@@@
@       @       @          @    @     @        @    @    @    @
@       @       @          @    @     @        @    @   @      @
@@@@@   @@@@@   @          @    @@@@@@          @@@@@   @      @
@       @       @          @    @     @             @   @      @
@       @       @          @    @     @             @    @    @
@       @@@@@@@ @@@@@@@   @@@   @@@@@@          @@@@      @@@@

The Finite Element Library 90


(C) 2005 Council for the Central Laboratory of the Research Councils 
       
         Rutherford Appleton Laboratory  (CCLRC)

=======================================================================

The Finite Element Library (FELIB) is a program and subroutine 
library for the solution of partial differential equations using 
the finite element method.

This software is for the use of the academic community in their
personal research and not for commercial exploitation. It comes
with no guarantees in any form concerning it correctness or fitness
for purpose.

Anyone wishing to make commercial use the software should contact
the Head of the Software Engineering Group.

=======================================================================

Please address all bugs or problems to author : 
 
   Prof Chris Greenough
   Head, Software Engineering Group
   Rutherford Appleton Laboratory
   Chilton
   DIDCOT OX11 0QX

   Tel: +44 (0) 1235 445307
   Fax: +44 (0) 1235 446626
   Email: c.greenough@rutherford.ac.uk

=======================================================================

1.Introduction

FELIB is very simple to install and only requires a Fortran 90/95 compiler.

The library has the following directories:

* globals90 - global parameters
* utils90 - utilities routines 
* machine90 - machine dependent routines (defaults should work ok)
* lib90 - basic library routines
* prog90 - example programs
* data - example data
* results - results generated on Redhat 7.2 using g77
* README.txt - this file

2.Installation

The version on the web has been compiled and tested under Linux 
(Redhat 7.2/Redhat WS 3/Fedora Core 2) using Lahey Fortran 95 compiler.

* Copy the tar file to a suitable area
* Extract the contents - this will produce a directory felib90
* Change directory to felib90
* You may need to modify the machine dependent routines in directory
  machine with the necessary values for your system.
* Set up the Makefile.compilers file and the Makefile.<compiler> if required 
* Type make clean, then make FC=<compiler> - the default is g95
* Change directory to programs 
* Type ./test-run.sh - this will run all the programs and compare the
  results with those provide with 'diff'. If the files are not empty
  have a look at the differences.

3.Help

You can always drop me an Email and I will do what I can.

Best of luck

Chris Greenough

