# Project Title
Loki is a parallel code for the solution of the coupled 2D-2V Vlasov-Poisson or
Vlasov-Maxwell equations.  Inclusion of the Fokker-Planck term for particle
collisions is a user option.  Loki consists of 2 separate programs.  The first
is the actual parallel simulation code and is named vlasovPoisson4D and the
second is a serial post-processing tool named vp4DPostProcess.  Despite the
name, the simulation code is able to solve either the Vlasov-Poisson or the
Vlasov-Maxwell system depending on the user input deck specification.

## Getting Started
Loki depends on several external packages.  Complete instructions for how to
build these packages, how to build Loki itself, and how to run Loki may be
found in docs/Loki_Tutuorial.pdf.  This file also includes detailed information
about the input deck options and syntax.  There are several regression tests in
the test directory which cover several different basic problem types and
configurations.  These are a good starting point to understand how to set up a
Loki input deck.

## Authors
Jeffrey Banks banksj3@rpi.edu (Rensselaer Polytechnic Institute, Amos Eaton 301, 110 8th St., Troy, NY 12180)  
Jeffrey Hittinger hittinger1@llnl.gov (LLNL, P.O Box 808, Livermore, CA 94551)  
William Arrighi arrighi2@llnl.gov (LLNL, P.O Box 808, Livermore, CA 94551)  
Richard Berger berger5@llnl.gov (LLNL, P.O Box 808, Livermore, CA 94551)  
Thomas Chapman chapman29@llnl.gov (LLNL, P.O Box 808, Livermore, CA 94551)  
Stephan Brunner stephan.brunner@epfl.ch (Ecole Polytechnique Federale de Lausanne, EPFL SB SPC-TH, PPB 312, Station 13, CH-1015 Lausanne, Switzerland)  

## License
Loki is distributed under the terms of the MIT license.  See LICENSE for
details.

SPDX-License-Identifier: MIT

LLNL-CODE-840411
