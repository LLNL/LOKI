###
# example command file for driven vlasovPoisson problem
#   '#' or '#' are comment characters
#   a blank line will stop reading from the command file and move to terminal input
###
#
# we can use PERL to do some math
#
$pi = 3.1415926535897932384626;
$xa = -2*$pi;
$xb =  2*$pi;
$ya = -2*$pi;
$yb =  2*$pi;
$vmax =  7;
$vmin = -7;
#
#
# now change the default quantities as desired
verbosity = 1
#
# Configuration space Domain
#
domain_limits = $xa $xb $ya $yb 
N = 16 16
#N = 32 32
#N = 64 64
periodic_dir = true true
#periodic_dir = false false
#
# Time step controls
#
cfl = 1.0
final_time = 5.
save_times = .1
sequence_write_times = .1
max_step = 1000000
#
# Plasma definition
#
number_of_species = 1
#
kinetic_species.1.name = "electron"
kinetic_species.1.velocity_limits = $vmin $vmax $vmin $vmax
#kinetic_species.1.Nv = 32 32
kinetic_species.1.Nv = 64 64
kinetic_species.1.mass = 1.0
kinetic_species.1.charge = -1.0
kinetic_species.1.tz.name = "TrigTZSource"
kinetic_species.1.tz.amp = 1
#
# external driver for electrons
#
#
# Diagnostic Controls
#
#
save_data = true
#
start_from_restart = false
restart.time_interval = 1
restart.write_directory = "EPWTZ"
