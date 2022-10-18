###
# example command file for driven vlasovPoisson problem
#   '#' or '#' are comment characters
#   a blank line will stop reading from the command file and move to terminal input
###
#
# we can use PERL to do some math
#
#
$klde = 1/3;
#
$IMass_by_EMass = 100.0;
$ETemp_by_ITemp = 10;
$IAlpha = sqrt($ETemp_by_ITemp)*sqrt($IMass_by_EMass);
$IonZ = 1;
#
$pi = 3.1415926535897932384626;
$xa = -$pi/$klde;
$xb =  $pi/$klde;
$ya = -78*$pi/$klde;
$yb =  78*$pi/$klde;
#
$vEmin = -7;
$vEmax = 7;
#
$vImax =  10/$IAlpha;
$vImin = -10/$IAlpha;
#
#
$xwidth = ($xb-$xa)/2.0;
$ywidth = 300*$pi;
#
# Probe locations are given as fractions of the domain
#
$probe_x1 = 0.5;
$probe_y1 = 0.5;
$probe_x2 = 0.25;
$probe_y2 = 0.5;
$probe_x3 = 0.75;
$probe_y3 = 0.5;
$probe_x4 = 0.95;
$probe_y4 = 0.5;
$probe_x5 = 0.5;
$probe_y5 = 0.25;
$probe_x6 = 0.5;
$probe_y6 = 0.75;
#
#
# now change the default quantities as desired
verbosity = 1
#
# Configuration space Domain
#
domain_limits = $xa $xb $ya $yb 
N = 10 10
periodic_dir = true true
#periodic_dir = false false
#
# Time step controls
#
cfl = 1.0
final_time = 5.
save_times = 1.
sequence_write_times = .1
max_step = 1000000
spatial_solution_order = 6
temporal_solution_order = 6
#
# Plasma definition
#
number_of_species = 2
#
kinetic_species.1.name = "electron"
kinetic_species.1.velocity_limits = $vEmin $vEmax $vEmin $vEmax
kinetic_species.1.Nv = 16 10
kinetic_species.1.mass = 1.0
kinetic_species.1.charge = -1.0
kinetic_species.1.ic.name = "Perturbed Maxwellian"
kinetic_species.1.ic.tx = 1.0 
kinetic_species.1.ic.ty = 1.0
kinetic_species.1.ic.vx0 = 0.0
kinetic_species.1.ic.vy0 = 0.0
kinetic_species.1.ic.A = 0.0
kinetic_species.1.ic.B = 0.0
kinetic_species.1.ic.C = 0.0
kinetic_species.1.ic.kx1 = $klde
kinetic_species.1.ic.ky1 = $klde
#
# external driver for electrons
#
kinetic_species.1.num_external_drivers = 1
kinetic_species.1.external_driver.1.name = "Shaped Ramped Cosine Driver"
kinetic_species.1.external_driver.1.xwidth = $xwidth
kinetic_species.1.external_driver.1.ywidth = $ywidth
kinetic_species.1.external_driver.1.shape = 0.0
kinetic_species.1.external_driver.1.omega = .0381
kinetic_species.1.external_driver.1.E_0 = 0.1
kinetic_species.1.external_driver.1.t_ramp = 1.0
kinetic_species.1.external_driver.1.t_off = 2.0
kinetic_species.1.external_driver.1.x_shape = 0.0
kinetic_species.1.external_driver.1.lwidth = 50
kinetic_species.1.external_driver.1.x0 = 0.0
#
kinetic_species.2.name = "ion"
kinetic_species.2.velocity_limits = $vImin $vImax $vImin $vImax
kinetic_species.2.Nv = 16 10
kinetic_species.2.mass = $IMass_by_EMass
kinetic_species.2.charge = $IonZ
kinetic_species.2.ic.name = "Perturbed Maxwellian"
kinetic_species.2.ic.tx = 0.1
kinetic_species.2.ic.ty = 0.1
kinetic_species.2.ic.vx0 = 0.0
kinetic_species.2.ic.vy0 = 0.0
kinetic_species.2.ic.A = 0.0
kinetic_species.2.ic.B = 0.0
kinetic_species.2.ic.C = 0.0
kinetic_species.2.ic.kx1 = $klde
kinetic_species.2.ic.ky1 = $klde
#
# Diagnostic Controls
#
#
number_of_probes = 6
probe.1.location = $probe_x1 $probe_y1
probe.2.location = $probe_x2 $probe_y2
probe.3.location = $probe_x3 $probe_y3
probe.4.location = $probe_x4 $probe_y4
probe.5.location = $probe_x5 $probe_y5
probe.6.location = $probe_x6 $probe_y6
#
plot_ke_fluxes = true
save_data = true
#
start_from_restart = false
restart.time_interval = 5
restart.write_directory = "planeIAW2D"
