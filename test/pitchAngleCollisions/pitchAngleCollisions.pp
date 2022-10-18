$klde = .3;
#
$pi = 3.1415926535897932384626;
$xa = -$pi/$klde;
$xb =  $pi/$klde;
$ya = -.1;
$yb =  .1;
#
$vEmin = -7;
$vEmax = 7;
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
verbosity = 0
#
# Configuration space Domain
#
domain_limits = $xa $xb $ya $yb 
N = 32 5
periodic_dir = true true
#periodic_dir = false false
#
# Time step controls
#
cfl = 0.8
final_time = 5.0
save_times = 1.0
sequence_write_times = .1
max_step = 1000000
#
# Plasma definition
#
number_of_species = 1
kinetic_species.1.name = "electron"
kinetic_species.1.velocity_limits = $vEmin $vEmax $vEmin $vEmax
kinetic_species.1.Nv = 64 64
kinetic_species.1.mass = 1.0
kinetic_species.1.charge = -1.0
kinetic_species.1.ic.name = "Perturbed Maxwellian"
kinetic_species.1.ic.tx = 1.0 
kinetic_species.1.ic.ty = 1.0
kinetic_species.1.ic.vx0 = 0.0
kinetic_species.1.ic.vy0 = 0.0
kinetic_species.1.ic.A = 0.0001
kinetic_species.1.ic.B = 0.0
kinetic_species.1.ic.C = 0.0
kinetic_species.1.ic.kx1 = $klde
kinetic_species.1.ic.ky1 = $klde
kinetic_species.1.num_collision_operators = 1
kinetic_species.1.collision_operator.1.name = "Pitch Angle Collision Operator"
kinetic_species.1.collision_operator.1.collision_vel_range_lo = -5.25 -5.25
kinetic_species.1.collision_operator.1.collision_vel_range_hi = 5.25 5.25
kinetic_species.1.collision_operator.1.collision_vfloor = 0.01
kinetic_species.1.collision_operator.1.collision_vthermal_dt = 1.0
kinetic_species.1.collision_operator.1.collision_nuCoeff = 0.1
kinetic_species.1.collision_operator.1.collision_conservative = 1
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
restart.write_directory = "pitchAngleCollisions"
