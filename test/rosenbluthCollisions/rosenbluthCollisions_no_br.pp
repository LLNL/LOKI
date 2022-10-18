$klde = 1.0/3.0;
#
$pi = 3.1415926535897932384626;
$xa = -.8;
$xb =  .8;
$ya = -.8;
$yb =  .8;
#
$vHEmin = -4;
$vHEmax = 4;
$vEmin = -8;
$vEmax = 8;
$m1 = 4.0;
#
# Probe locations are given as fractions of the domain
#
$probe_x1 = 0.5;
$probe_y1 = 0.5;
$probe_x2 = 0.25;
$probe_y2 = 0.5;
$probe_x3 = 0.75;
$probe_y3 = 0.5;
$probe_x4 = 0.0;
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
N = 16 16
periodic_dir = true true
#periodic_dir = false false
#
# Time step controls
#
cfl = 0.9
# final_time = 1.0
max_step=30
save_times = 0.01
coll_save_steps = 1
# coll_save_times = 0.01
coll_seq_write_steps = 1
# coll_seq_write_times = 0.01
sequence_write_times = 0.01
#
# Plasma definition
#
number_of_species = 2
kinetic_species.1.name = "heavyelectron"
kinetic_species.1.velocity_limits = $vHEmin $vHEmax $vHEmin $vHEmax
kinetic_species.1.Nv = 64 64
kinetic_species.1.mass = $m1
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
kinetic_species.1.ic.ky1 = 0.0
kinetic_species.1.ic.vflowinitx = 0.5
kinetic_species.1.ic.vflowinity = 0.5
kinetic_species.1.num_collision_operators = 1
kinetic_species.1.collision_operator.1.name = "Rosenbluth Collision Operator"
kinetic_species.1.collision_operator.1.collision_interspecies = true
kinetic_species.1.collision_operator.1.collision_nuCoeffInterspecies = 0.1
kinetic_species.1.collision_operator.1.collision_back_reaction = false
kinetic_species.1.collision_operator.1.collision_self = false
kinetic_species.1.collision_operator.1.collision_diagnostic = true
kinetic_species.1.collision_operator.1.kernel_alg = "primitive"
kinetic_species.1.collision_operator.1.collision_vel_range_lo = -3.25 -3.25
kinetic_species.1.collision_operator.1.collision_vel_range_hi = 3.25 3.25

#
kinetic_species.2.name = "electron"
kinetic_species.2.velocity_limits = $vEmin $vEmax $vEmin $vEmax
kinetic_species.2.Nv = 64 64
kinetic_species.2.mass = 1.0
kinetic_species.2.charge = -1.0
kinetic_species.2.ic.name = "Perturbed Maxwellian"
kinetic_species.2.ic.tx = 1.0
kinetic_species.2.ic.ty = 1.0
kinetic_species.2.ic.vx0 = 0.0
kinetic_species.2.ic.vy0 = 0.0
kinetic_species.2.ic.A = 0.0
kinetic_species.2.ic.B = 0.0
kinetic_species.2.ic.C = 0.0
kinetic_species.2.ic.kx1 = $klde
kinetic_species.2.ic.ky1 = 0.0
kinetic_species.2.ic.vflowinitx = -1.0
kinetic_species.2.ic.vflowinity = -1.0
kinetic_species.2.num_collision_operators = 1
kinetic_species.2.collision_operator.1.name = "Rosenbluth Collision Operator"
kinetic_species.2.collision_operator.1.collision_interspecies = true
kinetic_species.2.collision_operator.1.collision_nuCoeffInterspecies = 0.1
kinetic_species.2.collision_operator.1.collision_back_reaction = false
kinetic_species.2.collision_operator.1.collision_self = false
kinetic_species.2.collision_operator.1.collision_diagnostic = true
kinetic_species.2.collision_operator.1.kernel_alg = "primitive"
kinetic_species.2.collision_operator.1.collision_vel_range_lo = -6.5 -6.5
kinetic_species.2.collision_operator.1.collision_vel_range_hi = 6.5 6.5
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
save_coll_data = true
#
start_from_restart = false
restart.time_interval = 0.05
restart.write_directory = "rosenbluthCollisions_no_br"
