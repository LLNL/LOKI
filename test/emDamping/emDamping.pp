$omega = 3.16;
$clight = 22.36;
$klde = sqrt($omega**2-1)/$clight;
#
$Ey = 1.e-4;
$Bz = $klde*$Ey/$omega;
$uy = -$Ey/$omega;
#
$pi = 3.1415926535897932384626;
$pihalf = $pi/2.0;
#
$time = 10.0;
#
$xa = -$pi/$klde;
$xb =  $pi/$klde;
$ya = -10.;
$yb =  10.;
#
$vEmin = -7;
$vEmax = 7;
#
# Probe locations are given as fractions of the domain
#
$probe_x1 = 0.0;
$probe_x2 = 0.25;
$probe_x3 = 0.5;
$probe_x4 = 0.75;
$probe_x5 = 1.0;
$probe_y1 = 0.5;
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
final_time = $time
save_times = 0.2
sequence_write_times = .1
max_step = 1000000
light_speed = $clight
#
# System type
#
sys_type = "maxwell"
#
# Plasma definition
#
number_of_species = 1
#
kinetic_species.1.name = "electron"
kinetic_species.1.velocity_limits = $vEmin $vEmax $vEmin $vEmax
kinetic_species.1.Nv = 64 64
kinetic_species.1.mass = 1.0
kinetic_species.1.charge = -1.0
kinetic_species.1.ic.name = "Perturbed Maxwellian"
kinetic_species.1.ic.alpha = 1.0 
kinetic_species.1.ic.beta = 1.0
kinetic_species.1.ic.vx0 = 0.0
kinetic_species.1.ic.vy0 = $uy
kinetic_species.1.ic.A = 0.0
kinetic_species.1.ic.B = 0.0
kinetic_species.1.ic.C = 0.0
kinetic_species.1.ic.kx1 = 0.0
kinetic_species.1.ic.ky1 = 0.0
kinetic_species.1.ic.x_wave_number = $klde
kinetic_species.1.ic.y_wave_number = 0.0
kinetic_species.1.ic.phase = $pihalf
#
# Field solve parameters
#
maxwell.avWeak = 0.0
maxwell.avStrong = 0.1
maxwell.em_ic.1.name = "SimpleEMIC"
maxwell.em_ic.1.field = "E"
maxwell.em_ic.1.xamp = 0.0
maxwell.em_ic.1.yamp = $Ey
maxwell.em_ic.1.zamp = 0.0
maxwell.em_ic.1.x_wave_number = $klde
maxwell.em_ic.1.y_wave_number = 0.0
maxwell.em_ic.1.phase = 0.0
maxwell.em_ic.2.name = "SimpleEMIC"
maxwell.em_ic.2.field = "B"
maxwell.em_ic.2.xamp = 0.0
maxwell.em_ic.2.yamp = 0.0
maxwell.em_ic.2.zamp = $Bz
maxwell.em_ic.2.x_wave_number = $klde
maxwell.em_ic.2.y_wave_number = 0.0
maxwell.em_ic.2.phase = 0.0
maxwell.vel_ic.1.name = "SimpleVELIC"
maxwell.vel_ic.1.amp = 0.0
maxwell.vel_ic.1.x_wave_number = 0.0
maxwell.vel_ic.1.y_wave_number = 0.0
maxwell.vel_ic.1.phase = 0.0
#
# Diagnostic Controls
#
#
number_of_probes = 5
probe.1.location = $probe_x1 $probe_y1
probe.2.location = $probe_x2 $probe_y1
probe.3.location = $probe_x3 $probe_y1
probe.4.location = $probe_x4 $probe_y1
probe.5.location = $probe_x5 $probe_y1
#
plot_ke_fluxes = true
#
start_from_restart = false
restart.time_interval = 5
restart.write_directory = "emDamping"
