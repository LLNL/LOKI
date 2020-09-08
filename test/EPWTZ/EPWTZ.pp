$pi = 3.1415926535897932384626;
$xa = -2*$pi;
$xb =  2*$pi;
$ya = -1*$pi;
$yb =  1*$pi;
#
#
$EMass = 1 ; 
$vthE = sqrt(1.0/$EMass) ; 
$EAlpha = 1/$vthE ; 
#
#
$vExmax =  7;
$vExmin = -7;
$vEymax =  9;
$vEymin = -9;
#
#
verbosity = 0
#
#
temporal_solution_order = 6
spatial_solution_order = 6
$Nx = 16;
$Ny = 16;
$Nvx = 16;  
$Nvy = 16;  
#
#
domain_limits = $xa $xb $ya $yb 
N = $Nx $Ny
periodic_dir = true true 
#
#
cfl = 1.0
final_time = 1 
save_times = .2
sequence_write_times = .2
max_step = 1000000 
#
#
number_of_species = 1
#
#
kinetic_species.1.name = "electron" 
kinetic_species.1.velocity_limits = $vExmin $vExmax $vEymin $vEymax 
kinetic_species.1.Nv = $Nvx $Nvy
kinetic_species.1.mass = $EMass 
kinetic_species.1.charge = -1.0
kinetic_species.1.tz.name = "ElectronTrigTZSource"
kinetic_species.1.tz.amp = 0.1 
#
#
kinetic_species.1.ic.name = "Perturbed Maxwellian"
kinetic_species.1.ic.alpha = $EAlpha 
kinetic_species.1.ic.beta = $EAlpha 
kinetic_species.1.ic.vx0 = 0.0
kinetic_species.1.ic.vy0 = 0.0
kinetic_species.1.ic.A = 0.0
kinetic_species.1.ic.B = 0.0
kinetic_species.1.ic.C = 0.0
kinetic_species.1.ic.kx1 = 0.0
kinetic_species.1.ic.ky1 = 0.0
kinetic_species.1.ic.kx2 = 1.0
#
#
save_data = true
#
#
start_from_restart = false
restart.time_interval = .1
restart.write_directory = "TwilightZone"
#
#
