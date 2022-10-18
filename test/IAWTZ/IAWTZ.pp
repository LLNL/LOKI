###
#
#
#
$pi = 3.1415926535897932384626;
$xa = -2*$pi;
$xb =  2*$pi;
$ya = -1*$pi;
$yb =  1*$pi;
#
#
$IMass = 10;
$IT = 1.0; 
$EMass = 1 ;
$ET = 1.0; 
$IAlpha = sqrt($IMass) ;  
$EAlpha = sqrt($EMass) ; 
$vthI = 1.0/$IAlpha ; 
$vthE = 1.0/$EAlpha ; 
#
#
$vExmax =  7*$vthE;
$vExmin = -7*$vthE;
$vEymax =  9*$vthE;
$vEymin = -9*$vthE;
#
$vIxmax =  7*$vthI;
$vIxmin = -7*$vthI;
$vIymax =  9*$vthI;
$vIymin = -9*$vthI;
#
verbosity = 0
#
#
temporal_solution_order = 6
spatial_solution_order = 6
$Nx = 16 ;
$Ny = 16 ; 
$Nvx = 16 ; 
$Nvy = 16 ; 
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
number_of_species = 2
#
#
kinetic_species.1.name = "electron" 
kinetic_species.1.velocity_limits = $vExmin $vExmax $vEymin $vEymax 
kinetic_species.1.Nv = $Nvx $Nvy
kinetic_species.1.mass = $EMass 
kinetic_species.1.charge = -1.0
kinetic_species.1.tz.name = "TwoSpecies_ElectronTrigTZSource"
kinetic_species.1.tz.amp = 0.1 
kinetic_species.1.tz.electron_mass = $EMass
kinetic_species.1.tz.ion_mass = $IMass
#
#
kinetic_species.1.ic.name = "Perturbed Maxwellian"
kinetic_species.1.ic.tx = $ET 
kinetic_species.1.ic.ty = $ET 
kinetic_species.1.ic.vx0 = 0.0
kinetic_species.1.ic.vy0 = 0.0
kinetic_species.1.ic.A = 0.0
kinetic_species.1.ic.B = 0.0
kinetic_species.1.ic.C = 0.0
kinetic_species.1.ic.kx1 = 0.0
kinetic_species.1.ic.ky1 = 0.0
kinetic_species.1.ic.kx2 = 0.0
#
#
kinetic_species.2.name = "ion" 
kinetic_species.2.velocity_limits = $vIxmin $vIxmax $vIymin $vIymax 
kinetic_species.2.Nv = $Nvx $Nvy
kinetic_species.2.mass = $IMass 
kinetic_species.2.charge = +1.0
kinetic_species.2.tz.name = "TwoSpecies_IonTrigTZSource"
kinetic_species.2.tz.amp = 0.1 
kinetic_species.2.tz.electron_mass = $EMass
kinetic_species.2.tz.ion_mass = $IMass
#
#
kinetic_species.2.ic.name = "Perturbed Maxwellian"
kinetic_species.2.ic.tx = $IT 
kinetic_species.2.ic.ty = $IT 
kinetic_species.2.ic.vx0 = 0.0
kinetic_species.2.ic.vy0 = 0.0
kinetic_species.2.ic.A = 0.0
kinetic_species.2.ic.B = 0.0
kinetic_species.2.ic.C = 0.0
kinetic_species.2.ic.kx1 = 0.0
kinetic_species.2.ic.ky1 = 0.0
kinetic_species.2.ic.kx2 = 0.0
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
