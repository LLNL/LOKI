$C_frac = 5.0/3.0;
$mp_over_me = 1836.;
#electrons
$Nvx_e = 24;
$Nvy_e = 16;
# He
$Z_he   = 2.;
$A_he   = 4.;
$T_he   = 1.;
$m_he   = $A_he*$mp_over_me;
$vth_he = sqrt($T_he/$m_he);
$Nvx_he = 24;
$Nvy_he = 16;
# C
$Z_c    = 6.;
$A_c    = 12.;
$T_c    = 1.;
$m_c    = $A_c*$mp_over_me;
$vth_c  = sqrt($T_c/$m_c);
$Nvx_c  = 24;
$Nvy_c  = 16;
#
$Nx = 128;
$Ny = 7;
$xa = -62.5;
$xb =  62.5;
$dx =  ($xb-$xa)/$Nx;
$ya = -0.5*$Ny*$dx;
$yb =  0.5*$Ny*$dx;
#
$vxmin_e  = -7.5;
$vxmax_e  =  7.5;
$vymin_e  = -7.5;
$vymax_e  =  7.5;
$vxmin_he = -7.5*$vth_he;
$vxmax_he =  7.5*$vth_he;
$vymin_he = -7.5*$vth_he;
$vymax_he =  7.5*$vth_he;
$vxmin_c  = -7.5*$vth_c;
$vxmax_c  =  7.5*$vth_c;
$vymin_c  = -7.5*$vth_c;
$vymax_c  =  7.5*$vth_c;
#
#
# now change the default quantities as desired
verbosity = 1
#
# Configuration space Domain
#
domain_limits = $xa $xb $ya $yb 
N = $Nx $Ny
periodic_dir = true true
#
# Time step controls
#
cfl = 0.95
final_time = 5.
save_times = 1.
sequence_write_times = .1
max_step = 1000000
spatial_solution_order = 6
temporal_solution_order = 6
#
# Plasma definition
#
number_of_species = 3
#
kinetic_species.1.name = "electron"
kinetic_species.1.velocity_limits = $vxmin_e $vxmax_e $vymin_e $vymax_e
kinetic_species.1.Nv = $Nvx_e $Nvy_e
kinetic_species.1.mass = 1.0
kinetic_species.1.charge = -1.0
kinetic_species.1.ic.name = "Interpenetrating Stream"
kinetic_species.1.ic.syntax = "half plane"
kinetic_species.1.ic.tl = 1.0
kinetic_species.1.ic.tt = 1.0
kinetic_species.1.ic.theta = 0.0
kinetic_species.1.ic.d = 31.25
kinetic_species.1.ic.beta = 0.768
kinetic_species.1.ic.two_sided = true
kinetic_species.1.ic.floor = 0.05
kinetic_species.1.ic.frac = 10.0
kinetic_species.1.ic.frac2 = 10.0
#
kinetic_species.2.name = "He"
kinetic_species.2.velocity_limits = $vxmin_he $vxmax_he $vymin_he $vymax_he
kinetic_species.2.Nv = $Nvx_he $Nvy_he
kinetic_species.2.mass = $m_he
kinetic_species.2.charge = $Z_he
kinetic_species.2.ic.name = "Interpenetrating Stream"
kinetic_species.2.ic.syntax = "half plane"
kinetic_species.2.ic.tl = $T_he
kinetic_species.2.ic.tt = $T_he
kinetic_species.2.ic.theta = 0.0
kinetic_species.2.ic.d = 31.25
kinetic_species.2.ic.beta = -0.768
kinetic_species.2.ic.centered = true
kinetic_species.2.ic.floor = 0.0
kinetic_species.2.ic.frac = 0.025
#
kinetic_species.3.name = "C"
kinetic_species.3.velocity_limits = $vxmin_c $vxmax_c $vymin_c $vymax_c
kinetic_species.3.Nv = $Nvx_c $Nvy_c
kinetic_species.3.mass = $m_c
kinetic_species.3.charge = $Z_c
kinetic_species.3.ic.name = "Interpenetrating Stream"
kinetic_species.3.ic.syntax = "half plane"
kinetic_species.3.ic.tl = $T_c
kinetic_species.3.ic.tt = $T_c
kinetic_species.3.ic.theta = 0.0
kinetic_species.3.ic.d = 31.25
kinetic_species.3.ic.beta = 0.768
kinetic_species.3.ic.two_sided = true
kinetic_species.3.ic.floor = 0.0
kinetic_species.3.ic.frac = $C_frac
kinetic_species.3.ic.frac2 = $C_frac
#
# Diagnostic Controls
#
number_of_probes = 0
#
save_data = true
#
start_from_restart = false
restart.time_interval = 5
restart.write_directory = "InterpenetratingStreams"
