from pathlib import Path

from base import MolNames


def em(configs_path: Path):
    config = """
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 500.0         ; Stop minimization when the maximum force < 1000.0 J/mol/nm
emstep      = 0.001         ; Minimization step size fs
nsteps      = 5000           ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
rcoulomb        = 1.2       ; Short-range electrostatic cut-off
rvdw            = 1.2       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions

; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT

; Dispersion correction
DispCorr        = EnerPres  ; account for cut-off vdW scheme
"""
    with open(configs_path / 'em.mdp', 'w') as config_file:
        config_file.write(config)
    return configs_path / 'em.mdp'


def nvt(configs_path: Path, compounds: list[MolNames], T: int):
    tau_ts = ' '.join(['0.1' for compound in compounds])
    groups = ' '.join([compound.value for compound in compounds])
    Ts = ' '.join([str(T) for compound in compounds])

    config = f"""
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 10000     ; 1 * 1000 * 10 = 10 ps
dt                      = 0.001     ; 1 fs

; Output control
nstxout                 = 0        ; suppress bulky .trr file by specifying
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 1000      ; save energies every 1.0 ps
nstlog                  = 1000      ; update log file every 1.0 ps
nstxout-compressed      = 1000      ; save compressed coordinates every 2.0 ps
compressed-x-grps       = System    ; save the whole system

; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
nstlist                 = 20        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = {groups}       ; two coupling groups - more accurate
tau_t                   = {tau_ts}       ; time constant, in ps
ref_t                   = {Ts}           ; reference temperature, one for each group, in K

; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; group(s) for center of mass motion removal
nstcomm                  = 10
comm-grps                = {groups}

; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = {T}       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
"""
    with open(configs_path / 'nvt.mdp', 'w') as config_file:
        config_file.write(config)
    return configs_path / 'nvt.mdp'


def npt(configs_path: Path, compounds: list[MolNames], T: int):
    tau_ts = ' '.join(['0.1' for compound in compounds])
    groups = ' '.join([compound.value for compound in compounds])
    Ts = ' '.join([str(T) for compound in compounds])

    config = f"""
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 10000     ; 1 * 1000 * 10 = 10 ps
dt                      = 0.001     ; 1 fs

; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 1000      ; save energies every 1.0 ps
nstlog                  = 1000      ; update log file every 1.0 ps
nstxout-compressed      = 1000      ; save compressed coordinates every 1.0 ps
compressed-x-grps       = System    ; save the whole system

; Bond parameters
continuation            = yes       ; Restarting after NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
nstlist                 = 20        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = {groups}       ; two coupling groups - more accurate
tau_t                   = {tau_ts}       ; time constant, in ps
ref_t                   = {Ts}           ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; group(s) for center of mass motion removal
nstcomm                  = 10
comm-grps                = {groups}

; Velocity generation
gen_vel                 = no        ; Velocity generation is off
"""
    with open(configs_path / 'npt.mdp', 'w') as config_file:
        config_file.write(config)
    return configs_path / 'npt.mdp'


def md(configs_path: Path, compounds: list[MolNames], T: int):
    tau_ts = ' '.join(['0.1' for compound in compounds])
    groups = ' '.join([compound.value for compound in compounds])
    Ts = ' '.join([str(T) for compound in compounds])
    config = f"""
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 2000
dt                      = 0.002     ; 2 fs

; Output control
nstxout                 = 0      ; suppress bulky .trr file by specifying
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 1000      ; save energies every 2.0 ps
nstlog                  = 1000      ; update log file every 2.0 ps
nstxout-compressed      = 1000      ; save compressed coordinates every 2.0 ps
compressed-x-grps       = System    ; save the whole system

; Bond parameters
continuation            = yes       ; Restarting after NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
nstlist                 = 20        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = {groups}       ; two coupling groups - more accurate
tau_t                   = {tau_ts}       ; time constant, in ps
ref_t                   = {Ts}           ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1

; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC

; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme

; Velocity generation
gen_vel                 = no        ; Velocity generation is off

; group(s) for center of mass motion removal
nstcomm                  = 10
comm-grps                = {groups}
"""
    with open(configs_path / 'md.mdp', 'w') as config_file:
        config_file.write(config)
    return configs_path / 'md.mdp'
