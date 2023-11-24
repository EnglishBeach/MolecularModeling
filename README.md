# Diffusion molecular model
It is part of my diplom work. There are files for molecular modeling system oligocarbonate methacrylate OKM-2 and butan-1-ol.
This simulation helps refine diffusion coefficient and improve a prediction force of my diffusion-kinetic model.
See the general repo (it is process): https://github.com/EnglishBeach/Diplom-work

 * Using GROMACS for preparation and simulation, VMD for visualisation and analysis
 * Using oplsaa forcefield and tip4p for water forcefiled
 * Using ATB for generating .itp and .gro files of moleculas

# Simulating
Rihgt now there is only one type simulation - OKM-butanol\Simple_MD
But in future this list will be increased

 1. Install GROMACS, VMD
Install GROMACS: https://manual.gromacs.org/documentation/current/install-guide/index.html
Install VMD: https://mejk.github.io/moldy/chapters/visualization/vmd-install.html

 2. Get pdb file or generate .gro .itp files for molecules on web resources (Charmm-GUI,ATB etc. ) Keep attention with forcefields names

 3. Generate box from pdb or .gro file:
    ```gmx pdb2gmx -f butane.pdb -o butane.gro```
Change box size:
    ```gmx editconf -f butane.gro -box 4 4 4 -o butane.gro```

 4. Create .top file
 It has restrictions see documentations to do it correctly. This file is the box description (parametrs, atomtymes and so on)

 5. Insert molecules into 1 box and correct .top file by hands:
    ```gmx insert-molecules -f box.gro -ci butane.gro -nmol 10 -p config.top```

 6. Compile simulation and run it:
Run it from full done system directory.
    ```gmx grompp -p config.top -c box.gro -f params.mdp```
    ```gmx mdrun -v``
After waiting you will have trajectory file, but need to correct it

 7. Correcting:
    ```gmx trjconv traj_comp.xtc -o trajout.xtc -pbc mol``


# In process
 * Using python for simplify solves and preparations
 * Automatization calculation for different conditions
