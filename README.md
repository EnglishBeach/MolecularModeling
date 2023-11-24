# Diffusion molecular model
It is part of my diploma work. There are files for the molecular modeling system oligocarbonate methacrylate OKM-2 and butan-1-ol.
This simulation helps refine the diffusion coefficient and improve the prediction force of my diffusion-kinetic model.
See the general repo (it is in process): https://github.com/EnglishBeach/Diplom-work


 * Using GROMACS for preparation and simulation and VMD for visualization and analysis
 * Using oplsaa forcefield and tip4p for water forcefield
 * Using ATB to generate .itp and .gro files of moleculas


# Simulating
Right now, there is only one type of simulation: OKM-butanol\Simple_MD.
But in the future, this list will be increased.


 1. Install GROMACS and VMD.
Install GROMACS: https://manual.gromacs.org/documentation/current/install-guide/index.html
Install VMD:
https://mejk.github.io/moldy/chapters/visualization/vmd-install.html


 2. Get a pdb file or generate .gro .itp files for molecules on web resources (Charmm-GUI, ATB, etc.). Keep attention to forcefield names.


 3. Generate a box from a pdb or .gro file:
    ```gmx pdb2gmx -f butane.pdb -o butane.gro```
Change box size:
    ```gmx editconf -f butane.gro -box 4 4 4 -o butane.gro```


 4. Create a .top file.
 It has restrictions: see documentation to do it correctly. This file is the box description (parameters, atomtypes and so on)


 5. Insert molecules into 1 box and correct the .top file by hand:
    ```gmx insert-molecules -f box.gro -ci butane.gro -nmol 10 -p config.top```


 6. Compile the simulation and run it.
Run it from the fully prepared system directory.
    ```gmx grompp -p config.top -c box.gro -f params.mdp```
    ```gmx mdrun -v``
After waiting, you will have a trajectory file, but you need to correct it.


 7. Correcting:
    ```gmx trjconv traj_comp.xtc -o trajout.xtc -pbc mol``



# In process
 * Using python for simplifying solutions and preparations
 * Automated calculation for different conditions