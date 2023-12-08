## Code description

The code solves the 2D diffusion and coagulation equation. Compilation is done using gcc. The makefile should produce an executable called **tet.exe**. After the code is executed, pictures (in .ppm format) will be stored in **run_folder/imgs**.

Command-line arguments: ``./tet.exe 0.3 0.6 0 0.2 30 9``: ``0.3`` - x-coordinate of the monomer source (in %), ``0.6`` - y-coordinate of the monomer source (in %), ``0`` - coagulation on/off (0 -- off, 1 -- on), ``0.2`` - velocity value (in m/s), ``30`` - velocity direction (in degrees), ``9`` - the source intensity (in $cm^{-3}$ ).

The source of monomers can be specified at any point of the two-dimensional plane.

What can be modified:
* boundary conditions;
* advection;
* computational domain;
* image size;
* aggregation kernel;
* equation coefficients;
* source coordinates;
* discretization steps along all the axes.
