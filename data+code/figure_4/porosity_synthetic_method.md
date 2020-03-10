# Porosity Synthetic Code

1. We begin with the *ProportionPacking.mat* output, as output from Packing3D (Mollon amd Zhao, 2014). For reference purposes, we've included the *Packing3D_Main_Program.mat*, which contains the variables required to reconstruct a granular sample that has similar properties to the one used in this work. NOTE: Some files for this routine are too large for many online repositories, so they are not included here. They will be shared upon request. 

2. We run *points_to_csv.mat* twice, once for each of the decimated and undecimated datasets:

    ```MATLAB
        % Decimated volume
        points_to_csv(trains, VolumesGrains, '\Grains_Decimated', false, 1, true, 0.15);
        % Undecimated volume
        points_to_csv(Grains, VolumesGrains, '\All_grains', false, 1, false, 0.15);
    ```

3. Next, we run *generate_grains.gh*, a Grasshopper file, in Rhino 3D. Grasshopper is a plug-in for Rhino 3D that enables parametric operations using a graphical user interface. For performance purposes, we recommend disabling the Grasshopper solver (Solution --> Disable Solver in the menu bar) before opening *generate_grains.gh*. The Grasshopper code parses a directory (i.e., either Grains_Decimated or All_grains) and builds a 3D mesh for each grain. 

1. Here, we use the volume of grains modelled in *porosity_grains_model.3dm* (a Rhino 3D file). In addition to the grains, the file also includes layers corresponding to the centroid of the volume, as well as a cross section through the volume (made at z=0).
2. Within Rhino 3D, we run *grain_sectioning_porosity_calculation.gh*, which is a Grasshopper file. Grasshopper is a plug-in for Rhino 3D that enables parametric operations using a graphical user interface. For performance purposes, we recommend disabling the Grasshopper solver (Solution --> Disable Solver in the menu bar) before opening *grain_sectioning_porosity_calculation.gh*.
3. 