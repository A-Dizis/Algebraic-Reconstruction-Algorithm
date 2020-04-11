# Algebraic-Reconstruction-Algorithm
This is an application of ART algorithm in 2-D using multithreading and a new approach on the calculation of the projection matrices.

For the user to get a grasp of the internal workings of this algorithm, it is highly recommended to, initially skip through the paper "A tutorial on ART (Algebraic Reconstruction Technique)" by R. Gordon.

In this paper the projection matrices are not calculated with accuracy, as they are given bit values 0 or 1 in each position. Here we calculate the values with great accuracy, using lattice points for the estimation of the real area of each position as it can be seen in the image. 

To provide data (a sinogram) to the algorithm, and octave file output format was used. This is nothing more than values set in x - y axis format. Values should be normalized to one, as they represent the intensity of the beam received  after passing though the object that is being reconstucted. The format can be seen in the example file sinogram.dat.

Description of functions:

- get_info: Gets info about sinogram(available angles(angles), last angle position in degrees(thetamax), available projections(projections) and N the pixel width of the reconstructed image NxN).

- get_sinogram: Reads sinogram.dat file and loads it in memory.

- show_sinogram: Helper function, will print the projection matrix for error checking.

- get_density: Asks user for the density of lattice's points. The greater, the better the quality gets and slower the algorithm is being executed.

- create_points_position: Creates a struct that holds the point's x,y positions.

- write_points: Helper function, will print position of points created.

- get_thread_number: Asks the user the number of the threads that will be created. By rule, you should use as many as your CPU has available.

- initialize_image: Sets the pixels of the created image to 0.

- make_projection: Creates the projection matrices for the current angle of the calculation.

- normalize_projection: Normalizes the areas of the projections matrices to 1.

- print_projection_matrices: Helper function, will print all projection matrices for the angle to a file.

- contribution: Calculates the new corrections and refreshes the reconstructed image.

- ask_for_repetition: Asks for the continuation o f calculations by pass throught the sinogram data once more. By this time the algorithm will have already made a full pass through the sinogram data.

- normalize_image: If the calculations are concluded, this function will normalize the final image.

- output_image: Will output the final image to a file.

EXAMPLE: You can try reconstruct the example sinogram sinogram.dat into an image. The paramenters here are the there are 360 angles available, and 224 projections. Keep the file in the same directory with the executable for it to load. If you executed the reconstruction correct a final image file should appear under the name output_image.dat. This file contains all the appropriate data in order to be transformed back again into a bmp or png image using octave or appropriate software.

IDEAS FOR EXTRA FEATURES: 
- It would be a nice feature to add a live previewer of the created image using OpenGL.
- Autorepetition or autocancelation can be made using a measure of change between consecutive steps of correction.
