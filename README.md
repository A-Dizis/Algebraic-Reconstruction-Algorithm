# Algebraic-Reconstruction-Algorithm
This is an application of ART algorithm in 2-D using multithreading and a new approach on the calculation of the projection matrices.

For the user to get a grasp of the internal workings of this algorithm, it is highly recommended to, initially skip through the paper "A tutorial on ART (Algebraic Reconstruction Technique)" by R. Gordon.

In this paper the projection matrices are not calculated with accuracy, as they are given bit values 0 or 1 in each position. Here we calculate the values with great accuracy, using lattice points for the estimation of the real area of each position as it can be seen in the images. 

To provide data (a sinogram) to the algorithm, and octave file output format was used. This is nothing more than values set in x - y axis format. Values should be normalized to one, as they represent the intensity of the beam received  after passing though the object that is being reconstucted. The format can be seen in the file.

Description of functions:

 get info: Gets info about sinogram(available angles(angles), last angle position in degrees(thetamax), available projections(projections) and N the pixel width of the reconstructed image NxN).

 get sinogram: Reads sinogram.dat file and loads it in memory.

 show sinogram: Helper function, will print the projection matrix for error checking.

 get density: Asks user for the density of lattice's points. The greater, the better the quality gets and slower the algorithm is being executed.

 create points position: Creates a struct that holds the point's x,y positions.

 write points: Helper function, will print position of points created.

 get thread number: Asks the user the number of the threads that will be created. By rule, you should use as many as your CPU has available.

 initialize image: Sets the pixels of the created image to 0.

 make projection: Creates the projection matrices for the current angle of the calculation.

 normalize projection: Normalizes the areas of the projections matrices to 1.

 print projection matrices: Helper function, will print all projection matrices for the angle to a file.

 contribution: Calculates the new corrections and refreshes the reconstructed image.

 ask for repetition: Asks for the continuation o f calculations by pass throught the sinogram data once more. By this time the algorithm will have already made a full pass through the sinogram data.

 normalize image: If the calculations are concluded, this function will normalize the final image.

 output image: Will output the final image to a file.
