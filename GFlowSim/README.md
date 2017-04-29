Options for the 'driver' program

SIMULATION OPTIONS

-number	    	The number of particles to start the program with. By default, the simulation will have enough particles to fill the simulation volume by 50% (phi=0.5). This option instead starts the simulation with exactly the specified number of particles.

-width	    	The width of the simulation domain  (default 4).

-height	    	The height of the simulation domain (default 4).

-radius     	The radius (or max radius if there is dispersion) of the particle sigmas (default 0.05).

-bR         	The radius of the intruding ball (for loaded buoyancy runs) (default 0.2) .

-density    	The density of the intruding ball (for loaded buoyancy runs) (default 100).

-velocity   	The target velocity of particles in normal simulations, or the velocity of the intruding ball for loaded buoyancy runs (default 0.25).

-dispersion 	The radius dispersion for the particles. Denotes the max % deviation, e.g. if radius=1 and dispersion=0.15, the radii of the particles will be uniformly distributed over the interval [0.85, 1]. (Default 0.)

-temperature	The temperature of the system, uses brownian perturbations and viscous drag. (Default 0.)

-gravity       	The (y component of the) gravity of the system. (Default is -1.)

-drag 		The drag coefficient of the system (default 0.).

-coeff 		The coefficient of friction of the particles. The default value is whatever default value is set to in DefaultConstants.h.

-time 		Simulation run time. (Default 1 'second.')

-start 		When to start recording data (if applicable). (Default 0.)

-epsilon 	Time step. (Default 1e-4.)

-phi            The packing fraction, the simulation (when doing a normal run) calculates how many particles to use based off this and the simulation volume.

-skinDepth      The skin depth to use in creating neighbor lists.

-interaction	What kind of interaction the particles should have. 0 - Hard sphere, 1 - Lennard-Jones sphere, 2 - Hard triangle. Default is hard spheres (0).

-LJ 		Sets the interaction to be Lennard-Jones (same as -interaction=1).

-Tri 		Sets the interaction to be hard triangles (same as -interaction=2).

-interact 	Whether the particles should interact with one another. (Default true.)

-seedRand 	Whether the random number generators should be seeded before program execution (default true).

-quiet     	Level of quiet. -1 - Normal, 1 - Full quiet, 2 - Quiet if label!=0

-lattice   	What kind of lattice initialization to use for creating buoyancy simulations. 0 - None, use find packed solution, 1 - Hexagonal lattice, 2 - rectangular lattice. (Default 0)

ANIMATION OPTIONS

-mode   	The animation mode (default 0).

-noprint 	If true, don't print a summary (default false). 

-animate 	Whether to create an animation of position (default false).

-center   	Whether to center the animation window around the intruding ball (default false). 
		
-snapshot 	Whether to take a snapshot of the simulation at the end (default false).

-special 	Animate, coloring particles by what processor they are on (default false).

-forces  	Animate, coloring by the force or pressure felt by the particles.

-forceChoice 	Which forces to include in the force animation.

-typeChoice 	Whether to use particles and walls (0), particles (1), or walls (2) when making force animation.

-bubbles 	Whether to record the number and volume of "bubbles." (Default false.)

-visBubbles 	Create an animation of the "bubbles." (Default false.)
		
-writeFields 	Whether to write the resource and waste field data to files (for printing) (default false). Only meaningful when doing a 'bacteria' run.

-writeFitness	Whether to write the fitness "field" to files (for printing) (default false).

-writeAnimation If true, we write animation data to files rather then printing it out (default true).

-bulk  		Create a Bulk Animation (default false).

STAT FUNCTION OPTIONS (ALL DEFAULT FALSE)

-omega		Average Omega.

-KE      	Average kinetic energy.

-KEX     	Average x-direction kinetic energy.

-KEY    	Average y-direction kinetic energy.

-LKE     	Average linear kinetic energy.

-RKE     	Average rotational kinetic energy.

-cluster 	Average clustering.

-triAlign 	Average triangle alignment (don't yet have a good way to measure this).

-trackHeight  	Height of the large ball.

-trackX  	X - position of the large ball.

-GPE     	Gravitational potential energy of the system

-maxV    	Max velocity of any particle in the system

-maxVy   	Max y-velocity of any particle in the system.

-num     	The number of particles in the simulation

STAT PLOTTING OPTIONS (ALL DEFAULT FALSE)

-velDist        Plot the velocity frequency distribution.

-pressurePlot 	Plot the pressure vs. y distribution

MISCELANEOUS

-novid		Don't print video making commands

-printSectors 	Print sector occupancy

-fps 		Frames per second (data recording frequency, default 15)

-label		A label for data names

-scale		The scale for making videos

