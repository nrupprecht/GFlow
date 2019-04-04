# GFlow v 4.0 Configuration files.

## Getting Started.

Configurations files are used to define and set up simulations. A configurations files is essentially a text file with some basic commands, no particular file ending is required. I generally use .txt files for convenience. Configuration files can be used to define interparticle forces, the type of integrator to use, to fill areas with particles, to add individual particles, add modifiers to particles, and related tasks. The configuration file does not specify how long a simulation should run - this should be specified separately from the command line.

Configuration files are parsed by the class *Parser* class, located in *src/utility/parser.hpp*, and the parse tree is converted into a simulation by the class *FileParseCreator*, located in *src/creators/fileparsecreator.hpp*.

### Basic Usage.

Configuration files are loaded, parsed, and converted into simulations by the program *bin/driver*. A basic usage of this program would be

    ./bin/driver -load=[configuration file] -time=[time]
    
where the configuration file should be the path to the file relative to the working directory, and the amount of time should be in seconds. 
Another common example would be a run where kinetic energy data and particle trajectory data is recorded, to make a video of a specified length:

    ./bin/driver -load=[configuration file] -time=[time] -KE -animate -videoLength=[v time]

### Basic setup.

Configuration files contain lists of structures of the form

    [Head] : [param-1]=[arg-1], ... , [param-k]=[arg-k] { 
      ... [body] ...
    }

where the body also contains the same type of structures. Arguments of parameters are often unneccessary, e.g.

    [Head] : [param-1], [param-2], ... , [param-k]=[arg-k] { 
      ... [body] ...
    }

as are parameters or bodies, or, more occassionally, heads.

The parser collects all top level heads, and the file parse creator looks for specific heading in the order that *it* chooses to - consequently, there is no reason to write a configuration file in any particular order.


### Debugging.

If an invalid header is used (i.e. one that is not recognized by the file parse creator), a message will be printed to the screen of the form " -- Heading: [Head]".
If a fatal error occurs while parsing or constructing the simulation, an exception will be thrown, typically a *BadStructure* error, which will result in an error message being printed, e.g. "Caught Bad Structure error: Unrecognized force grid parameter." Following this will be a "Build Message" that summarizes what parts of the simulation construction were able to occur before the error was thrown. Finally, a file "ParseTrace.txt" will be written to your working directory that shows how the configuration file was parsed. You should use all this to discover what is wrong with your configuration file.

## Configuration file setup.

This section describes the various options available for creating simulations.

### Variables.

Variables can be defined for use in configuration files. Anywhere a parameter or value is expected, a variable can be used. A variable may be either a string or a number, the file parse creator will interpret it to be whatever it should be in the context in which it is used. Since variables can represent strings or numerical values, mathematical expressions will not be parsed for variables. Variables should be defined as follows:

    Var: [Name]=[Default value]
    
Any number of variable statements can be included in the configuration file. Variables must be in the the top level of the configuration file, i.e. not in a body.

The other useful thing about variables is that they can have their values set from the command line by adding the flag *-[Name]=[Value]*. For example, if you have defined a variable to specify the temperature of a simulation,

    Var: temperature=0.25
    
you can set the Temperature variable to be 1.0 from the command line by running the driver,

    ./bin/driver -load=my-configuration.txt -temperature=1.0
    
### Values.

Values are similar for variables, but cannot be set from the command line, and can be functions of variables and previously defines values. Another difference is a value **must** be a numerical value (it will evaluate as a float type). Therefore, any variables that the value depends on must evaluate as numbers. Arithmetic expressions with '+' '-' '*' '/' '(' and ')' may be used. For example, a value called "length" can be defined as

    Value: lenth = (width + 1)/2
  
and used to define subsequent values, e.g.

    Value: volume = length*length

### Comments.
Comments follow the C/C++ convention. Single line comments are started with a // and last until the end of the line. Multiline comments are of the form 

    /* ... Comment goes here ... */ 
    
and are only terminated by the closing */, not by endline characters. Multiline comments do not have to stretch to multiple lines, anything withing the starting and ending terminators will be ignored by the parser.

### Essential parts of a configuration file.

There are two top level heads that must be included in each configuration file: **Bounds**, and **Interactions**. Bounds specifies the simulation bounds, which currently must be rectangular. Bounds should be specified as follows,

Full specification: the bound's min/max in each dimension is explicitly given. The number of pairs of bounds must match the dimensionality of the simulation.

    Bounds: {
      : [min x], [max x]
      : [min y], [max y]
      ...
      : [min d_n], [max w]
    }
    
Box specification: bounds with min / max of 0, [length] in each dimension is created,

    Bounds: Box=[length]

There are several options as to how to specify interactions. The usual way is to specify each pair of particles types, and what interaction should exist between them, e.g.

    Force-grid: {
      : 0, 0, HardSphere
      : 0, 1, LennardJones
      : 1, 1, None
    }

Forces are by default reflexive, that is, it is unneccessary to also say *: 1, 0, LennardJones*. There is a way to make forces non-reflexive, but you shouldn't do that, so I'm not going to say what it is. The force options currently are summarized below.

| Token        | Description                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------------|
| HardSphere   | Hard sphere force - spring like repulsion.                                                            |
| HardSphereDs | Hard sphere force with dissipation                                                                    |
| LennardJones | Lennard Jones 6-12 force.                                                                             |
| Buckingham   | Buckingham force.                                                                                     |
| Detector     | Not a force, but an interaction that stops the simulation  if a particle with enough kinetic hits it. |
| None         | Explicit indication that there is no force.                                                           |

A second way to specify the force is to make all particles of all types interact with the same force. To do this, use

    Force-grid: [Interaction-token]

The final way to allocate forces is by the creation of a random force network. To choose this option, simply include

    Force-grid: Random
    
This is useful in the simulation of many-type-physics simulations. The file parse creator will choose the forces randomly to either be **HardSphere**, or **LennardJones**.

### Other important configuration data.

The following data is not *essential* to include in a configuration file (i.e. if the information is not included, default values will be assumed, and no exceptions will be thrown), but is often useful, or at least a good idea to include explicitly.

**Dimensions**: The number of dimensions the simulation should run in. In theory, any positive integer is valid, but high dimensional simulation will be very slow, and I have not tested them extensively. The number of dimensions is two by default, and can be specified as

    Dimensions: [D]

**Boundary**: This is how boundary conditions are specified. Boundary conditions are **Wrap** (harmonic) by default. The options are **Wrap** (harmonic), **Reflect** (reflecting), or **Repulse** (boundary repells particles via a hard sphere force). Similarly to **Bounds**, boundary conditions are specified by

    Boundary: {
      : [BC-x]
      : [BC-y]
      ...
      : [BC-w]
    }
    
  All the boundary conditions can be simultaneously be set to the same type by:
    
    Boundary: [BC]
  
**NTypes**: The (maximum) number of types of particles in the simulation. By default, there is only one type (type 0).

    NTypes: [ntypes]
    
**Integrator**: How to integrate the motion of the particles. The default integrator is the velocity verlet integrator. The valid types of integrators are **VelocityVerlet**, **OverdampedIntegrator**, **OverdampedLangevin**, and **LangevinIntegrator**. 

All integrators can have the following options specified in their bodies,

    Integrator: [Integrator type] {
      Delay: [delay] // Set the delay between checking if the timestep should be adjusted. 
      MinDT: [minDT] // Minimum allowable timestep.
      MaxDT: [maxDT] // Maximum allowable timestep.
      Adjust: [0/1]  // Whether the timestep should be adjusted on the fly. True by default.
      UseV: [0/1]    // Whether the maximum velocity should be used to determine the timestep. True by default.
      UseA: [0/1]    // Whether the maximum acceleration should be used to determine the timestep. False by default.
    }
    
Additionally, both Langevin type integrators can take the following options,

    Temperature: [temperature]
    Viscosity:   [viscosity]
    
which specify the simulation's temperature and viscosity.
  
### Particles and Modifiers.

Of course, to have a useful simulation, there should be particles. There are currently two ways to add particles: via a **Fill-area**, or by adding individual particles via **Particle**.

**Particle**: Adds a single particle to the simulation. You must specify the position and cutoff radius of the particle - other data entries are optional (have default values). Some modifiers can also be attached to the particle.

    Particle: {
      Position: [x], [y], ... , [w] // The number of coordinates must match the number of dimensions.
      Sigma: [r]
      
      // The velocity is optional, and zero by default.
      Velocity: [vx], [vy], ... , [vw] 
      
      // Specify the mass - pick *one* of the following three ways. Or none - by default, density is 1.
      Mass: [m] // To specify the mass of the particle.
      <or>
      Mass: Density=[density] // To specify the density of the particle.
      <or>
      Mass: inf // For an infinitely massive particle.
      
      Type: [type] // Type is 0 by default.
      
      // Any number of modifiers can be added. Try not to add irreconcilable modifiers.
      Modifier: CV // Constant velocity modifier.
      <or>
      Modifier: CV-D=[distance] // Move at constant velocity for a set distance.
    }
    
**Modifier**: Adds a modifier to the simulation. Currently, the only modifier is **WindTunnel**, which creates the wind tunnel modifier. This accelerates particles at the X edges of the simulation in the +X direction. When using this, you should be sure that the X boundary condition is **Wrap**.

    Modifier: WindTunnel=[Driving velocity]


**Attraction**: Really, this should be moved to be under modifier. This modifier is a constant accelleration attraction towards the center of the simulation. Specify by

    Attraction: [acceleration]

### Particle Filling.

This is the most useful, flexible, and complex way to create particles. Basically, some volume (area, hypervolume) has particles added to it. There are several types of fills. 

The first is **Fill: Area**. This fills an area with particles by picking positions uniformly. The types, velocities, radii, etc. for the particles can all be specified in different ways. The bounds of the area to be filled *must* be specified, as well as information on how many particles should be added. There should be no parameters. The body specifies how to create particles.

**Bounds**: The bounds can either be specified in the same format as the simulation bounds, or by 

    Bounds: Full

which sets the bounds to be the full simulation bounds. 

Additionally, Spherical bounds can be set by

    Bounds: Sphere {
      : [x], [y], ... , [w]
      : [radius]
    }

**Template**: The way to define particles is though particle templates, which specify how a particle should be generated. Any number of particle templates may be specified. Particle templates are their own block, and have the template name as their only parameter. The mass, radius, and type for the template must all be specified, either as a deterministic number, or using a random generator. 

If a specific value is requested, it should be specified by e.g.

    Mass: Density=1

There are several options for random values. In all of the below, *Head* can be either **Mass**, **Radius**, or **Type**. There is a discrete random generator:

    [Head] : Random {
      P: [value-1], [prob-1]
      P: [value-2], [prob-2]
      ...
      P: [value-k], [prob-k]
    }
    
so *value-j* occurs with (relative) probability *prob-j*. In other words, there is no need to normalize your probabilities.

The **Uniform** generator is defined by

    [Head] : Uniform {
      Min: [min]
      Max: [max]
    }
    
and generates random numbers uniformly between *min* and *max*. *Min* and *Max* can come in either order.

The **Normal** generator is defined by

    [Head] : Uniform {
      Ave: [average]
      Var: [variance]
    }

*Ave* and *Var* can come in either order.

Finally, there is a special option for type specification: **Equiprobable**. This will result in particles of every possible type, from 1 to *ntypes* being generated with equal probability. This is useful for many-type-physics simulations. Explicitly, you say

    Type: Equiprobable

**Number**: Particle number information is based on particle templates. Either a specific number of particles can be specified, or the simulation is filled to a specified volume density. The (relative) probability of particle templates should be included in the body of the **Number** header:

    Number: [number] {
    <or>
    Number: Phi=[density] {
      [Template A] : [p-A]
      ...
      [Template Z] : [p-Z]
    }
    
It is not neccessary to normalize the probabilities for the particle templates. If, for example, you want each particle template to occur with equal probability, you can just use "1" for each p-value.
    
**Velocity**: There are currently not many ways to specify the initial velocity of the particles. The default is to assign particles random normal velocities. To set all the particles to move at a specific velocity, use one of the following two headers

    Velocity: [vx], [vy], ... , [vw] // Specifies a velocity vector.
    <or>
    Velocity: Zero // All particles will start with zero velocity.

**Attraction**: This option adds a constant accelaration that attracts all particles towards the center of the simulation domain. This can be useful if there will be such an attraction included in the full simulation, as it will allow particles to rest naturally in place. This can be specified by:

    Attraction: [acceleration]
    
**Exclude**: Defines excluded regions within which particles will (probably) not be placed. If too many failures to pick a random point outside the exclude bounds occur, then particles may be placed within the exclude region. Also, particle relaxation could cause particles to be pushed into the excluded region. Regions may be defined via the same construction as the Fill: Area bounds (except for **Bounds: Full**).


The second type of filling is **Fill: Polymer**. This creates one or more random polymer chains within the simulation bounds. Each polymer is created via a random walk out of two kinds of particles, (usually large) "primary" particles, and (usually small) "chain" particles. There are harmonic bonds between adjacent particles in the chain. There are several parameters that can all be set with the construction

    [Heading]: [Value]

**Number**: The number of polymers to create.

**Length**: The length of each polymer.

**R**: The radius of the primary particles.

**r**: The radius of the chain particles.

**Phi**: The target density of the polymer, Phi = 2 * R / length.

**IdP**: The type of the primary particles.

**IdC**: The type of the chain particles.

**Correlation**: If true, a group correlation object is added to monitor correlations between the first polymer created and the rest of the particles in the simulation.

**Parallel**: If this is set to true (non-zero), then two straight, parallel polymers are constructed.

**H**: Only for use with **Parallel** (though no errors occur otherwise - it's just that nothing happens). This specifies how far apart (as a fraction of **R**) the polymers will be from one another.

In addition, a LineEntropicForce data object is added to the simulation, attached to the first polymer created.

### Particle Relaxation.

Since particle positions are generally chosen randomly, e.g. by a fill area command, particles likely overlap with one another. To rectify this, the simulation setup is run for some amount of time with an overdamped integrator to force particles to exclude one another. The default relaxation time is 0.5, but the amount of time can be specified,

    Relax: [time]
    
*After* the relaxation step, particles will have their velocities assigned.

### Particle Reconciliation.

If particles overlap each other by too much, they can be forcibly removed from the simulation via

    Reconcile: Remove=[factor]

where *factor* is the overlap factor between particles: (R1 + R2 - r) / min(R1, R2) where R1, R2 are the radii of the particles, and r is the distance between them. This is currently the only reconciliation option.

## Conclusion

As development continues, more advanced configurations will become possible. This guide is current as of the time of its writing, on 4/4/2019. In the directory *GFlow/GFlowSimulation/configurations*, there are a number of sample configuration files.
