# ARICHSim

Code to identify particles from ARICH data by using a fast optical photon Monte-Carlo simulation and a likelihood approach.

This code runs using ROOT 6.14/06. Some of the scripts in analysis.cpp rely on the outputs of DetectorSim, but these scripts are not essential.

The general principle:
We take in as input a distribution of photons detected in the EMPHATIC ARICH detector, and measurements of a particle's momentum, position, and direction from the particle trackers. Given these values, we want to determine what particle this is. This can be accomplished by running a Monte Carlo simulation of each particle hypothesis, and comparing the resulting photon distribution to the real experimental distrubtion. 


## Installation:
To build, simply check out this repository, source ROOT, then run *make* in the *ARICHSim* directory. 

This will create the executable *./bin/arichsim*


## Running:
ARICHSim is controlled using command line arguments. Running *./bin/arichsim* with no arguments will print out a list of possible modes to run in, and what their parameters are. These are listed below in more detail:

|Mode description| Arguments |
| --- | ---|
|Identify particle from photon histogram and particle information:       | `-id <eventfile> <histname> <mom> <xdir> <ydir> <xpos> <ypos>`  |
|Save expected photon distribution, given velocity `beta` of particle | `-b <outputdir> <beta> [xdir ydir xpos ypos]` |
|Save expected photon distribution, given particle type and momentum | `-p <outputdir> <pid> <mom> [xdir ydir xpos ypos]`  |
|Test particle identification by simulating an event | `-pid <pid> <mom> [xdir ydir xpos ypos]`  |
|Test multi-particle identification by simulating ndetected random events | `-mp <ndetected>`  |
|Given output from the Geant4-based DetectorSim, project photons and produce histogram of photon hits| `-gpdf <g4filename> <pid>`  |
|Run through a set of simulated event photon distributions, and evaluate the distribution of loglikelihood ratios for each particle hypothesis| `-s <analysisdir> <mom> [xdir ydir xpos ypos]`  |

|Parameter | Explanation |
| --- | ---|
|`eventfile` | ROOT file containing photon distribution information |
|`histname` | Name of ROOT histogram (`TH2D`) in eventfile containing a single event's photon data. The histogram should be 48 by 48 bins, representing the the array of PMT pixels in the ARICH detector. Each entry should represent a single photon hit on a pixel|
|`mom` | Momentum of particle, in GeV/c|
|`xdir` | X-component of unit direction vector of particle|
|`ydir` | Y-component of unit direction vector of particle|
|`xpos` | X-Position (in cm) of where particle enters aerogel, relative to its center |
|`ypos` | Y-Position (in cm) of where particle enters aerogel, relative to its center |
|`outputdir` | Output directory to save files to|
|`beta` | Velocity of particle (relative to *c*)|
|`pid` | Particle ID: 0 is a pion+, 1 is a kaon+, 2 is a proton|
|`ndetected` | Number of independent particles to simulate in multiparticle test|
|`g4filename` | Name of ROOT TTree file containing info on optical photons, produced with DetectorSim|
|`analysisdir` | Directory containing three DetectorSim-produced TTrees containing optical photon info: Pion.root, Kaon.root, Proton.root|



## Code Structure
### Analysis
The main file for particle identification.
Contains functions for running particle identification, as well as testing out particle identification with data generated from either ARICHSim itself, or DetectorSim.

### Arich
Class representing a simulation of the ARICH detector. Contains two Aerogel objects and a Detector object as members. Contains methods to simulate photon distributions, given a set of parameters including the velocity, position, and direction of a particle.

### Beam
Class that generates particles with randomly thrown positions and directions

### Aerogel
Class that represents the aerogel layers included in the ARICH detector, as well as the optical processes that occur inside. This is where a particle will generate set of photons, and where those photons will then scatter and refract around.

### Detector
Class representing the PMT array. Contains information about quantum efficiency, and method to project photons onto a histogram representing the array.

### Particle
A particle, with a position and direction.

### Photon
A photon, with a position, direction, and wavelength