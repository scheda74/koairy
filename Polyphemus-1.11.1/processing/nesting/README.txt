How to use the example for nested simulations:
==============================================

First, the configuration files are meant to work with the chemistry module
named "Decay" which allows to have any number of species (gaseous or aerosol)
involved. Consequently, use the program "polair3d-decay". As for each example
file, launch the simulation from directory "processing/nesting".


Nested simulations are quite simple:

  * A first simulation is performed on a large domain. It is a normal
    simulation, except for what is saved. Indeed concentrations are only saved
    at the boundaries of the subdomain. The dimensions of the subdomain are
    given in the saver configuration file. All files regarding this simulation
    have "nesting" in their names, i.e. launch the simulation with
    "processing/nesting/polair3d-nesting.cfg".

  * A second simulation is performed on the subdomain. The concentrations
    saved during the first simulation are used as boundary conditions for the
    second simulation. All files regarding this simulation have "nested" in
    their names, i.e. launch the simulation with
    "processing/nesting/polair3d-nested.cfg".