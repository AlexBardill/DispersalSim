# PopEcoSim
A population ecology simulation I wrote to look at how the permeability ('stickyness') of a landscape affects individuals' investment in dispersal vs fecundity. 

The model has birth (and movement) and death events happening in series.  Individuals have only one ‘variable’ which is the investment in fecundity (F, limited to the range 0-1) with the remainder invested in dispersal (D = 1 – F).  Individuals inherit this strategy from their parent with a probability of mutation (1% here).

Birth events are density dependent with the probability of a birth occurring declining with density as well as depending on the investment strategy. There are two probability hurdles to pass – the density dependence with a probability of 1-(local density/K*2) and the fecundity investment level with a probability equal to the strategy. A new individual disperses at birth in a random direction with a distance drawn from a negative exponential distribution with a mean of the dispersal investment (i.e. 1- fecundity investment). This distance is then multiplied by a landscape permeability factor.
Death events occur with an even probability (here 0.1 unless otherwise stated) which is not density dependent.
Age is recorded, but has no influence on birth, movement or death events.
There are 2N events in each timestep.  
Simulations are initiated with 40000 individuals placed randomly on the landscape. Individuals start with a random value for F (0 - 1) and a random age (0 – 10).
