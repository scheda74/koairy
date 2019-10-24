
-- The number of species.
Number_species = 1

-- Particle density in g/cm3
Particle_density = 1.4

default = {
   Number_monte_carlo = 100000,

   diameter = { min = 0.01, max = 10.0, Nb = 5 , Nf = 0},

   --section_definition_auto = false,

   composition_section_0 = { S0 = {0, 1.0} },

   -- The number of bins is given either by number of size bins
   -- times composition bins or by the number of GS_* sections.

   general_section_0 = { 0, 0},
}

--im = {
  --Number_monte_carlo = 10000,

   --diameter = { min = 0.001, max = 10.0, Nb = 10},

   --composition_section_0 = { }
--}
