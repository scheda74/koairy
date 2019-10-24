import os, sys
sys.path.insert(0, '../')
from run.ensemble_generation import config, EnsembleProgram
sys.path.pop(0)

####################################################
# An example to launch several Polyphemus programs #
####################################################

# Polyphemus directory.
polyphemus_dir = '../../..'
polyphemus_dir = os.path.abspath(polyphemus_dir)

# Reads the configuration files.
epr = EnsembleProgram('parameter.cfg', 'program.cfg',
                      only_preprocessing = True)

# Sets the Polyphemus directory.
epr.SetPolyphemusDirectory(polyphemus_dir)

# Prints the model number.
print "Nmodel: ", epr.parameter.Nmodel

# Generates a new ensemble.
epr.parameter.GenerateIdEnsemble()

# Sets the ConfigReplacement object.
config_replacement = config.ConfigPolair3D()
epr.SetConfigReplacement(config_replacement)

# Gets the dictionary with all variables wich will be replaced in the generic
# configuration files for the first model.
d = epr.GetGeneralDict(0)

# Gets an instance 'Polyphemus' and creates all directories where you will
# have the results.
ens = epr.GetEnsemble(group = 'polyphemus')

# The list of programs names.
program_name_list = []
for program in ens.program_list:
    program_name_list.append(program.name.split('/')[-1])

print "\nNprogram: ", len(ens.program_list)

# Gets the available hosts from your hosts list (in the file
# '$HOME/.dsh/group/polyphemus').
load_averages = ens.net.GetAvailableHosts()

# You can launch all programs.
print "\nYou can launch all programs with the method:"
print "ens.RunNetwork()"
#ens.RunNetwork()
