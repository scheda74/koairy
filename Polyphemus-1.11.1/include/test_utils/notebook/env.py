from glob import glob
import os
from os import path

import matplotlib as mpl
import matplotlib.pyplot as plt

import atmopy
from atmopy.display import *

from test_utils import *
from test_utils.config_helpers import *

from IPython import get_ipython
ipython = get_ipython()
ipython.magic("matplotlib inline")

POLYPHEMUS = path.realpath(os.environ["POLYPHEMUS"])
assert path.isdir(POLYPHEMUS + "/processing")
current_path = os.getcwd()

print "Polyphemus is at", POLYPHEMUS
print "Atmopy is at", path.dirname(atmopy.__file__)
print "This notebook is at", os.getcwd()

# This symbolic link is used in the notebook html
# to create urls to Polyphemus:
os.system("ln -sf " + POLYPHEMUS + " polyphemus")
