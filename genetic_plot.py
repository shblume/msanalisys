''' Numpy general plotter. '''

# Libraries imported 

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects
from plot_ease import *
from collections import OrderedDict as odict

# User defined variables

gen = 6000;
output = 'image_result.png'

# Changing working directory

os.chdir('/home/earaujo/Repositories/HLA/workbench/gen_a/best.mi')
# Functions defined

# Program execution
gen_dict = odict()
for index in range(0, gen):
    name = 'best.{}.npy'.format(str(index))
    best = np.load(name)
    gen_dict[index] = best

keys = list(gen_dict.keys()); values = list(gen_dict.values());

fig, ax = splot(values, None, lw=1300,
                xname='Generations', yname='MI', tname='Transinformation per generation')


plt.savefig(output)
