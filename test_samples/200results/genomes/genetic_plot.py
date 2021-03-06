''' Numpy general plotter. '''

# Libraries imported 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects
from plot_ease import *

# User defined variables

gen = 200; gen = int(gen)
#output = entry+'.jpg'

# Program execution

plot_list = []; name_list = []
for index in range(0, gen):
     a_arr_name = 'genome.'+str(index)+'.npy'; name_list.append(str(index))
     b_arr_name = 'genome.'+str(index+1)+'.npy'
     a_array = np.load(a_arr_name); b_array = np.load(b_arr_name)
     plot_list.append(a_array-b_array)

multiplot(name_list, plot_list, xname='generation', yname='fitness',note='bh',a=0.7,y_min=0)


plt.show()
