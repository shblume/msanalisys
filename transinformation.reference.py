""" This program analyses a Multiple Sequence Alignment file, to give us informations about the conservation and correlation of the present sequences using Shannon's Information Theory. This program results in 5 graphic plots and a .xlsx file for analisys. """

# Libraries imported:
try:
    import numpy as np
    import statistics as st
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import matplotlib.patheffects
    from Bio import AlignIO
    from collections import OrderedDict as odict
except ImportError:
    print(''' One or more of the libraries needed to the execution of this program are not yet installed or malfunctioning.
     Check if the following Python extensions are installed and perfectly working:
     Numpy, Matplotlib.Pyplot, BioPython, Statistics, Pandas, Collections
     If the functioning of this script is not clear and for tutorials on how to install those libraries, go to;
     "https://github.com/shblume". ''')
    raise

# Constant defined:

ALPHABET = odict()
ALPHABET = {"A": 0, "R": 1, "N": 2, "D": 3, "Q": 4,
          "E": 5, "G": 6, "H": 7, "L": 8, "K": 9,
          "M":10, "F":11, "S":12, "T":13, "W":14,
          "Y":15, "C":16, "I":17, "P":18, "V":19,
          "-":20, ".":20, "B": 2, "Z": 4, "X":20, "J":20}
Q = 21
VLIST = list(ALPHABET.values())

# Functions defined:

def reader(fname):
    """ Transforms a given clustalx file onto a matrix based on a user
        defined alphabet. """
    psdmtx = []
    for record in AlignIO.read(fname+'.aln','clustal'):
        sequence = record.seq
        line = list(sequence) 
        psdline = []
        for item in line:
            psdline.append(ALPHABET[item])
        psdmtx.append(psdline)
    matrix = np.array(psdmtx).reshape(len(psdmtx),len(psdmtx[0]))
    return matrix

def accountant(col,mtx):
    """ Counts the frequence proportion in a column for a given matrix. """
    count = np.zeros((Q)); #print(count)
    for index in range(0, len(mtx)):
        line = mtx[index]
        aa = line[col]; aa = int(aa) 
        count[aa] += 1
    return count/len(mtx)

def entropy(prob):
    """ Calculates the shannon entropy for a given value. """
    entropy = 0
    for item in prob:
        if item != 0:
            entropy += -1*(item*np.log2(item))
        elif item == 0:
            entropy += 0
    return entropy

def conjprob(col_i,col_j,mtx):
    """ Counts the joint proportion in two columns for a given matrix. """
    count = np.zeros((Q,Q));
    for ind in range(0,len(mtx)):
        line = mtx[ind];  
        aa_i = line[col_i]; aa_i = int(aa_i)
        aa_j = line[col_j]; aa_j = int(aa_j)
        for ind1 in range(0,Q):
            for ind2 in range(0,Q):
                if (aa_i, aa_j) == (ind1, ind2):
                    count[ind1][ind2] += 1
#    print('count:\n{}'.format(count))

    return count/(len(mtx))
    
def mutualinfo(prob_i,prob_j,prob_ij):
    """ Calculates the transinformation between two proportions and the
        joint proportion. """
    MI = 0
    for index1 in range(0, len(prob_i)):
        for index2 in range(0, len(prob_j)):
            if prob_i[index1]*prob_j[index2] != 0:
                if prob_ij[index1][index2] != 0:
                    MI += prob_ij[index1][index2]*np.log2(prob_ij[index1][index2]/(prob_i[index1]*prob_j[index2]))
                elif prob_ij[index1][index2] == 0:
                    MI += 0
    return MI

def normalizer(info,entropy):
    """ Normalizes a value for another, considering the divisions by zero. """
    if entropy == 0:
        if info == 0:
            normal = 0
        elif info > 0: 
            normal = 'positive infinity'
        elif info < 0:
            normal = 'negative infinity'
    elif entropy != 0:
        normal = (info/entropy)
    return normal

def multiplot(x_axis,y_axis, xname=None, yname=None, tname=None, note='ko', y_max=None, x_max=None,subplot=None,coord=None, y_scale='linear'):
    """ Pyplot plottage for any instance. """
    # For simple plottage
    if subplot == None:
        fig, ax = plt.subplots()
        ax.plot(x_axis,y_axis,note)
        ax.set(xlabel=xname, ylabel=yname, title=tname)
        ax.set_yscale(y_scale)
        if y_max != None:
            plt.ylim(ymax=y_max)
        if x_max != None:
            plt.xlim(xmax=x_max)
    
    # For multiplotage
    else:
        print(' Warning: for subplots you must call "fig, ax = plt.subplots()" at least to run the code.')
        if coord == None:
            axs[subplot].plot(x_axis,y_axis,note)
            axs[subplot].set(xlabel=xname, ylabel=yname, title=tname)
            axs[subplot].set_yscale(y_scale)
            if y_max != None:
                plt.ylim(ymax=y_max)
            if x_max != None:
                plt.xlim(xmax=x_max)
        else:
            axs[subplot][coord].plot(x_axis,y_axis,note)
            axs[subplot][coord].set(xlabel=xname, ylabel=yname, title=tname)
            axs[subplot][coord].set_yscale(y_scale)
            if y_max != None:
                plt.ylim(ymax=y_max)
            if x_max != None:
                plt.xlim(xmax=x_max)

def center_spines(ax=None, centerx=0, centery=0,grdval=True):
    """Centers the axis spines at <centerx, centery> on the axis "ax", and
    places arrows at the end of the axis spines."""
    if ax is None:
        ax = plt.gca()

    # Set the axis's spines to be centered at the given point
    # (Setting all 4 spines so that the tick marks go in both directions)
    ax.spines['left'].set_position(('data', centerx))
    ax.spines['bottom'].set_position(('data', centery))
    ax.spines['right'].set_position(('data', centerx - 1))
    ax.spines['top'].set_position(('data', centery - 1))

    # On both the x and y axes...
    for axis, center in zip([ax.xaxis, ax.yaxis], [centerx, centery]):
        # Turn on minor and major gridlines and ticks
        axis.set_ticks_position('both')
        axis.grid(grdval, 'major', ls='solid', lw=0.5, color='gray')
        axis.grid(grdval, 'minor', ls='solid', lw=0.1, color='gray')
        axis.set_minor_locator(mpl.ticker.AutoMinorLocator())

def point_show(result_a,result_b,a_cut=None,b_cut=None):
    higher = [];
    k_list = list(result_a.keys())
    if a_cut == None and b_cut == None:
        a_val = list(result_a.values()); b_val = list(result_b.values()) 
        a_cut = st.median(a_val); b_cut = st.median(b_val)
        for key in k_list:
            if result_a[key] >= a_cut:
                if result_b[key] >= b_cut:
                    higher.append(key) 
    else:
        a_cut = float(a_cut); b_cut = float(b_cut)
        for key in k_list:
            if result_a[key] >= a_cut: 
                if result_b[key] >= b_cut: 
                    higher.append(key)
    return higher


# Execution:

print(' ====================================================================================================') 
print('\n     Multiple Sequence Alignment analyzer')
print('     Developed by: "shblume"')
print('     https://github/shblume')
print('\n ====================================================================================================')

try:
    # Opens the alignment and transforms into a matrix.
    aln_file = input('\n Please, insert the alignment in clustal file format (".aln"): ')
    matrix = reader(aln_file); 
except:
    print(''' The file you just entered is not valid or is not in the same path as this program.
      The correct file must be a clustal file, with extension ".aln". You must input it's name ONLY,
      as a example for a file named "alignment.aln" just input "alignment".
      If the use of this program is not clear to you, visit the link "https://github/shblume" to
      know how to use the program. ''')

# Sets the reference channel as the last column and calculates all useful values
column_j = -1
count_j = accountant(column_j,matrix);
entropy_j = entropy(count_j)
print(' The entropy of the reference column is: {} bit.'.format(entropy_j))
entro = odict (); transin = odict(); normali = odict()

# Calculates all values of entropy and transinformation for all other columns
for column_i in range(0, len(matrix[0])):
    count_i = accountant(column_i,matrix)
    entropy_i = entropy(count_i)
    conj_prob = conjprob(column_i,column_j,matrix); 
    MI = mutualinfo(count_i,count_j,conj_prob); 
    normalized = normalizer(MI,entropy_i) 
    print(' Column {} out of {}: H = {}bit, I = {}bit.'.format(column_i,(len(matrix[0])-1),entropy_i,MI))

    # Appends the results to lists who will be used to generate plots
    entro[column_i+1] = entropy_i
    transin[column_i+1] = MI 
    normali[column_i+1] = normalized

# In the plottages:
plot_inp = input('\n Plot the normalized transinformation by entropy centered by median of all the data? [y/n]: ')
plot_names = list(entro.keys())
normal_plot = list(normali.values())
trans_plot = list(transin.values())
entropy_plot = list(entro.values())

# Calculates the average in the entropy and normalized transinformation result 

# biased on the median
if plot_inp == 'y' or plot_inp == 'Y':
    median_norma = st.median(normal_plot); median_entro = st.median(entropy_plot);
    plt.figure(1)
    multiplot(entropy_plot,normal_plot,xname='Entropy',yname='Normalized transinformation',tname='Graph for normalized transinformation and entropy',note='k+')
    center_spines(centerx = median_entro, centery = median_norma)
    plt.savefig(aln_file+'.entropy_trans.png')
    accept = point_show(entro,normali)
    print('\n The points accepted between the cutoffs are: \n "{}"'.format(accept))    

# Or calculates it biased on some user's input
else:
    center_entro = input(' Enter the center for entropy: '); center_entro = float(center_entro)
    center_norma = input(' Enter the center for normalized transinformation: '); center_norma = float(center_norma)
    plt.figure(1)
    multiplot(entropy_plot,normal_plot,xname='Entropy',yname='Normalized transinformation',tname='Graph for normalized transinformation and entropy',note='k+')
    center_spines(centerx = center_entro, centery = center_norma)
    plt.savefig(aln_file+'.entropy_trans.png')
    accept = point_show(entro,normali,a_cut=center_entro,b_cut=center_norma)
    print('\n The points accepted between the cutoffs are: \n "{}"'.format(accept))


# Generate plots for each result and a plot for mean normalized and entropy
plt.figure(2)
fig, axs = plt.subplots(2, 2, tight_layout=False)

multiplot(plot_names,entropy_plot,tname='entropy',xname='Columns',yname='Entropy',note='go',subplot=0,coord=0)
multiplot(plot_names,trans_plot,tname='transinformation',xname='Columns',yname='Transinformation',note='yo',subplot=1,coord=0)
multiplot(plot_names,normal_plot,tname='normalized transinformation',xname='Columns',yname='Normalized transisnformation',note='bo',subplot=0,coord=1)
multiplot(entropy_plot,normal_plot,tname='normalized vs entropy',xname='Entropy',yname='Normalized Transinformation',note='ko',subplot=1,coord=1)

plt.savefig(aln_file+'.result.png')

# Generates the result data frame and exports it as a ".xlsx" file
df_plot = []; ac_plot = []
col_names = ['entropy (H)','transinformation (I)', 'I/H']
for item in plot_names:
    df_plot.append([entro[item],transin[item],normali[item]])

for item in accept:
    ac_plot.append([entro[item],transin[item],normali[item]])

writer = pd.ExcelWriter(aln_file+'.result.xlsx')
df1 = pd.DataFrame(df_plot,index=plot_names, columns=col_names)
df2 = pd.DataFrame(ac_plot,index=accept, columns=col_names)
df1.to_excel(writer, 'Results')
df2.to_excel(writer, 'Above cutoff values')
writer.save()
writer.close()

print('\n Program succesful! Check your file.')
print(' Program closing...')

""" End of the code. """
