from   matplotlib import pyplot
from   textwrap   import wrap
import numpy
import matplotlib.ticker as ticker

#-------------------------------------READING_OF_FILE_AND_PREPARE_DATA_FOR_GRAPH--------------------------------------

whole_data    = []
num_of_elems  = []
lzw_ratio_no_dict  = []
lzw_ratio_dict     = []
zstd_ratio_no_dict  = []
zstd_ratio_dict     = []

with open('data/total_comp.txt') as file:
    for line in file:
        whole_data.append([float(x) for x in line.split()])

for x in (whole_data):
    num_of_elems.append(x[0])
    lzw_ratio_no_dict.append(x[1])
    lzw_ratio_dict.append(x[2])
    zstd_ratio_no_dict.append(x[3])
    zstd_ratio_dict.append(x[4])

#--------------------------------------------PIC_SETTINGS------------------------------------------------------------

fig, ax = pyplot.subplots (figsize = (16, 10), dpi = 100)

#-----------------------------------------SET_AXIS_LOCATORS-----------------------------------------------------------

ax.axis([0, max(num_of_elems) + max(num_of_elems)* 0.1, 0, 100])


#Set name of axis
ax.set_ylabel ("Ratio of compressed file to input, %", size = 20)
ax.set_xlabel ("Size of data, Kbts"      ,    size = 20)

#--------------------------------------SET_BUILDING_AREA_SET----------------------------------------------------------

ax.set_title ("\n".join(wrap('Comparing ratio of lzw and zstd on literature files', 60)), loc = 'center', size = 30)

ax.grid(which = 'major', color = 'gray')
ax.grid(which = 'minor', color = 'gray', linestyle = ':')
ax.minorticks_on ()


ax.plot(num_of_elems, lzw_ratio_no_dict, color = '#21421E', linewidth = 3, label = 'lzw ratio no dict(n)', marker = 'D', markevery = 0.005, mfc = '#21421E')
ax.plot(num_of_elems, lzw_ratio_dict, color = '#E32636', linewidth = 3, label = 'lzw ratio with dict(n)', marker = 'D', markevery = 0.05, mfc = '#E32636')
ax.plot(num_of_elems, zstd_ratio_no_dict, color = '#00FFFF', linewidth = 3, label = 'zstd ratio no dict(n)', marker = 'D', markevery = 0.005, mfc = '#00FFFF')
ax.plot(num_of_elems, zstd_ratio_dict, color = '#00FF00', linewidth = 3, label = 'zstd ratio with dict(n)', marker = 'D', markevery = 0.05, mfc = '#00FF00')

#Set legend
ax.legend(shadow = False, loc = 'upper right', fontsize = 20)

#--------------------------------------SAVING_OF_GRAPH------------------------------------------------------------

fig.savefig('pic/total_comp.png')
