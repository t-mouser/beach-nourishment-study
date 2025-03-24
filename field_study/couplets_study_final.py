""" 
This module is used to analyze the field data collected for the coastlines
of New Jersey (NJ). The data for this analysis is loaded from the file
field_study_couplets.csv and is available in a more the readable Microsoft
Excel file format including equations in the file field_study_couplets.xlsx.

To run this code, the data csv must be in the same directory as the code,
or you will need to change the csv read location in the "load the data" 
section.

Running this will create four plots exploring the difference in nourishment
efforts and efficiencies for couplets along the NJ coastline. If you wish to
save these plots to local folder, please use the section entitled "Save Files"
and set the variables save_files4 and save_files67 to True. You will also need
to create a directory structure to house your outputs as such:
    current_directory/figures/eps
    current_directory/figures/png

This code automates the save file name, but you can override this if desired 
in the same section as above.
"""

#%% Import Libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#%% Plot Styles / Sizes
# Sizes
m_size=200 #marker size
m_size_big=350
m_size75=75 #marker size 75
m_size100=100 #marker size 100
title_size = 28 # Title Size
axes_size = 24 #Axes Label Size
axes_size2=36
tick_size = 18
tick_size2=20
legendfont = 14
annote_size=14
fig_w_size=10 #width of figure
fig_h_size = 10 #height of figure

#colors
navy = 'navy'
crimson = 'crimson'
dark_gray='#2c2c2c'
blue='#332288'
d_green='#117733'
l_green = '#44AA99'
l_blue='#88CCEE'
slate_b='#446677'
yellow='#DDCC77'
pink='#CC6677'
rose='#AA4499'
marker1='o'
markerrest='o'
# Colormap for communities
colormap=np.array([rose,pink,l_green,blue])
shapes=np.array([marker1,markerrest,markerrest,markerrest])
categories=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,2,3,3,3])
mapping={0:'X',1:'v',2:'v',3:'v'}

#%%Load the data
df = pd.read_csv('field_study_couplets.csv')
# Grab the data
Municipality = df['RegionName']
couplets = df['Couplet']
House_wealth_diff = df['Housing_Wealth_Diff']
Sand_place_diff = df['Sand_Placed_m_Diff']
Beach_w_20_diff = df['Beach_Width_20_Diff']
Beach_w_00_diff = df['Beach_Width_00_Diff']
efficiency_diff = df['Efficiency_diff']
efficiency_diff_99_20=df['Eff_diff_99_20']
b_width_diff = df['beach_change_diff']
study_len = len(Municipality)
couplet_names = pd.DataFrame(list(couplets)).transpose().to_numpy()

#%% Variables
variable_names = ['Couplet', #0
                  'Difference in Mean House Value \n relative to lower wealth', #1
                  """Difference in Nourishment
Volume Per Meter ($m^3$/m)
Relative to Lower Wealth ($$-$)""", #2 \n relative to lower wealth
                  'Difference in Beach Width in 2000', #3 \n relative to lower wealth
                  'Difference in Beach Width in 2020', #4 \n relative to lower wealth
                  'Difference in Efficiency\nRelative to Lower Wealth ($$-$)', #5 \n relative to lower wealth
                  'Difference in Change of Beach Width'] #6 \n relative to lower wealth

#%% Change variables here
# 2000 to 2020 Plots
xvar4 = efficiency_diff
yvar4=Sand_place_diff
n_xvar4=5
n_yvar4=2
# 1999 to 2020 plots
xvar6=efficiency_diff_99_20
yvar6=Sand_place_diff
n_xvar6=5
n_yvar6=2

#%% Save Files
# The nonsequence of numerical identifiers is due to other plots
# being made in exploration of the data.
save_files4 = True # 2000 to 2020
save_files67 = False # 1999 to 2020

def file_name_formatter(orig):
    new=orig.replace(' ','_')
    new=new.replace('\n','')
    new=new.replace('(','')
    new=new.replace(')','')
    new=new.replace('$','')
    new=new.replace('-','')
    new=new.replace('^','')
    new=new.replace('/','')
    return new
    
#File names
# 2000 to 2020
fileName4_raw=variable_names[n_xvar4]+"_"+variable_names[n_yvar4]
fileName4 = file_name_formatter(fileName4_raw)

# 1999 to 2020
fileName67=variable_names[n_xvar4]+"_"+variable_names[n_yvar4]
fileName45 = file_name_formatter(fileName67)

    
#%% 2000 to 2020: Wealth and Nourishment Volume
plt.figure()
fig,ax1 = plt.subplots()
fig.set_figwidth(fig_w_size)
fig.set_figheight(fig_h_size)

# Towns / Communities Plots
# Turn this on for unique markers based on study region
for i in range(len(xvar4)):
    plt.scatter(xvar4[i],yvar4[i],m_size,marker='none')

# Set coordinates for zoomed version
# Do this here to set location box on main plot
zoom_coord_xmin=-0.05
zoom_coord_xmax=0.03
zoom_coord_ymin=-75
zoom_coord_ymax=75
# Calculate limits to make shading boxes
x_lims4 = ax1.get_xlim() # Gets x limits
y_lims4 = ax1.get_ylim() # Gets y limits
y_range = y_lims4[1]-y_lims4[0] # calculates the y range
ymin_calc = abs(y_lims4[0])/y_range # Gets the y min for shading box
boxmax_calc = (zoom_coord_ymax-y_lims4[0])/y_range # calculate shading as percentage of plot
boxmin_calc = (zoom_coord_ymin-y_lims4[0])/y_range # calculate shading as percentage of plot
ax1.set_xlim(x_lims4) # Set x limits to ensure that boxes fit well
plt.axvspan(0,x_lims4[1],ymin=ymin_calc,color='#dad9d9') # Plot shaded boxes
plt.axvspan(x_lims4[0],0,y_lims4[0],ymin_calc,color='#dad9d9') # Plot shaded boxes
plt.axvspan(zoom_coord_xmin,zoom_coord_xmax,boxmin_calc,
            boxmax_calc,facecolor="None",edgecolor='#ca0043') # Plot zoom box
# Plot horizontal and vertical vertex lines
plt.axhline(y=0, color=slate_b, linestyle = 'dotted')
plt.axvline(x=0, color=slate_b, linestyle = 'dotted')

# Plot data
# Plot markers
for i in range(len(xvar4)):
    plt.scatter(xvar4[i],yvar4[i],m_size,marker=mapping[categories[i]],color=colormap[categories[i]])
#Plot abbreviations
for i, txt in enumerate(couplets):
    plt.annotate(txt, (xvar4[i],yvar4[i]),color=dark_gray,fontsize=annote_size
                 ,ha='right',va='center',rotation=75)#horizontalalignment='center',verticalalignment='center')

# Set plot design
# Ticks
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
# labels
plt.xlabel(variable_names[n_xvar4], fontsize=axes_size)
plt.ylabel(variable_names[n_yvar4], fontsize=axes_size)
#Title
plt.title("Years: 2000 to 2020", fontsize=title_size)

# Save File if Desired
if save_files4 is True:
    plt.savefig('figures/eps/CPL4_'+fileName4+'.eps',format='eps')
    plt.savefig('figures/png/CPL4_'+fileName4+'.png',format='png', bbox_inches="tight")
plt.show()

#%% ZOOMED: 2000 to 2020: Wealth and Nourishment Volume
plt.figure()
fig,ax1 = plt.subplots()
fig.set_figwidth(fig_w_size)
fig.set_figheight(fig_h_size)

# Towns / Communities Plots
# Turn this on for unique markers based on study region
for i in range(len(xvar4)):
    plt.scatter(xvar4[i],yvar4[i],m_size,marker='none')

# Set zoomed limits based on box in previous plot
plt.xlim(zoom_coord_xmin,zoom_coord_xmax)
plt.ylim(zoom_coord_ymin, zoom_coord_ymax)

# Find limits to create new shading boxes
x_lims4 = ax1.get_xlim()
y_lims4 = ax1.get_ylim()
ymin_calc = abs(y_lims4[0])/(y_lims4[1]-y_lims4[0])
ax1.set_xlim(x_lims4)
# Shading boxes
plt.axvspan(0,x_lims4[1],ymin=ymin_calc,color='#dad9d9')
plt.axvspan(x_lims4[0],0,y_lims4[0],ymin_calc,color='#dad9d9')
# Vertical and Horizontal Axes lines
plt.axhline(y=0, color=slate_b, linestyle = 'dotted')
plt.axvline(x=0, color=slate_b, linestyle = 'dotted')

# Plot data as markers
for i in range(len(xvar4)):
    plt.scatter(xvar4[i],yvar4[i],750,marker=mapping[categories[i]],color=colormap[categories[i]])
# Annotate markers
for i, txt in enumerate(couplets):
    plt.annotate(txt, (xvar4[i],yvar4[i]),color=dark_gray,fontsize=18
                 ,ha='center',va='top',rotation=75)#horizontalalignment='center',verticalalignment='center')

# Plot Design
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
# Plot labels
plt.xlabel(variable_names[n_xvar4], fontsize=axes_size)
plt.ylabel(variable_names[n_yvar4], fontsize=axes_size)
#Title
plt.title("Years: 2000 to 2020", fontsize=title_size)

# Save Plot if Desired
if save_files4 is True:
    plt.savefig('figures/eps/CPL4_zoomed_'+fileName45+'.eps',format='eps')
    plt.savefig('figures/png/CPL4_zoomed_'+fileName45+'.png',format='png', bbox_inches="tight")
plt.show()
    

#%% 1999 - 2020: Difference Efficiency and Volume
plt.figure()
fig,ax1 = plt.subplots()
fig.set_figwidth(fig_w_size)
fig.set_figheight(fig_h_size)

# Towns / Communities Plots
# Turn this on for unique markers based on study region
for i in range(len(xvar6)):
    plt.scatter(xvar6[i],yvar6[i],m_size,marker='none')

# Set coordinates for zoomed version
# Do this here to set location box on main plot
zoom_coord_xmin=-0.06
zoom_coord_xmax=0.04
zoom_coord_ymin=-90
zoom_coord_ymax=90
# Calculate limits to make shading boxes
x_lims6 = ax1.get_xlim() # Gets x limits
y_lims6 = ax1.get_ylim() # Gets y limits
y_range = y_lims6[1]-y_lims6[0] # calculates the y range
ymin_calc = abs(y_lims6[0])/y_range # Gets the y min for shading box
boxmax_calc = (zoom_coord_ymax-y_lims6[0])/y_range # calculate shading as percentage of plot 
boxmin_calc = (zoom_coord_ymin-y_lims6[0])/y_range# calculate shading as percentage of plot
ax1.set_xlim(x_lims6) # Set x limits to ensure that boxes fit well
plt.axvspan(0,x_lims6[1],ymin=ymin_calc,color='#dad9d9') # Plot shaded boxes
plt.axvspan(x_lims6[0],0,y_lims6[0],ymin_calc,color='#dad9d9') # Plot shaded boxes
plt.axvspan(zoom_coord_xmin,zoom_coord_xmax,boxmin_calc,
            boxmax_calc,facecolor="None",edgecolor='#ca0043') # Plot zoom box
# Plot horizontal and vertical vertex lines
plt.axhline(y=0, color=slate_b, linestyle = 'dotted')
plt.axvline(x=0, color=slate_b, linestyle = 'dotted')

# Plot data
# Plot markers
for i in range(len(xvar6)):
    plt.scatter(xvar6[i],yvar6[i],m_size,marker=mapping[categories[i]],color=colormap[categories[i]])
#Plot abbreviations
for i, txt in enumerate(couplets):
    plt.annotate(txt, (xvar6[i],yvar6[i]),color=dark_gray,fontsize=annote_size
                 ,ha='right',va='center',rotation=75)#horizontalalignment='center',verticalalignment='center')

# Set plot design
# Ticks
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
# labels
plt.xlabel(variable_names[n_xvar6], fontsize=axes_size)
plt.ylabel(variable_names[n_yvar6], fontsize=axes_size)
#Title
plt.title("Years: 1999 to 2020", fontsize=title_size)

# Save File if Desired
if save_files67 is True:
    plt.savefig('figures/eps/Y1999_CPL4_'+fileName4+'.eps',format='eps')
    plt.savefig('figures/png/Y1999_CPL4_'+fileName4+'.png',format='png', bbox_inches="tight")
plt.show()

#%% 1999 - 2020: Difference Efficiency and Volume ZOOMED IN
plt.figure()
fig,ax1 = plt.subplots()
fig.set_figwidth(fig_w_size)
fig.set_figheight(fig_h_size)

# Towns / Communities Plots
# Turn this on for unique markers based on study region
for i in range(len(xvar6)):
    plt.scatter(xvar6[i],yvar6[i],m_size,marker='none')

# Set zoomed limits based on box in previous plot
plt.xlim(zoom_coord_xmin,zoom_coord_xmax)
plt.ylim(zoom_coord_ymin, zoom_coord_ymax)

# Find limits to create new shading boxes
x_lims6 = ax1.get_xlim()
y_lims6 = ax1.get_ylim()
ymin_calc = abs(y_lims6[0])/(y_lims6[1]-y_lims6[0])
ax1.set_xlim(x_lims6)
# Shading boxes
plt.axvspan(0,x_lims6[1],ymin=ymin_calc,color='#dad9d9')
plt.axvspan(x_lims6[0],0,y_lims6[0],ymin_calc,color='#dad9d9')
# Vertical and Horizontal Axes lines
plt.axhline(y=0, color=slate_b, linestyle = 'dotted')
plt.axvline(x=0, color=slate_b, linestyle = 'dotted')

# Plot data as markers
for i in range(len(xvar6)):
    plt.scatter(xvar6[i],yvar6[i],750,marker=mapping[categories[i]],color=colormap[categories[i]])
# Annotate markers
for i, txt in enumerate(couplets):
    plt.annotate(txt, (xvar6[i],yvar6[i]),color=dark_gray,fontsize=18
                 ,ha='center',va='top',rotation=75)#horizontalalignment='center',verticalalignment='center')

# Plot Design
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
# Plot labels
plt.xlabel(variable_names[n_xvar6], fontsize=axes_size)
plt.ylabel(variable_names[n_yvar6], fontsize=axes_size)
#Title
plt.title("Years: 1999 to 2020", fontsize=title_size)

# Save Plot if Desired
if save_files67 is True:
    plt.savefig('figures/eps/Y1999_CPL4_zoomed_'+fileName45+'.eps',format='eps')
    plt.savefig('figures/png/Y1999_CPL4_zoomed_'+fileName45+'.png',format='png', bbox_inches="tight")
plt.show()
