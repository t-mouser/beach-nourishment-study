# Model Documentation
This document outlines how to use the model to create regime diagrams, single run analysis, and future predictive analyses. 

All data files are stored in the mat-data directory. This is where the code will save .mat data files and load data files from as well.

## Modeling Scripts (use to generate model data)
### Parametric Analyses
Use parametric_analyses.m to create regime diagrams based on different parameter inputs. 

**How to Run:** Input the parameter values in the top half of the script and run the script.

### Single Run Analysis
Use single_run.m to plot beach widths and benefits over time for specific community pairings. These pairings must be selected from a .mat file created by parametric_analyses.m, or the function inputs must be supplied by the user creating a standalone data file for this purpose (using the same input variables created in parametric_analyses.m). Please note that using individually created data files will not contain the optimization data output from parametric_analyses.m.

**How to Run:** Select the .mat data file you want to use to create the subplots. Load the property value vectors (C1_vec and C2_vec) to pick which property values you would like to use. Use their index number and hit run. This will output beach widths and net benefits over time.

### Future Analyses
Use future_analysis.m to create regime diagrams that look at specific community pairings with rising sand costs and erosion rates. This may be tweaked to explore increases in other coastal factors/stressors.  

**How to Run:** Input the parameter values in the top half of the script and run the script.

## Visualization Scripts (use to generate figures to create plots of data)
### Parametric Figures (no efficiency metric) 
Use parametric_figures.m to create regime diagram pseudocolor plots of behaviors, beach widths, and nourishment volumes. Other data outputs not pertaining to efficiency can be plotted here as well. There is no efficiency metric calculated in maincode due to our requirement of determining the Reference Value (W<sub>R</sub>) based on multiple model runs. 

Note: It is possible to make future figures in this file as well, though it would have to be slightly rewritten.

### Future Figures (no efficiency metric) 

This code creates same figures as parametric_analyses.m, specifically for future analyses in data-file future_analysis.m.

### Efficiency
Use efficiency.m to calculate nourishment efficiency and plot regime diagrams relating to nourishment efficiency. 

Note: This script can produce figures for future efficiency as well, but x and y labels will have to changed, as will vertex scaling.

### Future Efficiency
Use future_efficiency.m to do the same as in efficiency.m, except for future projections. It will calculate nourishment efficiency and plot regime diagrams relating to nourishment efficiency. 

## File Explanations
### Maincode
This code is a script containing a function called by other matlab files (parametric_analyses.m, future_analysis.m). Maincode runs each model run based on the input parameters from the other file. This function contains the equations calculating the shoreline and shore toe locations over time, benefits, behaviors, nourishment volumes, and nourishment widths. This code also performs the optimizations for the three different coordination levels: coordination, conservative non-coordination, and risky non-coordination. It returns matrices of values to build regime diagrams. 

### Parametric Analyses
This code creates data for regime diagrams based on a range of different baseline property values. It runs the maincode function based on the parameter inputs contained in this script. Matrices of the outputs from maincode are saved based on the three optimization types: coordination, conservative non-coordination, and risky non-coordination.

### Future Analysis
This code creates data for regime diagrams based on increasing erosion and sand cost rates. It also runs the maincode based on parameter inputs and outputs matrices with data pertaining to beach widths, net benefits, and rotation intervals. It selects which model runs to save based on three optimization types: coordination, conservative non-coordination, and risky non-coordination.

### Single Runs
This code runs single runs based on .mat data files. You select the community Property Value per community you would like to plot, and the maincode function (a duplicate contained within this file) runs and plots beach widths and net benefits over time. Parametric Analyses need to be run prior to utilizing this tool, or existing .mat files must be utilized.

### Figures Code (parametric, future)
This code is used to plot regime diagrams relating to nourishment volume and community behaviors. Future analysis data can be plotted in parametric_figures.m, but the x and y axes need to be adjusted to the new vectors (sand cost and erosion rate). For this reason, future_figures was also added as it is ready to use.

### Efficiency (parametric, future)
The efficiency codes use data files produced by analysis (parametric or future) and calculates efficiency and plots efficiency regime diagrams. While parametric_efficiency.m can analyze future_analysis.m data files, the x and y labels and values have to be changed (similar to figures code). As a result, future_efficiency.m was added so that it can be used without altering the code.

## Directories
### mat-data
The mat-data directory contains output data produced by the model. If you run code with this model, it will save the code in the mat-data directory, so it is necessary that you set up a directory structure in the same way.

### colormaps
The colormaps folder contains arrays that are used to create unique colormaps to shade the pseudocolor plots. These colormaps are loaded anywhere figures are created (figures and efficiency codes), and therefore need to be included in your directory structure for this code to operate correctly. 

Additionally, there is a create_colormaps.m file that is included in this folder which will allow you to create new colormaps should you wish to do so.

