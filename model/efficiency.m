%{
This code is used to calculate nourishment efficiencies from output 
matrices created by parametric_analyses.m. 

You will need to select and load the .mat data file you wish to analyze in 
the two sections below. From there, the script will calculate the 
nourishment efficiency and plot them in a regime diagram.

The reference value utilized in the nourishment efficiency equation is
calculated by analyzing inputs from all of the data files listed in the
variable "efficiency_files". 
%}

%% Data Files
dir_loc='mat-data/'; % the subdirectory where .mat files are located
% Place your mat files in the efficiency_files variable as a vector
% Vector of data files to be studied
efficiency_files=["parametric_analysis_cs_bs_v1.mat",... % #1
    "parametric_analysis_cs_ba_v1.mat",... % #2
    "parametric_analysis_cs_ba_v2.mat",... % #3
    "parametric_analysis_ca_bs_v2.mat",... % #4
    "parametric_analysis_ca_ba_v1.mat",... % #5
    "sensitivity_analysis_cs_ba_v1_5yr_avg"]; % #6

%% IMPORTANT Pick your file name by changing the variable below
file_number=2;

%% Determine Reference Value from Data Sources
nFiles = length(efficiency_files); % Number of data sources being evaluated for entire study
filenames=string.empty(0,nFiles); % Will hold strings of your files names
gama_data=NaN(nFiles); % Container to house gamma values from each data file
tmax_data=NaN(nFiles); % Container to house tmax data from each data file
w_init_data=NaN(nFiles); % Container to house w_init_data from each data file

% This loop goes through each data file and pulls the necessary data to
% calculate the reference width
for iData=1:nFiles     
    filenames(iData) = string(strcat(dir_loc,efficiency_files(iData))); % Get the file name and make sure its a string
    load(filenames(iData),"gamma_array","tmax","w_init"); % Load the datafile
    gama_data(iData)=max(gamma_array); % Gets the maximum gamma value for that dataset (is an array in the dataset)
    tmax_data(iData)=tmax; % Gets the tmax value for that dataset (one value in dataset)
    w_init_data(iData)=min(w_init); % Gets the min w_init for that dataset (is an array in the dataset)
end
ref_gamma=max(max(gama_data)); % Finds the max gamma value from all datasets
ref_tmax=max(max(tmax_data)); % Finds the max time horizon for all datasets
ref_w_init=min(min(w_init_data)); % Finds the minimum initial width for all datasets

%% Load the dataset 
load(filenames(file_number)); % Loads your file from the array of filenames used above to calculate reference for efficiency
ref_eff = abs(ref_w_init - ref_gamma * ref_tmax); % Caculates the reference point for the efficiency metric

%% Load Colormaps
% Colors used for behaviors
colorfile = matfile("colormaps/behaviors.mat"); 
colorize = colorfile.custommap;
colorizerev = colorfile.custommaprev;
% Colors used for subtractive analysis - midpoint is white
opeqcolorfile = matfile("colormaps/contrast_equal.mat"); 
opeq_color = opeqcolorfile.op_equal_cmap;
opeq_colorrev = opeqcolorfile.op_equal_cmap_rev;
% Colors representing Community 1
com1_colorfile=matfile("colormaps/community_1.mat"); 
com1_color= com1_colorfile.com1_cmap;
com1_colorrev= com1_colorfile.com1_cmap_rev;
% Colors representing Community 2
com2_colorfile=matfile("colormaps/community_2.mat"); 
com2_color= com2_colorfile.com2_cmap;
com2_colorrev= com2_colorfile.com2_cmap_rev;

%% General Plot Variables
% font sizes
axis_font=16; % axes font size
title_font=16; % Title Font Size
gca_font=14; % General plot font size
sg_font=18; % SG title font size
legendboxes=16; % Legend font size
MarkSize=40; % Marker size

% x and y axes labels
xlabel_c2vec='Community 2 Property Values ($1M)';
ylabel_c1vec='Community 1 Property Values ($1M)';

% Make transect lines
% Based on full property value
xtransect=[C2_vec(1) C2_vec(end)];
ytransect=[C1_vec(1) C1_vec(end)];
% Based on scientific notation PV
xtransect_e6=[C2_vec(1)/1e6 C2_vec(end)/1e6];
ytransect_e6=[C1_vec(1)/1e6 C1_vec(end)/1e6];

%% Calculation Storage Containers
% Changes in beach widths from reference
dw1_ref_cord=NaN(nn,mm);
dw2_ref_cord=NaN(nn,mm);
dw1_ref_cons=NaN(nn,mm);
dw2_ref_cons=NaN(nn,mm);
dw1_ref_risk=NaN(nn,mm);
dw2_ref_risk=NaN(nn,mm);
% Change in beach widths from w_init
dw1_init_cord=NaN(nn,mm);
dw2_init_cord=NaN(nn,mm);
dw1_init_cons=NaN(nn,mm);
dw2_init_cons=NaN(nn,mm);
dw1_init_risk=NaN(nn,mm);
dw2_init_risk=NaN(nn,mm);
% Efficiencies
E1_cord=NaN(nn,mm);
E2_cord=NaN(nn,mm);
E1_cons=NaN(nn,mm);
E2_cons=NaN(nn,mm);
E1_risk=NaN(nn,mm);
E2_risk=NaN(nn,mm);
% Differences of efficiencies relative to other community 
E1_E2_cord_diff=NaN(nn,mm);
E1_E2_cons_diff=NaN(nn,mm);
E1_E2_risk_diff=NaN(nn,mm);
% Differences of efficiencies relative to Non-Coordination
E1_cord_cons_diff=NaN(nn,mm);
E2_cord_cons_diff=NaN(nn,mm);
E1_cord_risk_diff=NaN(nn,mm);
E2_cord_risk_diff=NaN(nn,mm);
% Differences in Nourishment Volumes relative to other community 
V1_V2_cord_diff=NaN(nn,mm);
V1_V2_cons_diff=NaN(nn,mm);
V1_V2_risk_diff=NaN(nn,mm);
% Differences in Nourishment Volumes relative to Non-Coordination
V1_cord_cons_diff=NaN(nn,mm);
V2_cord_cons_diff=NaN(nn,mm);
V1_cord_risk_diff=NaN(nn,mm);
V2_cord_risk_diff=NaN(nn,mm);

%% Binary Storage Containers
% Efficiency Relative to other community 
E1_E2_binary_cord=NaN(nn,mm);
E1_E2_binary_cons=NaN(nn,mm);
E1_E2_binary_risk=NaN(nn,mm);
% Efficiency Relative to Non-Coordination
E1_binary_cord_cons=NaN(nn,mm);
E2_binary_cord_cons=NaN(nn,mm);
E1_binary_cord_risk=NaN(nn,mm);
E2_binary_cord_risk=NaN(nn,mm);
% Volume Relative to other community 
V1_V2_binary_cord=NaN(nn,mm);
V1_V2_binary_cons=NaN(nn,mm);
V1_V2_binary_risk=NaN(nn,mm);
% Volume Relative to Non-Coordination
V1_binary_cord_cons=NaN(nn,mm); 
V2_binary_cord_cons=NaN(nn,mm);
V1_binary_cord_risk=NaN(nn,mm); 
V2_binary_cord_risk=NaN(nn,mm);
% Vol & Eff Analysis relative to other community
C1_C2_EV_cord=NaN(nn,mm);
C1_C2_EV_cons=NaN(nn,mm);
C1_C2_EV_risk=NaN(nn,mm);
% Vol & Eff Analysis relative to Non-Coordination
E1_V1_cord_cons=NaN(nn,mm);
E2_V2_cord_cons=NaN(nn,mm);
E1_V1_cord_risk=NaN(nn,mm);
E2_V2_cord_risk=NaN(nn,mm);

%% Loops to calculate efficiency
for ii=1:nn 
    for jj=1:mm       
        %% Calculate beach gains or losses in relation to initial width
        dw1_init_cord(ii,jj)=avg_w1_cord_pa(ii,jj)-w_init1;
        dw2_init_cord(ii,jj)=avg_w2_cord_pa(ii,jj)-w_init2;
        dw1_init_cons(ii,jj)=avg_w1_cord_pa(ii,jj)-w_init1;
        dw2_init_cons(ii,jj)=avg_w2_cord_pa(ii,jj)-w_init2;
        dw1_init_risk(ii,jj)=avg_w1_cord_pa(ii,jj)-w_init1;
        dw2_init_risk(ii,jj)=avg_w2_cord_pa(ii,jj)-w_init2;
        %% Calculate change in beach widths with relation to reference value (Wr)
        dw1_ref_cord(ii,jj)=avg_w1_cord_pa(ii,jj)-ref_eff;
        dw2_ref_cord(ii,jj)=avg_w2_cord_pa(ii,jj)-ref_eff;
        dw1_ref_cons(ii,jj)=avg_w1_cons_pa(ii,jj)-ref_eff;
        dw2_ref_cons(ii,jj)=avg_w2_cons_pa(ii,jj)-ref_eff;
        dw1_ref_risk(ii,jj)=avg_w1_risk_pa(ii,jj)-ref_eff;
        dw2_ref_risk(ii,jj)=avg_w2_risk_pa(ii,jj)-ref_eff;
        %% Calculate Efficiencies
        % Calculate Efficiency for each model run

        %% Coordination
        % Community 1 
        if w1_nrsh_cord_pa(ii,jj)==0 && dw1_init_cord(ii,jj) <=0
            % There is no efficiency if a community does not nourish and loses beach width
            E1_cord(ii,jj)=NaN; 
        elseif w1_nrsh_cord_pa(ii,jj)==0 && dw1_init_cord(ii,jj) >0
            % Placeholder - will replace when max efficiency values are calculated
            E1_cord(ii,jj)=-250; 
        else
            % Calculate efficiency E = (delta_w + Wr)/(Nourish_W + Wr)
            E1_cord(ii,jj)=(dw1_init_cord(ii,jj)+ref_eff)/(w1_nrsh_cord_pa(ii,jj)+ref_eff);
        end

        % Community 2
        if w2_nrsh_cord_pa(ii,jj)==0 && dw2_init_cord(ii,jj) <=0
            % There is no efficiency if a community does not nourish and loses beach width
            E2_cord(ii,jj)=NaN;
        elseif w2_nrsh_cord_pa(ii,jj)==0 && dw2_init_cord(ii,jj) >0
            % Placeholder - will replace when max efficiency values are calculated
            E2_cord(ii,jj)=-250; 
        else
            % Calculate efficiency E = (delta_w + Wr)/(Nourish_W + Wr)
            E2_cord(ii,jj)=(dw2_init_cord(ii,jj)+ref_eff)/(w2_nrsh_cord_pa(ii,jj)+ref_eff);
        end

        %% Conservative Non-Coordination
        % Community 1 
        if w1_nrsh_cons_pa(ii,jj)==0 && dw1_init_cons(ii,jj) <=0
            % There is no efficiency if a community does not nourish and loses beach width
            E1_cons(ii,jj)=NaN;
        elseif w1_nrsh_cons_pa(ii,jj)==0 && dw1_init_cons(ii,jj) >0
            % Placeholder - will replace when max efficiency values are calculated
            E1_cons(ii,jj)=-250; 
        else
            % Calculate efficiency E = (delta_w + Wr)/(Nourish_W + Wr)
            E1_cons(ii,jj)=(dw1_init_cons(ii,jj)+ref_eff)/(w1_nrsh_cons_pa(ii,jj)+ref_eff);
        end

        % Community 2
        if w2_nrsh_cons_pa(ii,jj)==0 && dw2_init_cons(ii,jj) <=0
            % There is no efficiency if a community does not nourish and loses beach width
            E2_cons(ii,jj)=NaN;
        elseif w2_nrsh_cons_pa(ii,jj)==0 && dw2_init_cons(ii,jj) >0
            % Placeholder - will replace when max efficiency values are calculated
            E2_cons(ii,jj)=-250; 
        else
            % Calculate efficiency E = (delta_w + Wr)/(Nourish_W + Wr)
            E2_cons(ii,jj)=(dw2_init_cons(ii,jj)+ref_eff)/(w2_nrsh_cons_pa(ii,jj)+ref_eff);
        end

        %% Risky Non-Coordination
        % Community 1 
        if w1_nrsh_risk_pa(ii,jj)==0 && dw1_init_risk(ii,jj) <=0
            % There is no efficiency if a community does not nourish and loses beach width
            E1_risk(ii,jj)=NaN;
        elseif w1_nrsh_risk_pa(ii,jj)==0 && dw1_init_risk(ii,jj) >0
            % Placeholder - will replace when max efficiency values are calculated
            E1_risk(ii,jj)=-250; 
        else
            % Calculate efficiency E = (delta_w + Wr)/(Nourish_W + Wr)
            E1_risk(ii,jj)=(dw1_init_risk(ii,jj)+ref_eff)/(w1_nrsh_risk_pa(ii,jj)+ref_eff);
        end

        % Community 2
        if w2_nrsh_risk_pa(ii,jj)==0 && dw2_init_risk(ii,jj) <=0
            % There is no efficiency if a community does not nourish and loses beach width
            E2_risk(ii,jj)=NaN;
        elseif w2_nrsh_risk_pa(ii,jj)==0 && dw2_init_risk(ii,jj) >0
            % Placeholder - will replace when max efficiency values are calculated
            E2_risk(ii,jj)=-250; 
        else
            % Calculate efficiency E = (delta_w + Wr)/(Nourish_W + Wr)
            E2_risk(ii,jj)=(dw2_init_risk(ii,jj)+ref_eff)/(w2_nrsh_risk_pa(ii,jj)+ref_eff);
        end
    end
end

%% Find max efficiency to calculate passive efficiency (110% of max efficiency)
max_E1_cord=max(max(E1_cord));
max_E2_cord=max(max(E2_cord));
max_E1_cons=max(max(E1_cons));
max_E2_cons=max(max(E2_cons));
max_E1_risk=max(max(E1_risk));
max_E2_risk=max(max(E2_risk));
max_eff=max([max_E1_cord,max_E2_cord,max_E1_cons,max_E2_cons,max_E1_risk,max_E2_risk]);
passive_eff=max_eff*1.1;
%% Replace Placeholder Values with passive efficiency calculated above
E1_cord(E1_cord==-250)=passive_eff;
E2_cord(E2_cord==-250)=passive_eff;
E1_cons(E1_cons==-250)=passive_eff;
E2_cons(E2_cons==-250)=passive_eff;
E1_risk(E1_risk==-250)=passive_eff;
E2_risk(E2_risk==-250)=passive_eff;

%% Second loops
for ii=1:nn 
    for jj=1:mm
        %% Calculate the difference of efficiencies
        % Relative to other community (comm1-comm2)
        E1_E2_cord_diff(ii,jj)=E1_cord(ii,jj)-E2_cord(ii,jj);
        E1_E2_cons_diff(ii,jj)=E1_cons(ii,jj)-E2_cons(ii,jj);
        E1_E2_risk_diff(ii,jj)=E1_risk(ii,jj)-E2_risk(ii,jj);
        % Relative to Non-Coordination (Coordination - Non-Coordination)
        E1_cord_cons_diff(ii,jj)=E1_cord(ii,jj)-E1_cons(ii,jj);
        E2_cord_cons_diff(ii,jj)=E2_cord(ii,jj)-E2_cons(ii,jj);
        E1_cord_risk_diff(ii,jj)=E1_cord(ii,jj)-E1_risk(ii,jj);
        E2_cord_risk_diff(ii,jj)=E2_cord(ii,jj)-E2_risk(ii,jj);
        %% Calculate Difference in Nourishment Volumes
        % Relative to other community (comm1-comm2)
        V1_V2_cord_diff(ii,jj)=V1_nrsh_cord_pa(ii,jj)-V2_nrsh_cord_pa(ii,jj);
        V1_V2_cons_diff(ii,jj)=V1_nrsh_cons_pa(ii,jj)-V2_nrsh_cons_pa(ii,jj);
        V1_V2_risk_diff(ii,jj)=V1_nrsh_risk_pa(ii,jj)-V2_nrsh_risk_pa(ii,jj);
        % Relative to Non-Coordination
        V1_cord_cons_diff(ii,jj)=V1_nrsh_cord_pa(ii,jj)-V1_nrsh_cons_pa(ii,jj);
        V2_cord_cons_diff(ii,jj)=V2_nrsh_cord_pa(ii,jj)-V2_nrsh_cons_pa(ii,jj);
        V1_cord_risk_diff(ii,jj)=V1_nrsh_cord_pa(ii,jj)-V1_nrsh_risk_pa(ii,jj);
        V2_cord_risk_diff(ii,jj)=V2_nrsh_cord_pa(ii,jj)-V2_nrsh_risk_pa(ii,jj);
        
        %% Define Binaries: Volume Relative to Other Community
        % Coordination
        if V1_V2_cord_diff(ii,jj) > 0 % Community 1 Nourishes More
            V1_V2_binary_cord(ii,jj) = 1;
        elseif V1_V2_cord_diff(ii,jj) < 0 % Community 2 Nourishes More
            V1_V2_binary_cord(ii,jj) = 2;
        elseif V1_V2_cord_diff(ii,jj) == 0 % Nourish the same
            V1_V2_binary_cord(ii,jj) = 1.5;
        end
        % Non-Coordination - Conservative
        if V1_V2_cons_diff(ii,jj) > 0 % Community 1 Nourishes More
            V1_V2_binary_cons(ii,jj) = 1;
        elseif V1_V2_cons_diff(ii,jj) < 0 % Community 2 Nourishes More
            V1_V2_binary_cons(ii,jj) = 2;
        elseif V1_V2_cons_diff(ii,jj) == 0 % Same
            V1_V2_binary_cons(ii,jj) = 1.5;
        end
        % Non-Coordination - Risky
        if V1_V2_risk_diff(ii,jj) > 0 % Community 1 Nourishes More
            V1_V2_binary_risk(ii,jj) = 1;
        elseif V1_V2_risk_diff(ii,jj) < 0 % Community 2 Nourishes More
            V1_V2_binary_risk(ii,jj) = 2;
        elseif V1_V2_risk_diff(ii,jj) == 0 % Same
            V1_V2_binary_risk(ii,jj) = 1.5;
        end
        %% Define Binaries: Efficiency Relative to Other Community 
        % Coordination
        if E1_E2_cord_diff(ii,jj) > 0 % Community 1 is more efficient
            E1_E2_binary_cord(ii,jj) = 1;
        elseif E1_E2_cord_diff(ii,jj) < 0 % Community 2 is more efficient
            E1_E2_binary_cord(ii,jj) = 2;
        elseif E1_E2_cord_diff(ii,jj) == 0 % Same
            E1_E2_binary_cord(ii,jj) = 1.5;
        end
        % Non-Coordination - Conservative
        if E1_E2_cons_diff(ii,jj) > 0 % Community 1 is more efficient
            E1_E2_binary_cons(ii,jj) = 1;
        elseif E1_E2_cons_diff(ii,jj) < 0 % Community 2 is more efficient
            E1_E2_binary_cons(ii,jj) = 2;
        elseif E1_E2_cons_diff(ii,jj) == 0 % Same
            E1_E2_binary_cons(ii,jj) = 1.5;
        end
        % Non-Coordination - Risky
        if E1_E2_risk_diff(ii,jj) > 0 % Community 1 is more efficient
            E1_E2_binary_risk(ii,jj) = 1;
        elseif E1_E2_risk_diff(ii,jj) < 0 % Community 2 is more efficient
            E1_E2_binary_risk(ii,jj) = 2;
        elseif E1_E2_risk_diff(ii,jj) == 0 % Same
            E1_E2_binary_risk(ii,jj) = 1.5;
        end
        %% Categorical Model / Field Study Analysis - Relationship between communities
        % Coordination
        % Com 1 Nourishes More & Com 1 more efficient
        if V1_V2_binary_cord(ii,jj)== 1 && E1_E2_binary_cord(ii,jj) == 1 
            C1_C2_EV_cord(ii,jj)=1; 
        % Com 1 Nourishes More & Com 2 more efficient
        elseif V1_V2_binary_cord(ii,jj)== 1 && E1_E2_binary_cord(ii,jj) == 2 
            C1_C2_EV_cord(ii,jj)=2; 
        % Com 2 Nourishes More & Com 1 more efficient
        elseif V1_V2_binary_cord(ii,jj)== 2 && E1_E2_binary_cord(ii,jj) == 1 
            C1_C2_EV_cord(ii,jj)=3; 
        % Com 2 Nourishes More & Com 2 more efficient
        elseif V1_V2_binary_cord(ii,jj)== 2 && E1_E2_binary_cord(ii,jj) == 2 
            C1_C2_EV_cord(ii,jj)=4; 
        end
        % Conservative Non-Coordination
        % Com1 Nourish More & Com 1 more efficient
        if V1_V2_binary_cons(ii,jj)== 1 && E1_E2_binary_cons(ii,jj) == 1 
            C1_C2_EV_cons(ii,jj)=1; 
        % Com1 Nourish More & Com 2 more efficient
        elseif V1_V2_binary_cons(ii,jj)== 1 && E1_E2_binary_cons(ii,jj) == 2 
            C1_C2_EV_cons(ii,jj)=2; 
        % Com2 Nourish More & Com 1 more efficient
        elseif V1_V2_binary_cons(ii,jj)== 2 && E1_E2_binary_cons(ii,jj) == 1 
            C1_C2_EV_cons(ii,jj)=3; 
        % Com2 Nourish More & Com 2 more efficient
        elseif V1_V2_binary_cons(ii,jj)== 2 && E1_E2_binary_cons(ii,jj) == 2 
            C1_C2_EV_cons(ii,jj)=4; 
        end
        % Risky Non-Coordination
        % Com1 Nourish More & Com 1 more efficient
        if V1_V2_binary_risk(ii,jj)== 1 && E1_E2_binary_risk(ii,jj) == 1 
            C1_C2_EV_risk(ii,jj)=1; 
        % Com1 Nourish More & Com 2 more efficient
        elseif V1_V2_binary_risk(ii,jj)== 1 && E1_E2_binary_risk(ii,jj) == 2 
            C1_C2_EV_risk(ii,jj)=2; 
        % Com2 Nourish More & Com 1 more efficient
        elseif V1_V2_binary_risk(ii,jj)== 2 && E1_E2_binary_risk(ii,jj) == 1 
            C1_C2_EV_risk(ii,jj)=3; 
        % Com2 Nourish More & Com 2 more efficient
        elseif V1_V2_binary_risk(ii,jj)== 2 && E1_E2_binary_risk(ii,jj) == 2 
            C1_C2_EV_risk(ii,jj)=4; 
        end
        %% Define Binaries: Volume Relative to Non-Coordination
        % From maincode: volume diff = coordination - noncoordination
        % Community 1
        if V1_cord_cons_diff(ii,jj) == 0 % Nourish the same
            V1_binary_cord_cons(ii,jj)=NaN;
        end
        if V1_cord_cons_diff(ii,jj) > 0 % Undernourish in Non-Coord
            V1_binary_cord_cons(ii,jj)=1;
        end 
        if V1_cord_cons_diff(ii,jj) < 0 % Over-Nourish in Non-Coord
            V1_binary_cord_cons(ii,jj)=-1;
        end
        % Community 2
        if V2_cord_cons_diff(ii,jj) == 0 % Nourish the same
            V2_binary_cord_cons(ii,jj)=NaN;
        end
        if V2_cord_cons_diff(ii,jj) > 0 % Undernourish in Non-Coord
            V2_binary_cord_cons(ii,jj)=1;
        end 
        if V2_cord_cons_diff(ii,jj) < 0 % Over-Nourish in Non-Coord
            V2_binary_cord_cons(ii,jj)=-1;
        end

        %% Define Binaries: Efficiency Comparisons Relative to Non-Coordination
        % From maincode: Efficiency Diff = Coordination - NonCoordination
        % Community 1
        if E1_cord_cons_diff(ii,jj) > 0 % Coordination more efficient
            E1_binary_cord_cons(ii,jj) = 1;
        elseif E1_cord_cons_diff(ii,jj) < 0 % Non-Coordination more efficient 
            E1_binary_cord_cons(ii,jj) = -1;
        elseif E1_cord_cons_diff(ii,jj) == 0
            E1_binary_cord_cons(ii,jj) = NaN;
        end
        % Community 2
        if E2_cord_cons_diff(ii,jj) > 0 % Coordination more efficient
            E2_binary_cord_cons(ii,jj) = 1;
        elseif E2_cord_cons_diff(ii,jj) < 0 % Non-Coordination more efficient
            E2_binary_cord_cons(ii,jj) = -1;
        elseif E2_cord_cons_diff(ii,jj) == 0 % Equal
            E2_binary_cord_cons(ii,jj) = NaN;
        end

        %% Categorical Four Options per community: Relative to Non-Coordination
        % Community 1 - all relative to coordination
        % Overnourish & less efficient
        if V1_binary_cord_cons(ii,jj)==1 && E1_binary_cord_cons(ii,jj) == 1 
            E1_V1_cord_cons(ii,jj)=4; 
        % Overnourish & more efficient
        elseif V1_binary_cord_cons(ii,jj)==1 && E1_binary_cord_cons(ii,jj) == -1 
            E1_V1_cord_cons(ii,jj)=2; 
        % Undernourish & less efficient
        elseif V1_binary_cord_cons(ii,jj)==-1 && E1_binary_cord_cons(ii,jj) == 1 
            E1_V1_cord_cons(ii,jj)=3; 
        % Undernourish & more efficient
        elseif V1_binary_cord_cons(ii,jj)==-1 && E1_binary_cord_cons(ii,jj) == -1 
            E1_V1_cord_cons(ii,jj)=1; 
        end
        % Community 2 - all relative to coordination
        % Overnourish & less efficient
        if V2_binary_cord_cons(ii,jj)==1 && E2_binary_cord_cons(ii,jj) == 1 
            E2_V2_cord_cons(ii,jj)=4; 
        % Overnourish & more efficient
        elseif V2_binary_cord_cons(ii,jj)==1 && E2_binary_cord_cons(ii,jj) == -1 
            E2_V2_cord_cons(ii,jj)=2; 
        % Undernourish & less efficient
        elseif V2_binary_cord_cons(ii,jj)==-1 && E2_binary_cord_cons(ii,jj) == 1 
            E2_V2_cord_cons(ii,jj)=3; 
        % Undernourish & more efficient
        elseif V2_binary_cord_cons(ii,jj)==-1 && E2_binary_cord_cons(ii,jj) == -1 
            E2_V2_cord_cons(ii,jj)=1; 
        end
    end
end

%% Figures 2: Categorical Analysis
% clf(12)
figure (12) %figure 5e in paper
subplot(1,2,1)
sgtitle({'Efficiency and Nourishment'},'FontSize',sg_font,'FontWeight','bold')
pcolor(C1_vec/1e6,C2_vec/1e6,C1_C2_EV_cord)
colormap(opeq_color)
shading flat
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
hold on
p4 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#5D3A9B','MarkerEdgeColor','none');
p3 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#c2b5da','MarkerEdgeColor','none');
p2 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#f6c49f','MarkerEdgeColor','none');
p1 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#E66100','MarkerEdgeColor','none');
legend([p1, p2,p3,p4],{'Community 1 nourishes more & is more efficient',...
    'Community 1 nourishes more & is less efficient','Community 2 nourishes more & is less efficient',...
    'Community 2 nourishes more & is more efficient'},Location='southoutside',FontSize=legendboxes,Box='off')
pbaspect([1 1 1])
title('Coordination','FontSize',title_font,'FontWeight','normal')
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
clim([1 4])
% c1= colorbar;

subplot(1,2,2)
pcolor(C1_vec/1e6,C2_vec/1e6,C1_C2_EV_cons)
colormap(opeq_color)
shading flat
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
hold on
p4 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#5D3A9B','MarkerEdgeColor','none');
p3 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#c2b5da','MarkerEdgeColor','none');
p2 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#f6c49f','MarkerEdgeColor','none');
p1 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#E66100','MarkerEdgeColor','none');
legend([p1, p2,p3,p4],{'Community 1 nourishes more & is more efficient',...
    'Community 1 nourishes more & is less efficient','Community 2 nourishes more & is less efficient',...
    'Community 2 nourishes more & is more efficient'},Location='southoutside',FontSize=legendboxes,Box='off')
pbaspect([1 1 1])
title('Non-Coordination','FontSize',title_font,'FontWeight','normal')
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
clim([1 4])
% c2= colorbar;

%% Figures 3: Difference in Efficiency
max_cord1=max(max(E1_E2_cord_diff));
max_cord2=max(max(abs(E1_E2_cord_diff)));
cord_lims=max(max_cord1,max_cord2);
max_cons1=max(max(E1_E2_cons_diff));
max_cons2=max(max(abs(E1_E2_cons_diff)));
cons_lims=max(max_cons1,max_cons2);
lim_both=max(cord_lims*1.05,cons_lims*1.05);

figure (13) %figure 5e in paper
subplot(1,2,1)
sgtitle({'Difference in Efficiency';'Community 1 - Community 2'},'FontSize',sg_font,'FontWeight','normal','fontname','Arial')
pcolor(C1_vec/1e6,C2_vec/1e6,E1_E2_cord_diff)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
title('Coordination','FontSize',title_font,'FontWeight','bold')
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
clim([-1*lim_both lim_both])
pbaspect([1 1 1])
shading flat
colormap(opeq_colorrev)
colorbar('southoutside')

subplot(1,2,2)
pcolor(C1_vec/1e6,C2_vec/1e6,E1_E2_cons_diff)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
title('Non-Coordination','FontSize',title_font,'FontWeight','bold')
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
clim([-1*lim_both lim_both])
colormap(opeq_colorrev)
colorbar('southoutside')




