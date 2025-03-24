%{
This code is used to create subplots that show individual community beach
widths and net benefits line charts based on .mat model run data. 

To run this code, you will need to select your data file from the
"mat_data_files" variable using the "file_number" variable. You can also
add your own output data from the parametric_analyses.m file.

You will also need to select the the community baseline property values you
would like to plot from the C1_vec and C2_vec variables in your data file.

This code also calculates the efficiency based on the files available in
the "mat_data_files" variable. 
%}

%% Data Files
dir_loc='mat-data/'; % the subdirectory where .mat files are located
% Place your mat files in the efficiency_files variable as a vector
% Vector of data files to be studied
mat_data_files=["parametric_analysis_cs_bs_v1.mat",... % #1
    "parametric_analysis_cs_ba_v1.mat",... % #2
    "parametric_analysis_cs_ba_v2.mat",... % #3
    "parametric_analysis_ca_bs_v2.mat",... % #4
    "parametric_analysis_ca_ba_v1.mat"]; % #5 Vector of data files to be studied

%% IMPORTANT Pick your file name by changing the variable below
file_number=2;

%% IMPORTANT: Pick your community proprty values below
% This data needs to be determined from the PV matrices (C1_vec and C2_vec)
% Load your matrices first either by running this code up to "Load the
% Dataset" ~line40 or by adding your data to the Workspace manually
C1 = 9; % Community 1
C2 = 13; % Community 2

%% Determine Reference Value from Data Sources
nFiles = length(mat_data_files); % Number of data sources being evaluated for entire study
filenames=string.empty(0,nFiles); % Will hold strings of your files names
gama_data=NaN(nFiles); % Container to house gamma values from each data file
tmax_data=NaN(nFiles); % Container to house tmax data from each data file
w_init_data=NaN(nFiles); % Container to house w_init_data from each data file
for iData=1:nFiles    
    filenames(iData) = string(strcat(dir_loc,mat_data_files(iData))); % Get the file name and make sure its a string
    load(filenames(iData),"gamma_array","tmax","w_init"); % Load the datafile
    gama_data(iData)=max(gamma_array); % Gets the maximum gamma value for that dataset (is an array in the dataset)
    tmax_data(iData)=tmax; % Gets the tmax value for that dataset (one value in dataset)
    w_init_data(iData)=min(w_init); % Gets the min w_init for that dataset (is an array in the dataset)
end
ref_gamma=max(max(gama_data)); % Finds the max gamma value from all datasets
ref_tmax=max(max(tmax_data)); % Finds the max time horizon for all datasets
ref_w_init=min(min(w_init_data)); % Finds the minimum initial width for all datasets
ref_eff = abs(ref_w_init - ref_gamma * ref_tmax); % Caculates the reference point for the efficiency metric

%% Load the dataset
load(filenames(file_number)); 

%% Run the internal maincode Function for Coordination
[sr_cd_w,sr_cd_nb,tnb_t_cd,sr_cd_Eff_Theor1,sr_cd_Eff_Theor2,sr_Ben_cd,sr_Ben_disc_cd,sr_TBen_cd,sr_TC_cd]=...
    model_run(R1_cord_pa(C1,C2),R2_cord_pa(C1,C2),C1_vec(C1),C2_vec(C2),...
    nt,m,x_lot,w_init,D,Ka,Kc,theta_eq,beta1,beta2,disc,...
    dt,s,xN1,xN2,t,gamma,gamma_g4,gamma_g3,c,phi,tmax,...
    nlots_along1,nlots_along2,ref_eff);

%% Run the internal maincode Function for Conservative Non-Coordination
[sr_nc_w,sr_nc_nb,tnb_t_nc,sr_nc_Eff_Theor1,sr_nc_Eff_Theor2,sr_Ben_nc,sr_Ben_disc_nc,sr_TBen_nc,sr_TC_nc]=...
    model_run(R1_cons_pa(C1,C2),R2_cons_pa(C1,C2),C1_vec(C1),C2_vec(C2),...
    nt,m,x_lot,w_init,D,Ka,Kc,theta_eq,beta1,beta2,disc,...
    dt,s,xN1,xN2,t,gamma,gamma_g4,gamma_g3,c,phi,tmax,...
    nlots_along1,nlots_along2,ref_eff);

prop_line=zeros(nt,1); %use these to make horizontal lines
%% Plots: Used to make plots 
% find y limits 
% max nb
M_cd_nb=max(max(sr_cd_nb(:)));
M_nc_nb=max(max(sr_nc_nb(:)));
M_tnb_t_cd=max(max(tnb_t_cd(:)));
M_tnb_t_nc=max(max(tnb_t_nc(:)));
m_vec_nb=[M_cd_nb M_nc_nb M_tnb_t_cd M_tnb_t_nc];
Max_nb_x=max(m_vec_nb);
Max_nb=Max_nb_x+abs(Max_nb_x)*.1;
% min nb
M_cd_nb=min(min(sr_cd_nb(:)));
M_nc_nb=min(min(sr_nc_nb(:)));
m_vec_nb=[M_cd_nb M_nc_nb];
Min_nb_x=min(m_vec_nb);
Min_nb=Min_nb_x-abs(Max_nb_x)*.1;
% max w
M_cd_w=max(max(sr_cd_w(:)));
M_nc_w=max(max(sr_nc_w(:)));
m_vec_w=[M_cd_w M_nc_w];
Max_w_x=max(m_vec_w);
Max_w=Max_w_x+abs(Max_w_x)*.1;
% min w
M_cd_w=min(min(sr_cd_w(:)));
M_nc_w=min(min(sr_nc_w(:)));
m_vec_w=[M_cd_w M_nc_w];
Min_w_x=min(m_vec_w);
Min_w=Min_w_x-abs(Max_w_x)*.1;
% xlims
xlimvals=[0 tmax];


%% Plot colors (Beach Widths and Net Benefits)
% Colors
ghost1='#525151';
com1='#3e7fec';
com2='#ec3e7f';
ghost2='#bfbebe';
tnb_color='#656565';

%% Plot Legend Labels
pBPV1=num2str(floor(C1_vec(C1)));
pBPV2=num2str(floor(C2_vec(C2)));
ef_round=4;
r_round=1;
if isnan(R1_cord_pa(C1,C2))
    Rlabel1 = sprintf("R1: None - Efficiency:"+" "+round(sr_cd_Eff_Theor1,ef_round));
else
    Rlabel1 = sprintf('R1: %.3g - Efficiency: %2.11g',round(R1_cord_pa(C1,C2),r_round),round(sr_cd_Eff_Theor1,ef_round));
end
if isnan(R2_cord_pa(C1,C2))
    Rlabel2 = sprintf("R2: None - Efficiency:"+" "+round(sr_cd_Eff_Theor2,ef_round));
else
    Rlabel2 = sprintf('R2: %.3g - Efficiency: %2.11g',round(R2_cord_pa(C1,C2),r_round),round(sr_cd_Eff_Theor2,ef_round));
end
PVlabel1 = strcat('BPV1: $',pBPV1);
PVlabel2 = strcat('BPV2: $',pBPV2);
TNBlabel = ('Both Communities');
if isnan(R1_cons_pa(C1,C2))
    Rlabel1_nc=sprintf("R1: None - Efficiency:"+" "+round(sr_nc_Eff_Theor1,ef_round));
else
    Rlabel1_nc = sprintf('R1: %.3g - Efficiency: %2.11g',round(R1_cons_pa(C1,C2),r_round),round(sr_nc_Eff_Theor1,ef_round));
end
if isnan(R2_cons_pa(C1,C2))
    Rlabel2 = sprintf("R2: None - Efficiency:"+" "+round(sr_nc_Eff_Theor2,ef_round));
else
    Rlabel2_nc = sprintf('R2: %.3g - Efficiency: %2.11g',round(R2_cons_pa(C1,C2),r_round),round(sr_nc_Eff_Theor2,ef_round));
end

%% Plots
figure(16)
tiledlayout(2,2)
subplot(2,2,1)
hold on
box on
% axis square
% sgtitle('Two Community System','FontSize',20,'FontWeight','bold')
plot(t,sr_cd_w(:,1),color=com1,LineWidth=2) % Community 1
plot(t,sr_cd_w(:,2),color=com2,LineWidth=2) % Community 2
% plot(t,sr_cd_w(:,3),color=ghost1,LineWidth=2) % Ghost 1
% plot(t,sr_cd_w(:,4),color=ghost2,LineWidth=2) % Ghost 2
plot(t,prop_line(:,1),color='#7f7e7e',LineWidth=1,LineStyle='--')
ylim([Min_w Max_w])
xlim(xlimvals)
title({'Coordination';'Beach Width'},'FontSize',16,'FontWeight','normal');
xlabel('Years','FontSize',16,'FontWeight','normal');
ylabel('Beach Width (m)','FontSize',16,'FontWeight','normal');
legend({Rlabel1,Rlabel2},'Location','best')
set(gca,'XGrid','off','YGrid','on')

subplot(2,2,3)
hold on
box on
% axis square
plot(t,sr_cd_nb(:,1),color=com1,LineWidth=2) %use to plot NB
plot(t,sr_cd_nb(:,2),color=com2,LineWidth=2) %use to plot NB
plot(t,tnb_t_cd(:,1),color=tnb_color,LineWidth=2) %use to plot TNB
ylim([Min_nb Max_nb])
xlim(xlimvals)
title({'Coordation';'Net Benefits'},'FontSize',16,'FontWeight','normal');
xlabel('Years','FontSize',16,'FontWeight','normal');
ylabel('Net Benefits ($10M)','FontSize',16,'FontWeight','normal');
legend({PVlabel1,PVlabel2,TNBlabel},'Location','southeast')
set(gca,'XGrid','off','YGrid','on')

subplot(2,2,2)
hold on
box on
% axis square
sgtitle('Two Community System','FontSize',20,'FontWeight','bold')
plot(t,sr_nc_w(:,1),color=com1,LineWidth=2) % Beach Widths Community 1
plot(t,sr_nc_w(:,2),color=com2,LineWidth=2) % Beach Widths Community 2
% plot(t,sr_nc_w(:,1),color=ghost1,LineWidth=2) % Beach Widths Ghost 1
% plot(t,sr_nc_w(:,4),color=ghost2,LineWidth=2) % Beach Widths Ghost 2
plot(t,prop_line(:,1),color='#7f7e7e',LineWidth=1,LineStyle='--')
title({'Non-Coordination';'Beach Width'},'FontSize',16,'FontWeight','normal');
ylim([Min_w Max_w])
xlim(xlimvals)
xlabel('Years','FontSize',16,'FontWeight','normal');
ylabel('Beach Width (m)','FontSize',16,'FontWeight','normal');
legend({Rlabel1_nc,Rlabel2_nc},'Location','best')
set(gca,'XGrid','off','YGrid','on')

subplot(2,2,4)
hold on
box on
% axis square
plot(t,sr_nc_nb(:,1),color=com1,LineWidth=2) %use to plot NB
plot(t,sr_nc_nb(:,2),color=com2,LineWidth=2) %use to plot NB
plot(t,tnb_t_nc(:,1),color=tnb_color,LineWidth=2) %use to plot TNB
ylim([Min_nb Max_nb])
xlim(xlimvals)
title({'Non-Coordination';'Net Benefits'},'FontSize',16,'FontWeight','normal');
xlabel('Years','FontSize',16,'FontWeight','normal');
ylabel('Net Benefits ($10M)','FontSize',16,'FontWeight','normal');
legend({PVlabel1,PVlabel2,TNBlabel},'Location','southeast')
set(gca,'XGrid','off','YGrid','on')

hold off
% end

%% Model Run Function
function [w,nb,tnb_t,Eff_Theor1,Eff_Theor2,Ben,Ben_disc,TBen,TC]=...
    model_run(R1,R2,BPV1,BPV2,...
    nt,m,x_lot,w_init,D,Ka,Kc,theta_eq,beta1,beta2,disc,...
    dt,s,xN1,xN2,t,gamma,gamma_g4,gamma_g3,c,phi,tmax,...
    nlots_along1,nlots_along2,ref_eff)
    
PV1 = BPV1; PV2 = BPV2; % Sets First iteration of PV to Baseline Property Value
% Check the data values
formatSpecBPV = ['\nBPV1: %1.f   R1: %2.i ' ...
    '\nBPV2: %3.f,   R2: %4.i\n'];
fprintf(formatSpecBPV,BPV1,R1,BPV2,R2);
% formatSpecR = 'R1: %1.f and R2: %2.f \n';
% fprintf(formatSpecR,R1,R2);

alpha1=nlots_along1*BPV1; alpha2=nlots_along2*BPV2;
%% Storage Containers for each model run
theta=zeros(nt,m); qL=zeros(nt,m); qC=zeros(nt,m); fS=zeros(nt,m); fT=zeros(nt,m); 
xS=zeros(nt,m); xT=zeros(nt,m); w=zeros(nt,m);
Ben=zeros(nt,m); Ben_disc=zeros(nt,m); TBen=zeros(nt,m);
C=zeros(nt,m); TC=zeros(nt,m);
nb=zeros(nt,m);   
tnb_t=zeros(nt,1);
%containers to calculate time-specific property values
PropVal=zeros(nt,m); 

%% Initial Shoreface Conditions
xS(1,:)=x_lot(:)+w_init(:); % Sets the initial Shoreline
xT(1,:)=xS(1,:)+(D/theta_eq); % Set initial shore toe location
theta(1,:)=D./(xT(1,:)-xS(1,:)); 
w(1,:) = xS(1,:)-x_lot(1,:); % Set initial Width

%% Initial Benefits / Costs
% PROPERTY VALUE
PropVal(1,1)=BPV1*((w(1,1))/w(1,1))^beta1; % Property Value Community 1
PropVal(1,2)=BPV2*((w(1,2))/w(1,2))^beta2; % Property Value Community 2
% BENEFITS 
Ben(1,1)=alpha1*disc*((w(1,1))/w(1,1))^beta1; % Benefits Community 1
Ben(1,2)=alpha2*disc*((w(1,2))/w(1,2))^beta2; % Benefits Community 2
Ben_disc(1,:)=Ben(1,:);
TBen(1,:)=Ben(1,:)*dt; 
nb(1,:)=TBen(1,:);
tnb_t(1,1)=nb(1,1)+nb(1,2); 
%% Input Physical Parameters %%                
% R1 = R1_vector; %Rotation Length
nVol1 = 0.5*s(1)*D*xN1; nVol2 = 0.5*s(2)*D*xN2; % Nourishment Volume for Community 1
k1=1; k2=1; % Nourishment Counter

%% Main Code %%
for i=1:nt-1                       
    %% Nourishment Initiation + Volume
    % ff1 and ff2 determine if nourishment events occur
    ff1=round(t(i+1)-k1*R1,4);
    if ff1==0 
        k1=k1+1; 
        nourish1 = xN1;
    else
        nourish1=0;
    end
    ff2=round(t(i+1)-k2*R2,4);
    if ff2==0              
        k2=k2+1; 
        nourish2 = xN2;
    else
        nourish2 = 0;
    end
        
    %% Fluxes (Along/Cross-shore) and Shoreface Dynamics
    % alonshore fluxes
    qL(i,1)=2*Ka*((xS(i,4)-xS(i,1))/(s(4)+s(1))); 
    qL(i,2)=2*Ka*((xS(i,1)-xS(i,2))/(s(1)+s(2)));
    qL(i,3)=2*Ka*((xS(i,2)-xS(i,3))/(s(2)+s(3)));
    qL(i,4)=0; % No sediment flux between the boundary cells
    %s Shoreface Slopes
    theta(i,:)=D./(xT(i,:)-xS(i,:)); % Calculate the shoreface slope for all cells
    % Cross-shore fluxes
    qC(i,:)=Kc.*(theta(i,:)-theta_eq); % Caculcate the cross-shore fluxes for all cells
    
    %% ODE's 
    % Community Shoreface Toe
    fT(i,1)=(4*qC(i,1))/D-gamma; % Community 1
    fT(i,2)=(4*qC(i,2))/D-gamma; % Community 2
    % Boundary Cells Shoreface Toe
    fT(i,3)=(4*qC(i,3))/D-gamma_g3; % Boundary Cell 3
    fT(i,4)=(4*qC(i,4))/D-gamma_g4; % Boundary Cell 4
    % Shoreline Changes Communities - new code
    fS(i,1)=(2/s(1))*(qL(i,1)-qL(i,2))-(4*qC(i,1))/D-gamma; % Community 1
    fS(i,2)=(2/s(2))*(qL(i,2)-qL(i,3))-(4*qC(i,2))/D-gamma; % Community 2
    % Boundary Cells - new code
    fS(i,3)=(2/s(3))*(qL(i,3)-qL(i,4))-(4*qC(i,3))/D-gamma_g3; % Boundary 3
    fS(i,4)=(2/s(4))*(qL(i,4)-qL(i,1))-(4*qC(i,4))/D-gamma_g4; % Boundary 4
    
    %% Check if community loses property and stops nourishing
    if xS(i,1) <= x_lot(1)
        nourish1=0;
    end
    if xS(i,2) <= x_lot(2)
        nourish2=0;
    end

    %% Numerical Approximations 
    % Community Toe Location
    xT(i+1,1)=xT(i,1)+dt*fT(i,1); 
    xT(i+1,2)=xT(i,2)+dt*fT(i,2);
    % Boundary Cell Toe Location
    xT(i+1,3)=xT(i,3)+dt*fT(i,3); 
    xT(i+1,4)=xT(i,4)+dt*fT(i,4);
    
    % Community Shoreline Location
    xS(i+1,1)=xS(i,1)+dt*fS(i,1)+nourish1;
    xS(i+1,2)=xS(i,2)+dt*fS(i,2)+nourish2;
    % Boundary Cell Shoreline Location
    xS(i+1,3)=xS(i,3)+dt*fS(i,3);
    xS(i+1,4)=xS(i,4)+dt*fS(i,4);
    
    % Beach Width - used for economics not phyiscal properties
    w(i+1,1)=xS(i+1,1)-x_lot(1);
    w(i+1,2)=xS(i+1,2)-x_lot(2);  
    w(i+1,3)=xS(i+1,3)-x_lot(3);
    w(i+1,4)=xS(i+1,4)-x_lot(4);

    %% Benefit      
    %Calculate benefits for community 1
    if w(i+1,1) < 0
        Ben(i+1,1)=0;
        PropVal(i+1,1)=0;
    else
        Ben(i+1,1)=alpha1*disc*((w(i+1,1))/w(1,1))^beta1;
        PropVal(i+1,1)=PV1*((w(i+1,1))/w(1,1))^beta1;
    end
    % Calculate benefits for community 2
    if w(i+1,2) < 0
        Ben(i+1,2)=0;
        PropVal(i+1,2)=0;
    else
        Ben(i+1,2)=alpha2*disc*((w(i+1,2))/w(1,2))^beta2;
        PropVal(i+1,2)=PV2*((w(i+1,2))/w(1,2))^beta2;
    end

    %% Benefits with Discounting / Total Benefits
    % Future discounting of benefits
    Ben_disc(i+1,1)=Ben(i+1,1)*exp(-disc*t(i+1));
    Ben_disc(i+1,2)=Ben(i+1,2)*exp(-disc*t(i+1)); 
    % Summation of benefits over time
    TBen(i+1,1)=dt*(sum(Ben_disc(:,1)));
    TBen(i+1,2)=dt*(sum(Ben_disc(:,2)));

    %% Cost 
    % if statement determines if nourishment & cost occurs 
    if nourish1==0
        C(i+1,1)=0;
    else
        C(i+1,1)=(c+nVol1*phi)*exp(-disc*t(i+1));
    end
    if nourish2==0
        C(i+1,2)=0;
    else
        C(i+1,2)=(c+nVol2*phi)*exp(-disc*t(i+1));
    end
    TC(i+1,1)=C(i+1,1)+TC(i,1);
    TC(i+1,2)=C(i+1,2)+TC(i,2);

    %% marginal net benefit 
    nb(i+1,1)=TBen(i+1,1)-TC(i+1,1);
    nb(i+1,2)=TBen(i+1,2)-TC(i+1,2);
    tnb_t(i+1,1)=nb(i+1,1)+nb(i+1,2);

            
end

%% Calculate Efficiencies for Plots
% Temporal smoothing time period
analysis_avg_yrs=1; % The number of years the program should look at for behaviors
analysis_yr_start = (tmax-analysis_avg_yrs)/dt; % find index number for the start of the behavior period 

% Conditional to look for no nourishing
avg_width1=mean(w(analysis_yr_start:end,1));
avg_width2=mean(w(analysis_yr_start:end,2));
delta_w1=avg_width1+abs(ref_eff);
delta_w2=avg_width2+abs(ref_eff);
if isnan(R1) && avg_width1 <= w_init(1)
    % Eff_original=NaN;
    w_added1=0;   % total beach width added over course of model run
    Eff_Theor1='None';  % Theoretical efficiency set to NaN if no sand added
elseif isnan(R1) && avg_width1 > w_init(1)
    w_added1=0;   % total beach width added over course of model run
    Eff_Theor1='Maximum';  % Theoretical efficiency set to NaN if no sand added
else
    w_added1=(k1-1)*xN1;  % total beach width added over course of model run
    Eff_Theor1=delta_w1/(w_added1+abs(ref_eff)); % Original efficiency metric   
end
if isnan(R2) && avg_width2 <= w_init(2)
    w_added2=0;
    Eff_Theor2='None';
elseif isnan(R2) && avg_width2 > w_init(2)
    w_added2=0;
    Eff_Theor2='Maximum';
else
    w_added2=(k2-1)*xN2;
    Eff_Theor2=delta_w2/(w_added2+abs(ref_eff)); % Original Efficiency Metric
end
eff_spec1 = 'Delta Width1: %1.0f - Width Added1: %2.f - Eff1: %0.03f \n';
fprintf(eff_spec1,delta_w1,w_added1,Eff_Theor1);
eff_spec2 = 'Delta Width2: %1.0f - Width Added2: %2.f - Eff2: %0.03f \n';
fprintf(eff_spec2,delta_w2,w_added2,Eff_Theor2);

% End of Function
end

