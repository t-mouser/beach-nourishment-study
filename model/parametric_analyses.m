%{
This code is used to create the regime diagram data files that are used to
plot behaviors, nourishment volumes, and nourishment efficiencies. It runs
the maincode function. This script specifically evaluates communities 
optimizing their nourishment schedules under coordination, conservative 
non-coordination, and risky non-cordination. 

When running this code, please feel free to change the parameter inputs and 
file save names to generate new analyses. 

In this iteration, two vectors are created to analyze how different
baseline property values impact community nourishment decisions and 
outcomes.

Additionally, the erosion rates (gamma values) for the boundary cells can
be changed to force different inherent efficiencies for their adjacent
community.
%}

%% Code Timer Start
disp(datetime)
tic
time_last=0; % Used to time the model run

%% Use to SAVE the model run
save_run = true;

% File naming convention: s (symmetric), % a (asymmetric), 
% c (coastal dynamics), b (beta values)
% Example: parametric_analysis_cs_ba = 
% parametric analysis, coastal symmetric and beta asymmetric
file_title = 'parametric_analysis_cs_ba_v1_test';
dir_loc = 'mat-data/';
FileName = strcat(dir_loc,file_title);

%% Input parameters

%% Time Parameters
tmax=50; dt=0.1; t=0:dt:tmax; nt=length(t); 
analysis_avg_yrs=1; % The number of years the program should look at for behaviors

%% Rotation Schedules
Rmin = 1; Rmax = 25; dR = .25; % Rotation Parameters
R1_vec = Rmin:dR:Rmax; R2_vec = Rmin:dR:Rmax; % Creates rotation vectors
R1_vector = [NaN R1_vec]; R2_vector = [NaN R2_vec]; % Adds NaN to vectors

%% Analysis Data
BPVmin=0.1e6; BPV_max=2e6; BPV_num=20; % Set property value vectors
C1_vec=linspace(BPVmin,BPV_max,BPV_num); % Property Value 1 Vector
C2_vec=linspace(BPVmin,BPV_max,BPV_num); % Property Value 2 Vector

%% Economic Parameters
% BPV1=1.800e6; BPV2 = 1.500e6; % Community Baseline Property Values 
beta1=0.6; beta2=0.4; % Hedonic Price Parameter
disc=0.06; %discount rate

%% Community/Ghost Parameters
% BoundaryCellStationary=1; % Holds the boundary cell in place
w_init1 = 100; % Initial Beach Width for Ghost Cell 1
w_init2 = 100; % Initial Beach Width for Community 1
w_init3 = 100; % Initial Beach Width for Community 2
w_init4 = 100; % Initial Beach Width for Ghost Cell 2
s1=2000; s2 = 2000; % Community Alongshore Length (m)
a_lot1=40; a_lot2=a_lot1; % Alongshore Lot dimensions(m)
boundary3=2000;boundary4=boundary3; %Ghost cell widths - two ghost cells

% These are calculated from above
w_init=[w_init1 w_init2 w_init3 w_init4]; %initial beach width
s=[s1 s2 boundary3 boundary4]; %alongshore compartment lengths (m)
m=length(s); %number of cells
x_lot=[a_lot1 a_lot2 a_lot1 a_lot1]; % Cross-shore Lot dimensions(m)
nlots_along1=floor(s1/a_lot1); nlots_along2=floor(s2/a_lot2); % number of alongshore houses
rows_cross=1; %# of cross-shore property rows 2
prop_line=zeros(nt,1); %use these to make horizontal lines

%% Shoreface Parameters
theta_eq=0.02; %equilibrium shoreface slope
gamma_g4=4; % Erosion rate for ghost cell 1
gamma=4; %erosion rate for both communities (m/yr) 
gamma_g3=4; % erosion rate for ghost cell 2
gamma_array=[gamma_g4, gamma, gamma_g3];
Ka=250e3; % alongshore flux coeff
Kc=1000; %cross-shore flux coeff
D=10; %depth of closure (m) NJ=16

%% Nourishment Parameters
phi=10; %sand cost per community ($/m^3)
c=.25e6; %fixed nourishment cost
xN1=50; xN2=xN1; %Nourishment extent seaward (m)
nourish1=0; nourish2=0; % sets baseline nourishment parameters for time zero (non-nourishment)

%% VARIABLE INPUT ENDS - SCRIPT RUNS BELOW

%% Counters for the multi-run analyses
nn=length(C1_vec); mm=length(C2_vec); %vector lengths

%% Variable Storage
% Matrices to store data for all the model runs
% _pa after each variable simply stands for parametric analysis
%       Abbreviations:
%       Coordination = cord; 
%       Conservative Non-Coordination = cons; 
%       Risky Non-Coordination = risk
R1_cord_pa=NaN(nn,mm); R2_cord_pa=NaN(nn,mm); % chosen rotation intervals coordination
R1_cons_pa=NaN(nn,mm); R2_cons_pa=NaN(nn,mm); % chosen rotation intervals conservative non-coordination
R1_risk_pa=NaN(nn,mm); R2_risk_pa=NaN(nn,mm); % chosen rotation intervals risky non-coordination
TNB_cord_pa=NaN(nn,mm); NB1_cord_pa=NaN(nn,mm); NB2_cord_pa=NaN(nn,mm); % Total net benefits coordination
TNB_cons_pa=NaN(nn,mm); NB1_cons_pa=NaN(nn,mm); NB2_cons_pa=NaN(nn,mm); % Total net benefits cons. non coord
TNB_risk_pa=NaN(nn,mm); NB1_risk_pa=NaN(nn,mm); NB2_risk_pa=NaN(nn,mm); % Total net benefits risky non coord
Beh1_cord_pa=NaN(nn,mm); Beh2_cord_pa=NaN(nn,mm); Beh_cord_pa=NaN(nn,mm); % Behaviors coordination
Beh1_cons_pa=NaN(nn,mm); Beh2_cons_pa=NaN(nn,mm); Beh_cons_pa=NaN(nn,mm); % Behaviors cons. non coord
Beh1_risk_pa=NaN(nn,mm); Beh2_risk_pa=NaN(nn,mm); Beh_risk_pa=NaN(nn,mm); % Behaviors risky non coord
avg_w1_cord_pa=NaN(nn,mm); avg_w2_cord_pa=NaN(nn,mm); % Average beach widths coordination
avg_w1_cons_pa=NaN(nn,mm); avg_w2_cons_pa=NaN(nn,mm); % Average beach widths cons non coord
avg_w1_risk_pa=NaN(nn,mm); avg_w2_risk_pa=NaN(nn,mm); % Average beach widths risky non coord
V1_nrsh_cord_pa=NaN(nn,mm); V2_nrsh_cord_pa=NaN(nn,mm); % Nourishment volumes coordination
V1_nrsh_cons_pa=NaN(nn,mm); V2_nrsh_cons_pa=NaN(nn,mm); % Nourishment volumes cons non coord
V1_nrsh_risk_pa=NaN(nn,mm); V2_nrsh_risk_pa=NaN(nn,mm); % Nourishment volumes risky non coord
w1_nrsh_cord_pa=NaN(nn,mm); w2_nrsh_cord_pa=NaN(nn,mm); % Nourishment width added coordination
w1_nrsh_cons_pa=NaN(nn,mm); w2_nrsh_cons_pa=NaN(nn,mm); % Nourishment width added cons non coord
w1_nrsh_risk_pa=NaN(nn,mm); w2_nrsh_risk_pa=NaN(nn,mm); % Nourishment width added risky non coord

%% Uses Function maincode from maincode.m 
for ii=1:nn % BPV1 Loop, stored in rows
    %% vectors for data
    % Stores the data for each vector to be added to the end matrix
    R1_cord_vec=NaN(1,mm); R2_cord_vec=NaN(1,mm);
    R1_cons_vec=NaN(1,mm); R2_cons_vec=NaN(1,mm);
    R1_risk_vec=NaN(1,mm); R2_risk_vec=NaN(1,mm);
    TNB_cord_vec=NaN(1,mm); NB1_cord_vec=NaN(1,mm); NB2_cord_vec=NaN(1,mm);
    TNB_cons_vec=NaN(1,mm); NB1_cons_vec=NaN(1,mm); NB2_cons_vec=NaN(1,mm);
    TNB_risk_vec=NaN(1,mm); NB1_risk_vec=NaN(1,mm); NB2_risk_vec=NaN(1,mm);
    Beh1_cord_vec=NaN(1,mm); Beh2_cord_vec=NaN(1,mm); Beh_cord_vec=NaN(1,mm);
    Beh1_cons_vec=NaN(1,mm); Beh2_cons_vec=NaN(1,mm); Beh_cons_vec=NaN(1,mm);
    Beh1_risk_vec=NaN(1,mm); Beh2_risk_vec=NaN(1,mm); Beh_risk_vec=NaN(1,mm);
    avg_w1_cord_vec=NaN(1,mm); avg_w2_cord_vec=NaN(1,mm);
    avg_w1_cons_vec=NaN(1,mm); avg_w2_cons_vec=NaN(1,mm);
    avg_w1_risk_vec=NaN(1,mm); avg_w2_risk_vec=NaN(1,mm);
    V1_nrsh_cord_vec=NaN(1,mm); V2_nrsh_cord_vec=NaN(1,mm);
    V1_nrsh_cons_vec=NaN(1,mm); V2_nrsh_cons_vec=NaN(1,mm);
    V1_nrsh_risk_vec=NaN(1,mm); V2_nrsh_risk_vec=NaN(1,mm);
    w1_nrsh_cord_vec=NaN(1,mm); w2_nrsh_cord_vec=NaN(1,mm);
    w1_nrsh_cons_vec=NaN(1,mm); w2_nrsh_cons_vec=NaN(1,mm);
    w1_nrsh_risk_vec=NaN(1,mm); w2_nrsh_risk_vec=NaN(1,mm);

    %% Loop for vectors
    for jj=1:mm % BPV2 Loop, stored in Columns
        BPV1=C1_vec(ii); %
        BPV2=C2_vec(jj);%sets the gamma value for each loop
        alpha1=nlots_along1*BPV1; alpha2=nlots_along2*BPV2; % alpha (property value * # properties)
        %% maincode FUNCTION being called
        [R1_cord,R2_cord,TNB_cord,NB1_cord,NB2_cord,Beh1_cord,Beh2_cord,Beh_cord,...
        avg_w1_cord,avg_w2_cord,R1_cons,R2_cons,TNB_cons,NB1_cons,NB2_cons,Beh1_cons,...
        Beh2_cons,Beh_cons,avg_w1_cons,avg_w2_cons,R1_risk,R2_risk,TNB_risk,NB1_risk,...
        NB2_risk,Beh1_risk,Beh2_risk,Beh_risk,avg_w1_risk,avg_w2_risk,w1_nrsh_cord,...
        w2_nrsh_cord,w1_nrsh_cons,w2_nrsh_cons,w1_nrsh_risk,w2_nrsh_risk,V1_nrsh_cord,...
        V2_nrsh_cord,V1_nrsh_cons,V2_nrsh_cons,V1_nrsh_risk,V2_nrsh_risk]...
        =maincode(dt,nt,t,tmax,R1_vector,R2_vector,BPV1,BPV2,beta1,beta2,...
        disc,w_init,s,m,x_lot,alpha1,alpha2,theta_eq,gamma,gamma_g4,gamma_g3,Ka,Kc,D,...
        phi,c,xN1,xN2,analysis_avg_yrs);
        %% Vector Data Input (see variable storage section for descriptions)
        R1_cord_vec(jj)=R1_cord; R2_cord_vec(jj)=R2_cord;
        R1_cons_vec(jj)=R1_cons; R2_cons_vec(jj)=R2_cons;
        R1_risk_vec(jj)=R1_risk; R2_risk_vec(jj)=R2_risk;
        TNB_cord_vec(jj)=TNB_cord; NB1_cord_vec(jj)=NB1_cord; NB2_cord_vec(jj)=NB2_cord;
        TNB_cons_vec(jj)=TNB_cons; NB1_cons_vec(jj)=NB1_cons; NB2_cons_vec(jj)=NB2_cons;
        TNB_risk_vec(jj)=TNB_risk; NB1_risk_vec(jj)=NB1_risk; NB2_risk_vec(jj)=NB2_risk;        
        Beh1_cord_vec(jj)=Beh1_cord; Beh2_cord_vec(jj)=Beh2_cord; Beh_cord_vec(jj)=Beh_cord;
        Beh1_cons_vec(jj)=Beh1_cons; Beh2_cons_vec(jj)=Beh2_cons; Beh_cons_vec(jj)=Beh_cons;        
        Beh1_risk_vec(jj)=Beh1_risk; Beh2_risk_vec(jj)=Beh2_risk; Beh_risk_vec(jj)=Beh_risk;
        avg_w1_cord_vec(jj)=avg_w1_cord; avg_w2_cord_vec(jj)=avg_w2_cord;
        avg_w1_cons_vec(jj)=avg_w1_cons; avg_w2_cons_vec(jj)=avg_w2_cons;
        avg_w1_risk_vec(jj)=avg_w1_risk; avg_w2_risk_vec(jj)=avg_w2_risk;
        V1_nrsh_cord_vec(jj)=V1_nrsh_cord; V2_nrsh_cord_vec(jj)=V2_nrsh_cord;
        V1_nrsh_cons_vec(jj)=V1_nrsh_cons; V2_nrsh_cons_vec(jj)=V2_nrsh_cons;
        V1_nrsh_risk_vec(jj)=V1_nrsh_risk; V2_nrsh_risk_vec(jj)=V2_nrsh_risk;
        w1_nrsh_cord_vec(jj)=w1_nrsh_cord; w2_nrsh_cord_vec(jj)=w2_nrsh_cord;
        w1_nrsh_cons_vec(jj)=w1_nrsh_cons; w2_nrsh_cons_vec(jj)=w2_nrsh_cons;
        w1_nrsh_risk_vec(jj)=w1_nrsh_risk; w2_nrsh_risk_vec(jj)=w2_nrsh_risk;
        
        %% Time Counter Inner Loop
        formatSpec = ['\nRow %1.f - Column %2.f Loop Complete! ' ...
            '\n     Loop Time: %0.03g seconds ' ...
            '\n     Total Time:. %0.04g minutes\n'];
        
        total_time = round(toc/60,2);
        loop_time=round((toc/60-time_last)*60,2);
        fprintf(formatSpec,floor(ii),floor(jj),loop_time,total_time);
        time_last=total_time;
    end
    %% Adds each vector to the matrix to create the final regime matrices
    % These matrices are determined by the property value inputs in this
    % code
    R1_cord_pa(ii,:)=R1_cord_vec; R2_cord_pa(ii,:)=R2_cord_vec;
    R1_cons_pa(ii,:)=R1_cons_vec; R2_cons_pa(ii,:)=R2_cons_vec;
    R1_risk_pa(ii,:)=R1_risk_vec; R2_risk_pa(ii,:)=R2_risk_vec;
    TNB_cord_pa(ii,:)=TNB_cord_vec; NB1_cord_pa(ii,:)=NB1_cord_vec; NB2_cord_pa(ii,:)=NB2_cord_vec;
    TNB_cons_pa(ii,:)=TNB_cons_vec; NB1_cons_pa(ii,:)=NB1_cons_vec; NB2_cons_pa(ii,:)=NB2_cons_vec;
    TNB_risk_pa(ii,:)=TNB_risk_vec; NB1_risk_pa(ii,:)=NB1_risk_vec; NB2_risk_pa(ii,:)=NB2_risk_vec;    
    Beh1_cord_pa(ii,:)=Beh1_cord_vec; Beh2_cord_pa(ii,:)=Beh2_cord_vec; Beh_cord_pa(ii,:)=Beh_cord_vec;
    Beh1_cons_pa(ii,:)=Beh1_cons_vec; Beh2_cons_pa(ii,:)=Beh2_cons_vec; Beh_cons_pa(ii,:)=Beh_cons_vec;       
    Beh1_risk_pa(ii,:)=Beh1_risk_vec; Beh2_risk_pa(ii,:)=Beh2_risk_vec; Beh_risk_pa(ii,:)=Beh_risk_vec;
    avg_w1_cord_pa(ii,:)=avg_w1_cord_vec; avg_w2_cord_pa(ii,:)=avg_w2_cord_vec;
    avg_w1_cons_pa(ii,:)=avg_w1_cons_vec; avg_w2_cons_pa(ii,:)=avg_w2_cons_vec;
    avg_w1_risk_pa(ii,:)=avg_w1_risk_vec; avg_w2_risk_pa(ii,:)=avg_w2_risk_vec;
    V1_nrsh_cord_pa(ii,:)=V1_nrsh_cord_vec; V2_nrsh_cord_pa(ii,:)=V2_nrsh_cord_vec;
    V1_nrsh_cons_pa(ii,:)=V1_nrsh_cons_vec; V2_nrsh_cons_pa(ii,:)=V2_nrsh_cons_vec;
    V1_nrsh_risk_pa(ii,:)=V1_nrsh_risk_vec; V2_nrsh_risk_pa(ii,:)=V2_nrsh_risk_vec;
    w1_nrsh_cord_pa(ii,:)=w1_nrsh_cord_vec; w2_nrsh_cord_pa(ii,:)=w2_nrsh_cord_vec;
    w1_nrsh_cons_pa(ii,:)=w1_nrsh_cons_vec; w2_nrsh_cons_pa(ii,:)=w2_nrsh_cons_vec;
    w1_nrsh_risk_pa(ii,:)=w1_nrsh_risk_vec; w2_nrsh_risk_pa(ii,:)=w2_nrsh_risk_vec;

end

%% Save Data
if save_run == true
    save(FileName); %saves output data (from this external script only!) in your root directory
end
