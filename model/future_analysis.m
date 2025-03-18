%% For Processing Time Keeping Only
disp(datetime)
tic
time_last=0; % Used to time the model run

%% File Save Information - input your file name under FileName
% File naming convention: s (symmetric), % a (asymmetric), 
% c (coastal dynamics), b (beta values)
% Example: parametric_analysis_cs_ba = 
% parametric analysis, coastal symmetric and beta asymmetric

% Enter data here:
FileName = 'future_analysis_cs_ba_v1_dr'; % your file name without directory location
save_run = true; % Keep on true, saves the file
dir_loc = 'mat-data/'; % The subdirectory location
FullFileName = strcat(dir_loc,FileName); % The file to be saved

%% Input parameters

%% X and Y Parameters of Pseudoplot
% Choose the ranges below
gamma_min=4; gamma_max=16; gamma_num=20; % Range of gamma values
phi_min=10; phi_max=30;phi_num=gamma_num; % Range of sand costs

%Creates the test parameter vectors, do not change
C1_vec=linspace(gamma_min,gamma_max,gamma_num); % Erosion Rate Vector
C2_vec=linspace(phi_min,phi_max,phi_num); % Sand Cost Vector

%% Time Period
tmax=50; dt=0.1; t=0:dt:tmax; nt=length(t); 
analysis_avg_yrs=1; % The number of years the program should look at for behaviors

%% Rotation Schedules
Rmin = 1; Rmax = 25; dR = .25; % Rotation Parameters
R1_vec = Rmin:dR:Rmax; R2_vec = Rmin:dR:Rmax; % Creates rotation vectors
R1_vector = [NaN R1_vec]; R2_vector = [NaN R2_vec]; % Adds NaN to vectors

%% Economic Parameters
BPV1=1.3e6; BPV2 = 1.7e6; % Community Baseline Property Values 
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
ghost1=2000;ghost2=ghost1; %Ghost cell widths - two ghost cells

%% These are calculated from above
w_init=[w_init1 w_init2 w_init3 w_init4]; %initial beach width
s=[ghost1 s1 s2 ghost2]; %alongshore compartment lengths (m)
m=length(s); %number of cells
x_lot=[a_lot1 a_lot1 a_lot2 a_lot1]; % Cross-shore Lot dimensions(m)
nlots_along1=floor(s1/a_lot1); nlots_along2=floor(s2/a_lot2); % number of alongshore houses
rows_cross=1; %# of cross-shore property rows 2
prop_line=zeros(nt,1); %use these to make horizontal lines
alpha1=nlots_along1*BPV1; alpha2=nlots_along2*BPV2; % alpha (property value * # properties)

%% Shoreface Parameters
theta_eq=0.02; %equilibrium shoreface slope
% gamma_g1=4; % Erosion rate for ghost cell 1
% gamma=4; %erosion rate for both communities (m/yr) 
% gamma_g2=4; % erosion rate for ghost cell 2
gamma_array=[gamma_min, gamma_max];
Ka=250e3; % alongshore flux coeff
Kc=1000; %cross-shore flux coeff
D=10; %depth of closure (m) NJ=16

%% Nourishment Parameters
% phi=10; %sand cost per community ($/m^3)
c=.25e6; %fixed nourishment cost
xN1=50; xN2=xN1; %Nourishment extent seaward (m)
nourish1=0; nourish2=0; % sets baseline nourishment parameters for time zero (non-nourishment)

%% Code Storage and Processing

%% Output Storage
% Find vector lengths to make containers
nn=length(C1_vec); mm=length(C2_vec); %vector lengths
% Rotation Values
R1_cord_pa=NaN(nn,mm); R2_cord_pa=NaN(nn,mm);
R1_cons_pa=NaN(nn,mm); R2_cons_pa=NaN(nn,mm);
R1_risk_pa=NaN(nn,mm); R2_risk_pa=NaN(nn,mm);
% Net benefits
TNB_cord_pa=NaN(nn,mm); NB1_cord_pa=NaN(nn,mm); NB2_cord_pa=NaN(nn,mm);
TNB_cons_pa=NaN(nn,mm); NB1_cons_pa=NaN(nn,mm); NB2_cons_pa=NaN(nn,mm);
TNB_risk_pa=NaN(nn,mm); NB1_risk_pa=NaN(nn,mm); NB2_risk_pa=NaN(nn,mm);
% Behaviors
Beh1_cord_pa=NaN(nn,mm); Beh2_cord_pa=NaN(nn,mm); Beh_cord_pa=NaN(nn,mm);
Beh1_cons_pa=NaN(nn,mm); Beh2_cons_pa=NaN(nn,mm); Beh_cons_pa=NaN(nn,mm);
Beh1_risk_pa=NaN(nn,mm); Beh2_risk_pa=NaN(nn,mm); Beh_risk_pa=NaN(nn,mm);
% Beach Widths
avg_w1_cord_pa=NaN(nn,mm); avg_w2_cord_pa=NaN(nn,mm);
avg_w1_cons_pa=NaN(nn,mm); avg_w2_cons_pa=NaN(nn,mm);
avg_w1_risk_pa=NaN(nn,mm); avg_w2_risk_pa=NaN(nn,mm);
% Nourishment Volumes Added
V1_nrsh_cord_pa=NaN(nn,mm); V2_nrsh_cord_pa=NaN(nn,mm);
V1_nrsh_cons_pa=NaN(nn,mm); V2_nrsh_cons_pa=NaN(nn,mm);
V1_nrsh_risk_pa=NaN(nn,mm); V2_nrsh_risk_pa=NaN(nn,mm);
% Nourishment Width Added
w1_nrsh_cord_pa=NaN(nn,mm); w2_nrsh_cord_pa=NaN(nn,mm);
w1_nrsh_cons_pa=NaN(nn,mm); w2_nrsh_cons_pa=NaN(nn,mm);
w1_nrsh_risk_pa=NaN(nn,mm); w2_nrsh_risk_pa=NaN(nn,mm);

%% Maincode Function 
for ii=1:nn % BPV1 Loop, stored in rows
    %% Vector containers
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
        gamma=C1_vec(ii); gamma_g1=gamma; gamma_g2=gamma; %
        phi=C2_vec(jj);%sets the gamma value for each loop
        
        %% Maincode Function Being Called
        [R1_cord,R2_cord,TNB_cord,NB1_cord,NB2_cord,Beh1_cord,Beh2_cord,Beh_cord,...
        avg_w1_cord,avg_w2_cord,R1_cons,R2_cons,TNB_cons,NB1_cons,NB2_cons,Beh1_cons,...
        Beh2_cons,Beh_cons,avg_w1_cons,avg_w2_cons,R1_risk,R2_risk,TNB_risk,NB1_risk,...
        NB2_risk,Beh1_risk,Beh2_risk,Beh_risk,avg_w1_risk,avg_w2_risk,w1_nrsh_cord,...
        w2_nrsh_cord,w1_nrsh_cons,w2_nrsh_cons,w1_nrsh_risk,w2_nrsh_risk,V1_nrsh_cord,...
        V2_nrsh_cord,V1_nrsh_cons,V2_nrsh_cons,V1_nrsh_risk,V2_nrsh_risk]...
        =maincode(dt,nt,t,tmax,R1_vector,R2_vector,BPV1,BPV2,beta1,beta2,...
        disc,w_init,s,m,x_lot,alpha1,alpha2,theta_eq,gamma,gamma_g1,gamma_g2,Ka,Kc,D,...
        phi,c,xN1,xN2,analysis_avg_yrs);
        %% Storage Data Ouputs (See Data Output section for descriptions)
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
    
    %% End Data Processing
    % Fill in the storage matrices with vector values
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
time_elapsed=toc; %elapsed runtime for external script for computational accounting
if save_run == true
    save(FullFileName); %saves output data (from this external script only!) in your root directory
end
toc