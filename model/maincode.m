function [R1_cord,R2_cord,TNB_cord,NB1_cord,NB2_cord,Beh1_cord,Beh2_cord,Beh_cord,...
    avg_w1_cord,avg_w2_cord,R1_cons,R2_cons,TNB_cons,NB1_cons,NB2_cons,Beh1_cons,...
    Beh2_cons,Beh_cons,avg_w1_cons,avg_w2_cons,R1_risk,R2_risk,TNB_risk,NB1_risk,...
    NB2_risk,Beh1_risk,Beh2_risk,Beh_risk,avg_w1_risk,avg_w2_risk,w1_nrsh_cord,...
    w2_nrsh_cord,w1_nrsh_cons,w2_nrsh_cons,w1_nrsh_risk,w2_nrsh_risk,V1_nrsh_cord,...
    V2_nrsh_cord,V1_nrsh_cons,V2_nrsh_cons,V1_nrsh_risk,V2_nrsh_risk]...
    =maincode(dt,nt,t,tmax,R1_vector,R2_vector,BPV1,BPV2,beta1,beta2,...
    disc,w_init,s,m,x_lot,alpha1,alpha2,theta_eq,gamma,gamma_g4,gamma_g3,Ka,Kc,D,...
    phi,c,xN1,xN2,analysis_avg_yrs) % Function that runs this maincode - run from parametric analyses - see documentation

PV1 = BPV1; PV2 = BPV2; % Sets First iteration of PV to Baseline Property Value

%% Benefits/Behaviors Storage
nR1=length(R1_vector); nR2=length(R2_vector); %Gets length of rotation vectors
NB1_matrix=NaN(nR1,nR2); NB2_matrix=NaN(nR1,nR2); TNB_matrix=NaN(nR1,nR2); % Behavior storage matrices
Beh1_matrix=NaN(nR1,nR2); Beh2_matrix=NaN(nR1,nR2); Beh_matrix=NaN(nR1,nR2); % Behavior storage matrices
avg_w1_matrix=NaN(nR1,nR2); avg_w2_matrix=NaN(nR1,nR2); % End width storage matrices
Vadded1_matrix=NaN(nR1,nR2); Vadded2_matrix=NaN(nR1,nR2); % Volume Added Storage Containers
w_added1_matrix=NaN(nR1,nR2); w_added2_matrix=NaN(nR1,nR2); % Width Added storage matrices

%% Model Runs
for iR1=1:nR1 % Loop for Community 1, stored in rows
    NB1_vector=NaN(1,nR1); NB2_vector=NaN(1,nR2); TNB_vector=NaN(1,nR1); %Vectors to store Community 2 outcomes against each Community 1 R-value
    Beh1_vector=NaN(1,nR1); Beh2_vector=NaN(1,nR2); Beh_vector=NaN(1,nR1); % Behaviors storage
    avg_w1_vector=NaN(1,nR1); avg_w2_vector=NaN(1,nR2); % Width vectors
    Vadded1_vector=NaN(1,nR2); Vadded2_vector=NaN(1,nR2); % Volume Added Vectors
    w_added1_vector=NaN(1,nR2); w_added2_vector=NaN(1,nR2);

    for iR2=1:nR2 % Community 2 loop, stored in columns
        %% Model Run Community & Ghost Data:
        R1=R1_vector(iR1); % Selects the rotation length for this run (community 1)
        R2=R2_vector(iR2); % Selects the rotation length for this run (community 2)
        %% Storage Containers for each model run
        theta=zeros(nt,m); qL=zeros(nt,m); qC=zeros(nt,m); fS=zeros(nt,m); fT=zeros(nt,m); 
        xS=zeros(nt,m); xT=zeros(nt,m); w=zeros(nt,m);
        Ben=zeros(nt,m); Ben_disc=zeros(nt,m); TBen=zeros(nt,m);
        C=zeros(nt,m); TC=zeros(nt,m);
        nb=zeros(nt,m);        
        NB1=NaN(1); NB2=NaN(1); TNB=NaN(1);
        PropVal=zeros(nt,m); % Property Value at time of model run (unused in current model)
                
        %% Initial Shoreface Conditions
        xS(1,:)=x_lot(:)+w_init(:); % Sets the initial Shoreline
        xT(1,:)=xS(1,:)+(D/theta_eq); % Set initial shore toe location
        theta(1,:)=D./(xT(1,:)-xS(1,:)); % Set the intial shoreface slopes
        w(1,:) = xS(1,:)-x_lot(1,:); % Set initial Width
    
        %% Initial Benefits / Costs
        % PROPERTY VALUE
        PropVal(1,1)=BPV1*((w(1,1))/w(1,1))^beta1; % Property Value Community 1
        PropVal(1,2)=BPV2*((w(1,2))/w(1,2))^beta2; % Property Value Community 2
        % BENEFITS 
        Ben(1,1)=alpha1*disc*((w(1,1))/w(1,1))^beta1; % Benefits Community 1
        Ben(1,2)=alpha2*disc*((w(1,2))/w(1,2))^beta2; % Benefits Community 2
        Ben_disc(1,:)=Ben(1,:); % Benefits with discounting
        TBen(1,:)=Ben(1,:)*dt; % Total Benefits
        nb(1,:)=TBen(1,:); % Net benefits
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
            qL(i,1)=2*Ka*((xS(i,4)-xS(i,1))/(s(4)+s(1))); % See diagram
            qL(i,2)=2*Ka*((xS(i,1)-xS(i,2))/(s(1)+s(2))); % See diagram
            qL(i,3)=2*Ka*((xS(i,2)-xS(i,3))/(s(2)+s(3))); % See diagram
            qL(i,4)=0; % No sediment flux between the boundary cells
            %s Shoreface Slopes
            theta(i,:)=D./(xT(i,:)-xS(i,:)); % Calculate the shoreface slope for all cells
            % Cross-shore fluxes
            qC(i,:)=Kc.*(theta(i,:)-theta_eq); % Caculcate the cross-shore fluxes for all cells
            
            %% ODE's 
            % Shoreface Toe
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
            % Boundary Cell Toe 
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
            
            %% Benefit      
            %Calculate benefits for community 1
            if w(i+1,1) < 0 % if width is less than zero - property loss occurs
                Ben(i+1,1)=0;
                PropVal(i+1,1)=0;
            else % No property loss, effect of beach width on PV is calculated
                Ben(i+1,1)=alpha1*disc*((w(i+1,1))/w(1,1))^beta1;
                PropVal(i+1,1)=PV1*((w(i+1,1))/w(1,1))^beta1;
            end
            % Calculate benefits for community 2
            if w(i+1,2) < 0 % if width is less than zero - property loss occurs
                Ben(i+1,2)=0;
                PropVal(i+1,2)=0;
            else % No property loss, effect of beach width on PV is calculated
                Ben(i+1,2)=alpha2*disc*((w(i+1,2))/w(1,2))^beta2;
                PropVal(i+1,2)=PV2*((w(i+1,2))/w(1,2))^beta2;
            end

            %% Benefits with Discounting / Total Benefits
            % Discounting of future benefits
            Ben_disc(i+1,1)=Ben(i+1,1)*exp(-disc*t(i+1));
            Ben_disc(i+1,2)=Ben(i+1,2)*exp(-disc*t(i+1)); 
            % Summation of benefits over time
            TBen(i+1,1)=dt*(sum(Ben_disc(:,1)));
            TBen(i+1,2)=dt*(sum(Ben_disc(:,2)));
    
            %% Cost 
            % if statement determines if nourishment & cost occurs 
            if nourish1==0 % Community 1 does not nourish - no costs
                C(i+1,1)=0;
            else % Community 1 nourishes, calculate costs
                C(i+1,1)=(c+nVol1*phi)*exp(-disc*t(i+1));
            end
            if nourish2==0 % Community 2 does not nourish - no costs
                C(i+1,2)=0;
            else % Community 2 nourishes, calculate costs
                C(i+1,2)=(c+nVol2*phi)*exp(-disc*t(i+1));
            end
            TC(i+1,1)=C(i+1,1)+TC(i,1);
            TC(i+1,2)=C(i+1,2)+TC(i,2);
    
            %% Calculate Marginal Net Benefits 
            nb(i+1,1)=TBen(i+1,1)-TC(i+1,1);
            nb(i+1,2)=TBen(i+1,2)-TC(i+1,2);
    
            %% Calculate Net Benefits Per Community and Total Net Benefits - used for decision matrices
            NB1=nb(i+1,1);
            NB2=nb(i+1,2);
            TNB=NB1+NB2;
            
        end
        
        %% Behaviors
        % Min and Max of shoreline location
        analysis_yr_start = (tmax-analysis_avg_yrs)/dt; % find index number for the start of the behavior period 
        min_wS1_ss = min(w(analysis_yr_start:end,1)); % Minimum beach width for community 1
        max_wS1_ss = max(w(analysis_yr_start:end,1)); % Max beach width for community 1
        min_wS2_ss = min(w(analysis_yr_start:end,2)); % Minimum beach width for community 2
        max_wS2_ss = max(w(analysis_yr_start:end,2)); % Max beach width for community 2
        
        % Individual Community Behaviors
        %Community 1
        if isnan(R1) && min_wS1_ss <= 0
            behav_com1 = 9; % full retreat       
        elseif max_wS1_ss > w_init(1) + xN1
            behav_com1 = 0; %seaward growth
        elseif max_wS1_ss <= w_init(1) + xN1 && min_wS1_ss > 0
            behav_com1 = 3; % hold the line
        elseif ~isnan(R1) && min_wS1_ss < 0
            behav_com1 = 6; %slow retreat
        end
        %Community 2
        if isnan(R2) && min_wS2_ss <= 0
            behav_com2 = 9; % full retreat       
        elseif max_wS2_ss > w_init(2) + xN2
            behav_com2 = 0; %seaward growth
        elseif max_wS2_ss <= w_init(2) + xN2 && min_wS2_ss > 0
            behav_com2 = 3; % hold the line
        elseif ~isnan(R2) && min_wS2_ss < 0
            behav_com2 = 6; %slow retreat
        end
        %Combined Behavior
        if behav_com1 + behav_com2 == 0
            behavior_ss = 0; %seaward growth
        elseif behav_com1 + behav_com2 == 3
            behavior_ss = 1.5; % seaward growth / hold the line
        elseif behav_com1 == 3 && behav_com2 == 3 
            behavior_ss = 3; % hold the line
        elseif (behav_com1 == 0 && behav_com2 == 6) || (behav_com1 == 6 && behav_com2 == 0)
            behavior_ss = 4.5; %mixed: seaward growth / slow retreat
        elseif (behav_com1 == 3 && behav_com2 == 6) || (behav_com1 == 6 && behav_com2 == 3)
            behavior_ss = 4.5; %mixed: hold the line / slow retreat
        elseif (behav_com1 == 0 && behav_com2 == 9) || (behav_com1 == 9 && behav_com2 == 0)
            behavior_ss = 4.5; %seaward growth / retreat
        elseif (behav_com1 == 3 && behav_com2 == 9) || (behav_com1 == 9 && behav_com2 == 3)
            behavior_ss = 4.5; % mixed hold the line / retreat
        elseif (behav_com1 == 6 && behav_com2 == 6) 
            behavior_ss = 6; % slow retreat
        elseif behav_com1 + behav_com2 == 15
            behavior_ss = 7.5; %slow retreat / retreat
        elseif behav_com1 + behav_com2 == 18
            behavior_ss = 9; % Full retreat
        else
            behavior_ss = NaN; % not identified
        end
          
        %% Calculate Nourishment Volumes and Widths Added
        % Calculate the average width over observation time period (delta)
        avg_width1=mean(w(analysis_yr_start:end,1));
        avg_width2=mean(w(analysis_yr_start:end,2));

        % Calculate total volume and width added if community nourishes
        if isnan(R1) % If NaN, community 1 did not nourish at all
            w_added1=0;   
            V_added1=0;
        else % Community 1 nourished, calculate volumes and total width    
            w_added1=(k1-1)*xN1;  % total beach width added over course of model run
            V_added1=(k1-1)*0.5*xN1*s(1)*D;  % Total Volume added per community
        end
        if isnan(R2) % If NaN, community 2 did not nourish at all
            w_added2=0;
            V_added2=0;
        else % Community 2 nourished, calculate volumes and total width
            w_added2=(k2-1)*xN2; % total beach width added over course of model run
            V_added2=(k2-1)*0.5*xN2*s(2)*D; % Total Volume added per community
        end

        %% Data for Each Community 2 Rotation - Community 2 is stored in columns
        NB1_vector(iR2)=NB1; NB2_vector(iR2)=NB2; TNB_vector(iR2)=TNB; % Benefits storage inner loop vector
        Beh1_vector(iR2)=behav_com1; Beh2_vector(iR2)=behav_com2; Beh_vector(iR2)=behavior_ss; % Behaviors inner loop vector
        avg_w1_vector(iR2)=avg_width1; avg_w2_vector(iR2)=avg_width2; % Width inner loop vector
        Vadded1_vector(iR2)=V_added1; Vadded2_vector(iR2)=V_added2; % Volume inner loop vector
        w_added1_vector(iR2)=w_added1; w_added2_vector(iR2)=w_added2; % Width Added inner loop vector
    end
    %% Data for Each Community 1 Rotation - Community 1 is stored in rows
    NB1_matrix(iR1,:)=NB1_vector(1,:); NB2_matrix(iR1,:)=NB2_vector(1,:); TNB_matrix(iR1,:)=TNB_vector(1,:); %Benefits
    Beh1_matrix(iR1,:)=Beh1_vector(1,:); Beh2_matrix(iR1,:)=Beh2_vector(1,:); Beh_matrix(iR1,:)=Beh_vector(1,:); % Behaviors
    avg_w1_matrix(iR1,:)=avg_w1_vector(1,:); avg_w2_matrix(iR1,:)=avg_w2_vector(1,:); %Width storage matrix
    Vadded1_matrix(iR1,:)=Vadded1_vector(1,:); Vadded2_matrix(iR1,:)=Vadded2_vector(1,:); % Volume storage matrix
    w_added1_matrix(iR1,:)=w_added1_vector(1,:); w_added2_matrix(iR1,:)=w_added2_vector(1,:); % Stores width added for all model runs
end

%% Coordination - Economic - Optimizations
% Coordination find maximums
max_coord=max(max(TNB_matrix(:))); %max net benefit in matrix
[row_coord,col_coord]=find(TNB_matrix==max_coord); %finds location in matrix where maximum net benefit occurs
row_c=row_coord(1); col_c=col_coord(1); %ensures only one combination is output to external script in case optimization doesn't converge due to resolution issue (only occurs under specific parameter combinations and resolution issue is confirmed as the cause)

% Coordination Data Outputs
R1_cord=R1_vector(row_c); %optimal rotation length in community 1
R2_cord=R2_vector(col_c); %optimal rotation length in community 2
NB1_cord=NB1_matrix(row_c,col_c); %net benefit realized by community 1
NB2_cord=NB2_matrix(row_c,col_c); %net benefit realized by community 2
TNB_cord=TNB_matrix(row_c,col_c); %total net benefits under coordination
Beh1_cord=Beh1_matrix(row_c,col_c); % Community 1 Behavior under coordination
Beh2_cord=Beh2_matrix(row_c,col_c); % Community 2 Behavior under coordination
Beh_cord=Beh_matrix(row_c,col_c); %emergent behavior for combination of optimal rotation lengths under coordination
avg_w1_cord=avg_w1_matrix(row_c,col_c); % Width 1 under coordination
avg_w2_cord=avg_w2_matrix(row_c,col_c); % Width 2 under coordination
V1_nrsh_cord=Vadded1_matrix(row_c,col_c); % Volume Added to Comm 1 Coordination
V2_nrsh_cord=Vadded2_matrix(row_c,col_c); % Volume Added to Comm 2 Coordination
w1_nrsh_cord=w_added1_matrix(row_c,col_c); % Width added to Community 1 Coordination
w2_nrsh_cord=w_added2_matrix(row_c,col_c); % Width added to Community 2 Coordination

%% Non-Coordination: Conservative
% Find Maximums per community
% Community 1
Max1_cons=max(NB1_matrix(:,1)); % maximum community 1 net benefit assuming neighbor does nothing (conservative assumption)
[row1_c1_nc_consf,col1_c1_nc_consf]=find(NB1_matrix==Max1_cons); % location in net benefit matrix
row1_c1_nc_cons=row1_c1_nc_consf(1); col1_c1_nc_cons=col1_c1_nc_consf(1); %ensures only one selection occurs
% Communnity 2
Max2_nc_cons=max(NB2_matrix(1,:)); % maximum community 1 net benefit assuming neighbor does nothing (conservative assumption)
[row1_c2_nc_consf,col1_c2_nc_consf]=find(NB2_matrix==Max2_nc_cons); % location in net benefit matrix
row1_c2_nc_cons=row1_c2_nc_consf(1); col1_c2_nc_cons=col1_c2_nc_consf(1); %ensures only one selection occurs

% Data outputs
R1_cons=R1_vector(row1_c1_nc_cons); %optimal rotation length in community 1 under conservative non-coordination
R2_cons=R2_vector(col1_c2_nc_cons); %optimal rotation length for com.2 under conservative non-coordination
TNB_cons=TNB_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Finds total net benefits based on individual community decisions
NB1_cons=NB1_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Net benefits for community 1 under conservative non-coordination
NB2_cons=NB2_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Net benefits for community 2 under conservative non-coordination
Beh1_cons=Beh1_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Community 1 Behavior under conservative non-coordination
Beh2_cons=Beh2_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Community 2 Behavior under conservative non-coordination
Beh_cons=Beh_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Behavior for both communities under conservative non-coordination
avg_w1_cons=avg_w1_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Width 1 under conservative non-coordination
avg_w2_cons=avg_w2_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Width 2 under conservative non-coordination
V1_nrsh_cons=Vadded1_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Volume Added to Comm 1 conservative non-coordination
V2_nrsh_cons=Vadded2_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Volume Added to Comm 2 conservative non-coordination
w1_nrsh_cons=w_added1_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Width added to Community 1 Conservative Non-Coordination
w2_nrsh_cons=w_added2_matrix(row1_c1_nc_cons,col1_c2_nc_cons); % Width added to Community 2 Conservative Non-Coordination

%% Non-Coordination: Risky
% Find Maximums per community
% Community 1
Max1_nc_risk=max(NB1_matrix(:)); % maximum community 1 net benefit assuming neighbor does nothing (conservative assumption)
[row1_c1_nc_riskf,col1_c1_nc_riskf]=find(NB1_matrix==Max1_nc_risk); % location in net benefit matrix
row1_c1_nc_risk=row1_c1_nc_riskf(1); col1_c1_nc_risk=col1_c1_nc_riskf(1); %ensures only one selection occurs
% Communnity 2
Max2_nc_risk=max(NB2_matrix(:)); % maximum community 1 net benefit assuming neighbor does nothing (conservative assumption)
[row1_c2_nc_riskf,col1_c2_nc_riskf]=find(NB2_matrix==Max2_nc_risk); % location in net benefit matrix
row1_c2_nc_risk=row1_c2_nc_riskf(1); col1_c2_nc_risk=col1_c2_nc_riskf(1); %ensures only one selection occurs

% Data outputs
R1_risk=R1_vector(row1_c1_nc_risk); %optimal rotation length in community 1 under risky non-coordination
R2_risk=R2_vector(col1_c2_nc_risk); %optimal rotation length for com.2 under risky non-coordination
TNB_risk=TNB_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Finds total net benefits based on individual community decisions
NB1_risk=NB1_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Net benefits for community 1 under risky non-coordination
NB2_risk=NB1_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Net benefits for community 2 under risky non-coordination
Beh1_risk=Beh1_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Community 1 Behavior under risky non-coordination
Beh2_risk=Beh2_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Community 2 Behavior under risky non-coordination
Beh_risk=Beh_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Behavior for both communities under risky non-coordination
avg_w1_risk=avg_w1_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Width 1 under risky non-coordination
avg_w2_risk=avg_w2_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Width 2 under risky non-coordination
V1_nrsh_risk=Vadded1_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Volume Added to Comm 1 risky non-coordination
V2_nrsh_risk=Vadded2_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Volume Added to Comm 2 risky non-coordination
w1_nrsh_risk=w_added1_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Width added to Community 1 Risky Non-Coordination
w2_nrsh_risk=w_added2_matrix(row1_c1_nc_risk,col1_c2_nc_risk); % Width added to Community 2 Risky Non-Coordination