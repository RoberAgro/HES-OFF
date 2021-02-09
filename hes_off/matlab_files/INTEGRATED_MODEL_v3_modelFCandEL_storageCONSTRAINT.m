
global e workbook file eSheets eSheet1 GT_type

%Call "Excel_GT.m" file first
%That file while activate the excel file with the GT model
%The type of GT to be used is chosen there as well 
%GT_type = 1 --> LM2500
%GT_type = 2 --> LM6000
%For a faster run the type of model can be directly chosen here
%In this case, the calculations will be based on simplified correlations

%%
%Script for wind power

%Wind turbine to use (power curve)
Power_curve = Hywind;
WT_type = 2;
% Power_curve = NREL;
% WT_type = 1;

% %Data series of wind speed instances throughout one year (one minute resolution)
% %SLEIPNER
wdata = load('sleipnerwind.mat');
wspeed = wdata.wind(1).w;
wdata_source = 1;
% %SPYROS 2018
% wdata = load('ninja_wind_2018.mat');
% wspeed = wdata.ninja_wind_2018(1:8760);
% wdata_source = 2;
% %SPYROS 2017
% wdata = load('ninja_wind_2017.mat');
% wspeed = wdata.ninja_wind_2017(1:8760);
% wdata_source = 3;
% %SPYROS 2016
% wdata = load('ninja_wind_2016.mat');
% wspeed = wdata.ninja_wind_2016(1:8760);
% wdata_source = 4;

%Multiplication factor for wind speed (to adjust for bigger hub heigth etc)
wspeed_factor = (90/14)^(1/7); %1.25; 
  %90 = wind turbine hub height (90m); 
  %14 = reference height of wind speed data (it is likely 14 m); 
  %The exponent 1/7 (0.144) is commonly used to extrapolate the wind speed taken at a reference height. 
  %More info in whartonWindEnergy2010.pdf (pp. 15-16) the file is on the Sharepoint
 
%Normalized power. One instance per minute
PW_norm = windpower(Power_curve,wspeed,wspeed_factor);

%Considering the hourly average to decrease computational effort. One instance per hour
PW_norm_h = zeros(8760,1);
if wdata_source == 1
    for i=1:8760
       for j=1:60       
        PW_norm_h(i,1) = PW_norm_h(i,1)+PW_norm(j+60*(i-1));    
       end
       PW_norm_h(i,1)=PW_norm_h(i,1)/60;
    end
else
    for i=1:8760
        PW_norm_h(i,1)=PW_norm(i,1);
    end    
end



%%
%List of assumptions

%POWER DEMAND
%Three cases considered: peak demand, mid-life demand and tail demand
%Power demand (MW)
POWER = [43.6, 35.2, 32.9]; %[45, 35, 30]; %

%WIND FARM
%Wind capacity installed (MW)
WC = 36;    %-----> Decision variable <-----%

%GAS TURBINE
%GT type
% %If not using the excel model of GT then we need to specify GT type
GT_type = 1; %LM2500
% GT_type = 2; %LM6000
%Design GT power (MW) (power provided by the GT at design - 100% load)
if GT_type == 1
    PW_GT_d = 33.31; %LM2500
else
    PW_GT_d = 44.70; %LM6000
end
%Possibility to start a second GT %1 yes // 2 no
if GT_type == 1
    GT2=1;  %1 yes // 2 no
    if GT2==1
        GT2_only_peak=1;  %Is second GT avaialable only at peak years? 1 yes // 2 no
    else
        GT2_only_peak=2;  %Is second GT avaialable only at peak years? 1 yes // 2 no
    end   
else
    GT2=2; %1 yes // 2 no
    if GT2==1
        GT2_only_peak=1;  %Is second GT avaialable only at peak years? 1 yes // 2 no
    else
        GT2_only_peak=2;  %Is second GT avaialable only at peak years? 1 yes // 2 no
    end    
end
%Variable that keeps tracks of #GT operating at any instance
Num_GT = ones(8760,3);
%Max GT load
max_GT_l = 0.95; %1.00;     %-----> Decision variable <-----%
%Min GT load (GT cannot be operated at a lower load)
min_GT_l = 0.40;     %-----> Decision variable <-----%
%Max fraction of H2 (volumetric) that can be used in GT  %-----> Decision variable <-----%
if GT_type == 1
    max_H2_f = 0.032; %0.015 kg/s is about 10%vol H2 inlet at GT load 40 -- 0.018 kg/s about 12% -- 0.032 kg/s about 20%
else
    max_H2_f = 0.037; %0.017 kg/s is about 10%vol H2 inlet at GT load 40 -- 0.02 kg/s about 12% -- 0.037 kg/s about 20%
end

%ELECTROLYSER
%Electrolyser performance (kgH2/MJ) - if constant efficiency is used
perf_el = 0.8*1/120; %0.005; %
%Max size electrolsyer
max_size_el_peak = 10; %MW  %-----> Decision variable <-----%
%Size elctrolyser
size_electrolyser = 0; %MW %Initialised %Check the size of the EL actually used
%Size single EL in a stack
P_ELY_module = 0.500;   %MW; 
%Coefficients for H2_produced (Marcin's model)
C3=-1.41175e-20;
C2 = 1.59861e-14;
C1 = - 8.29304e-09; 
C0 = 6.82299e-03; 
%Max number of module in stacks
max_modules_el = max_size_el_peak/P_ELY_module;
%Max H2 produced based on max size electrolsyer
no_EL_on = max_size_el_peak/P_ELY_module; %Number of cells of the EL stack to activate
no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
no_EL_on = ceil(no_EL_on); %Round to the next integer number
if no_EL_on == 0
    PW_surplus_cell = 0; %[W]
else
    PW_surplus_cell = max_size_el_peak/no_EL_on*10^6; %[W]
end
H2_max_size_el_peak = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
H2_max_size_el_peak = H2_max_size_el_peak*no_EL_on;

%FUEL CELL
%Fuel cell performance (MJ/kgH2) - if constant efficiency is used
perf_fc = 66; %53.4; %0.5*120;
%Size fuel cell (initialized)
size_FC = 3.56;
%Size single FC of a stack
P_FCS_module = 0.125; %MW
%Constant -> (Total active area in cm2/(2*Faraday's constant))*(mol/s->kg/h)
A=(313800/(2*96485))*(2/1000); %Total active area based on 2 MW FCS unit made by NEDSTACK
%Coefficients for H2_usage
C2_FC = 7.767e-12;
C1_FC = 4.47e-06; 
C0_FC = 0.003237;  

%STORAGE
%Max storage (kg)
limit_storage = 5000; %-----> Decision variable <-----%
max_storage_surplus = limit_storage/2; % 5000;  %-----> Decision variable <-----%
max_storage_lack = limit_storage/2; %5000;  %-----> Decision variable <-----%
%Size energy storage
size_storage = 0; %kg %Initialised

%HYDROGEN IN GAS TURBINE
%Auxiliary variable to check if we used H2 in GT
H2_in_GT = zeros(8760,3);  
H2_GT = zeros(1,3);
%Lower heating values
LHV_H2 = 119.960; % Mj/kg
LHV_NG = 50.047; % Mj/kg  %CH4=50.047 / NG=47.383
%Molecular weights
MW_H2 = 2.01588; %kg/kmol
MW_NG = 16.0425; %kg/kmol  %CH4=16.0425 / NG=22.0998

%AUXILIARY VALUES TO KEEP TRACK OF:
%wind power
PW_wind_h = zeros(8760,1);
%wind power dissipated
PW_wind_h_dissipated = zeros(8760,3);
%Variable that keeps tracks of actual GT loads
GT_load_act = zeros(8760,3);
%Variable that keeps tracks of actual GT power
PW_GT_load_act = zeros(8760,3);
%Variable that keeps tracks of actual GT emissions
CO2_act = zeros(8760,3);
%Variable that keeps tracks of energy input
en_in = zeros(8760,3);
%the hydrogen produced by electrolyser
H2_produced_h = zeros(8760,3); 
%the hydrogen used in FC
H2_used_fc_h = zeros(8760,3);    
%the storage level
H2_stored_h = zeros(8760,3);
H2_stored_h_aux = zeros(8760,1);
%power used from fuel cell
PW_fc_h = zeros(8760,3);
%extra power to harvest from fuel cell
extra_FC= zeros(8760,3);

%CO2 EMISSIONS
%Cumulative CO2 emissions
CO2_tot = 0; %Initialised

%%
%ITERATIONS throghout the year for the three power demands
for j=1:3 %Three cases considered: peak demand, mid-life demand and tail demand
  
    %Define some auxiliary variables
    H2_stored = 0; %9999999;       
    H2_stored_aux = 0;
    H2_stored_aux2 = 0;
    H2_stored_aux2_aux=0;
    H2_stored_aux_adj = 0;
    
    min_level_storage = 0;
    max_level_storage = 0;    
    
    %FIRST ITERATION throughout the year - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %In this iteration we evaluate the H2 needed (because of power deficits) or the H2 produced (because of power surplus)
    %It is not evaluated how to produce or consume H2    
    for i=1:8760
        
        %Contribution of wind power at the hour i
        PW_wind_h(i,1)=PW_norm_h(i,1)*WC;
        %Remaining power to meet total demand at the hour i
        PW_remaining = POWER(1,j) - PW_wind_h(i,1);
        
        %Assumption that all remaining power should be covered by GT
        PW_GT_1 = PW_remaining;
        GT_load_1 = PW_GT_1/PW_GT_d; %related GT load to generate such power
        
        if GT_load_1 > max_GT_l  %A single GT is not able to cover the power demand
                                 %We need contribution of fuel cell or start a second GT
                       
            if GT2 == 1 & j == 1 %Only at peak conditions and with the LM25000 we can start a second GT
                GT_load_2 = PW_GT_1/PW_GT_d/2; %the two GTs shares equally the load
                PW_GT_2 = PW_GT_d*GT_load_2;
                H2_stored_h(i,j) = H2_stored;   %no H2 is produced or consumed. The storage level is constant                             
                GT_load_act(i,j) = GT_load_2;
                PW_GT_load_act(i,j) = PW_GT_2*2;
                Num_GT(i,j) = 2; 
            else %We do not activate a second GT. We need fuel cell to meet power demand.
                GT_load_2 = max_GT_l; %GT produces as much power as possible
                PW_GT_2 = PW_GT_d*GT_load_2;          
                PW_deficit = PW_remaining - PW_GT_2; %Deficit of power
                
                %If we use polynomial interpolation from Marcin's model
                no_FC_on = PW_deficit/P_FCS_module; %Number of cells of the FC stack to activate
                no_FC_on = no_FC_on - 0.001;  %To avoid wrong rounding to next integer
                no_FC_on = ceil(no_FC_on); %Round to the next integer number
                if no_FC_on >0
                    PW_deficit_cell = PW_deficit/no_FC_on*10^6; %[W]  
                else
                    PW_deficit_cell = 0;
                end     
                H2_used_fc = A*(C2_FC*PW_deficit_cell^2+C1_FC*PW_deficit_cell+C0_FC); %kg/s
                H2_used_fc = H2_used_fc*no_FC_on;
                H2_used_fc_h(i,j) = H2_used_fc;
                
%                 %If we use constant FC efficiency
%                 H2_used_fc = PW_deficit/perf_fc; %H2 used to produce power with fuel cell
                
                H2_stored = H2_stored - H2_used_fc*3600; %H2 stored (kg) --> H2 is consumed in FC and the storage level decreases
                H2_stored_h(i,j) = H2_stored;
                PW_fc_h(i,j) = PW_deficit;
                GT_load_act(i,j) = GT_load_2;
                PW_GT_load_act(i,j) = PW_GT_2;                
            end
            
        else if GT_load_1 < min_GT_l  %This is the case where we have excess of power from wind. We need to produce some H2.
                
                GT_load_2 = min_GT_l; %GT produces the minimum possible amount of power
                PW_GT_2 = PW_GT_d*GT_load_2;
                PW_surplus = PW_GT_2-PW_remaining; %Surplus power to be converted in H2
                
                %If we use polynomial interpolation from Marcin's model
                no_EL_on = PW_surplus/P_ELY_module; %Number of cells of the EL stack to activate
                no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                no_EL_on = ceil(no_EL_on); %Round to the next integer number
                if no_EL_on == 0
                    PW_surplus_cell = 0; %[W]
                else
                    PW_surplus_cell = PW_surplus/no_EL_on*10^6; %[W]
                end
                H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                H2_produced = H2_produced*no_EL_on;                
                
%                 %If we use constant EL efficiency
%                 H2_produced = PW_surplus*perf_el; %H2 produced from surplus power in electrolyser (kg/s)       
                
                
                if H2_produced > H2_max_size_el_peak  %We would need an electrolyse stack larger than what designed
                    H2_produced_old = H2_produced;
                    
                    %If we use polynomial interpolation from Marcin's model
                    no_EL_on = max_size_el_peak/P_ELY_module; %Number of cells of the EL stack to activate
                    no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                    no_EL_on = ceil(no_EL_on); %Round to the next integer number
                    if no_EL_on == 0
                        PW_surplus_cell = 0; %[W]
                    else
                        PW_surplus_cell = max_size_el_peak/no_EL_on*10^6; %[W]
                    end
                    H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                    H2_produced = H2_produced*no_EL_on;
                    PW_wind_h_dissipated(i,j) = -(POWER(1,j) - PW_wind_h(i,1) - PW_GT_2 + max_size_el_peak);                    
                
%                     %If we use constant EL efficiency
%                     H2_produced = max_size_el_peak*perf_el;                    
%                     PW_wind_h_dissipated(i,j) = (H2_produced_old - H2_produced)/perf_el;
                
                end

                H2_stored = H2_stored + H2_produced*3600- H2_in_GT(i,j)*3600; %H2 stored (kg) --> H2 is produced in EL and the storage level increases
                H2_stored_h(i,j) = H2_stored;
                H2_produced_h(i,j) = H2_produced;
                GT_load_act(i,j) = GT_load_2;
                PW_GT_load_act(i,j) = PW_GT_2;
                                
            else  %This is the case where we can produce power simply with GT and wind
                
                H2_stored_h(i,j) = H2_stored;   %no H2 is produced or consumed. The storage level is constant             
                GT_load_act(i,j) = GT_load_1;
                PW_GT_load_act(i,j) = PW_GT_1;
                
            end
                        
        end
                   
    end    
    %End of first run throughout the year %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check the size of the FC which will be defined in peak years where the power deficit is maximum
    size_FC_aux = max(PW_fc_h(:,j)); %Size of FC in MW
    if size_FC_aux > size_FC
       size_FC = size_FC_aux; %+2;
    end
    
    %SECOND ITERATION throughout the year - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %In this iteration we evaluate the primary strategies to produce H2 or consume H2
    %First attempt to set the storage balances to zero
    H2_stored
    %Auxiliary variables for H2 storage
    H2_stored_aux = H2_stored;  %H2_stored_aux keeps tracks of the H2 still to produce or consume to reach zero balance
                                %The starting point is the net balance throughout a year calculated in first iteration (H2_stored)
    %H2_stored_h(:,j) = H2_stored_h_aux(:,1);
    H2_stored_aux2 = 0;  %H2_stored_aux2 keeps tracks of the H2 storage level with the strategy used in this second iteration
                         %The starting point is zero
    if H2_stored < 0  %During the first iteration we registered a net negative balance of H2
                      %Lack of power (typical of peak years). Need to produce the H2 to be used in the fuel cell
                      %The primary approach is to increase GT load when possible and to use extra power in electrolyser to produce H2
     
     %Define storage factor (SF)
     %Check if I can avoid having this!!!!!!!!!!!!!!!!!!!!!!!!!!
     if H2_stored < (-10^5)
         SF=10;
     else
         SF=1;
     end
     
        for i=1:8760
             
             if H2_stored_aux < 0  %This to keep track that we haven't yet produced all the H2 needed

                if GT_load_act(i,j) == max_GT_l  %This are the cases when we do not have enough power and we need to use fuel cell
                                                 %Here we cannot produce H2
                    GT_load_act(i,j) = GT_load_act(i,j);
                    
                    %If we use polynomial interpolation from Marcin's model
                    no_FC_on = PW_fc_h(i,j)/P_FCS_module; %Number of cells of the FC stack to activate
                    no_FC_on = no_FC_on - 0.001;  %To avoid wrong rounding to next integer
                    no_FC_on = ceil(no_FC_on); %Round to the next integer number
                    if no_FC_on >0
                        PW_cell = PW_fc_h(i,j)/no_FC_on*10^6; %[W]  
                    else
                        PW_cell = 0;
                    end 
                    H2_stored_aux2 = H2_stored_aux2 - A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC)*3600*no_FC_on; %kg
                    
%                     %If we use constant FC efficiency
%                     H2_stored_aux2 = H2_stored_aux2 - PW_fc_h(i,j)/perf_fc*3600;
                    
                    %Tracking max and min level of storage
                    if H2_stored_aux2 > max_level_storage
                        max_level_storage = H2_stored_aux2;
                    else if H2_stored_aux2 < min_level_storage
                            min_level_storage = H2_stored_aux2;
                        end                    
                    end
                                      
                    %If we need to constraint the maximum size of the storage
                    %But this will require to activate a second GT                    
                    if (max_level_storage-H2_stored_aux2) > limit_storage & GT2 == 1 %We reached max dimension allowed to storage system & we allow using 2nd GT
                        Num_GT(i,j) = 2;
                        GT_load_act(i,j) = 0.5*(GT_load_act(i,j)*PW_GT_d+PW_fc_h(i,j))/PW_GT_d;
                        %If we use polynomial interpolation from Marcin's model
                        H2_stored_aux = H2_stored_aux + A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC)*3600*no_FC_on; %We are not using FC anymore thus the requirement for H2 decreases
                        H2_stored_aux2 = H2_stored_aux2 + A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC)*3600*no_FC_on; %We undo the calculation above (we do not any longer consume H2 in FC)
%                         %If we use constant FC efficiency
%                         H2_stored_aux = H2_stored_aux + PW_fc_h(i,j)/perf_fc*3600; %We are not using FC anymore thus the requirement for H2 decreases
%                         H2_stored_aux2 = H2_stored_aux2 + PW_fc_h(i,j)/perf_fc*3600; %We undo the calculation above (we do not any longer consume H2 in FC)
                        H2_stored_h(i,j) = H2_stored_aux2;
                        PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j)*Num_GT(i,j);
                        PW_fc_h(i,j) = 0;
                        H2_used_fc_h(i,j)=0;
                    end
                                        
                else  %GT_load_act(i,j) ~= max_GT_l %Here there is enough power and we can produce some H2 to be later used
                    
                    H2_produced_already=0;
                    PW_EL_already = 0;
                    if GT_load_act(i,j) == min_GT_l %We are already producing H2 not to waste wind
                        H2_produced_already=H2_produced_h(i,j);
                        PW_EL_already = PW_wind_h(i,1)+GT_load_act(i,j)*PW_GT_d-POWER(1,j)-PW_wind_h_dissipated(i,j);
                    end
                    
                    GT_old = GT_load_act(i,j);  %Original GT load                   
                    GT_load_act(i,j) = max_GT_l;  %If we set GT at max load...
                    
                    %If we use polynomial interpolation from Marcin's model
                    no_EL_on = (PW_GT_d*(GT_load_act(i,j)-GT_old)+PW_EL_already)/P_ELY_module; %Number of cells of the EL stack to activate
                    no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                    no_EL_on = ceil(no_EL_on); %Round to the next integer number
                    if no_EL_on==0
                        PW_surplus_cell = 0; %[W]
                    else
                        PW_surplus_cell = PW_GT_d*(GT_load_act(i,j)-GT_old)/no_EL_on*10^6; %[W]
                    end                    
                    H2_produced_bis = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                    H2_produced_bis = H2_produced_bis*no_EL_on; %...the H2 produced would be this %+H2_produced_already
%                     %If we use constant EL efficiency
%                     H2_produced_bis = PW_GT_d*(GT_load_act(i,j)-GT_old)*perf_el+H2_produced_already; %...the H2 produced would be this
%                     %H2_produced_bis = -(POWER(1,j)-PW_wind_h(i,1)-GT_load_act(i,j)*PW_GT_d)*perf_el;
                    
                    H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_bis*3600; %the storage tank would fill until this level
                    
                    %Check we are not producing more H2 than what allowed by electrolyser stack
                    %If we use polynomial interpolation from Marcin's model
                    if no_EL_on > max_modules_el  %we would need bigger electrolyser stack...
                        %GT_load_act(i,j)=(max_size_el_peak+GT_old*PW_GT_d)/PW_GT_d; %...thus we reduced H2 produced                        
                        GT_load_act(i,j)=(POWER(1,j)-PW_wind_h(i,1)+max_size_el_peak)/PW_GT_d; %...thus we reduced H2 produced
                        if GT_load_act(i,j)<min_GT_l
                            PW_wind_h_dissipated(i,j)=PW_GT_d*(min_GT_l-GT_load_act(i,j));                            
                            GT_load_act(i,j)=min_GT_l;
                        end                        
                        %no_EL_on = PW_GT_d*(GT_load_act(i,j)-GT_old)/P_ELY_module; %Number of cells of the EL stack to activate
                        no_EL_on =  max_modules_el; %-(POWER(1,j)-PW_wind_h(i,1)-GT_load_act(i,j)*PW_GT_d+PW_wind_h_dissipated(i,j))/P_ELY_module; %Number of cells of the EL stack to activate
                        no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                        no_EL_on = ceil(no_EL_on); %Round to the next integer number
                        if no_EL_on==0
                            PW_surplus_cell = 0;
                        else
                            PW_surplus_cell = P_ELY_module*10^6; %-(POWER(1,j)-PW_wind_h(i,1)-GT_load_act(i,j)*PW_GT_d+PW_wind_h_dissipated(i,j))/no_EL_on*10^6; %[W]
                        end                    
                        H2_produced_bis = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                        H2_produced_bis = H2_produced_bis*no_EL_on;
                        H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_bis*3600; %the storage tank would fill until this level
                    end                                    
                    
                    %Check we are not overpassing maximum allowed storage size                     
                    if H2_stored_aux2_aux>max_storage_surplus  %we reached max storage level thus we need to decrease the production of H2
                        H2_produced_bis=0;
                        GT_load_act(i,j)=(POWER(1,j)-PW_wind_h(i,1))/PW_GT_d; 
                        delta_PW=0;
                        if GT_load_act(i,j)<min_GT_l %we have to run GT at least at minimum load --> extra power
                            GT_load_act(i,j)=min_GT_l;
                            delta_PW=PW_GT_d*GT_load_act(i,j)+PW_wind_h(i,1)-POWER(1,j); %Extra power whether used in EL or dissipated
                            %Since we cannot store any more H2, we calculate how much H2 from EL we can send to GT (so we do not dissipate it)
                            H2_while = 0;
                            PW_EL_toH2GT = 0;
                            while H2_while < max_H2_f %Check how much we can use the electrolyser to produce H2 that is used in the GT
                                PW_EL_toH2GT = PW_EL_toH2GT+1000; %W
                                no_EL_on = PW_EL_toH2GT/10^6/P_ELY_module;
                                no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                                no_EL_on = ceil(no_EL_on); %Round to the next integer number
                                if no_EL_on==0
                                    PW_surplus_cell = 0;
                                else
                                    PW_surplus_cell = PW_EL_toH2GT/no_EL_on; %[W]
                                end                    
                                H2_while = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                                H2_while = H2_while*no_EL_on;                                
                            end
                           
                            PW_EL_toH2GT = PW_EL_toH2GT - 1000; %W
                            PW_max_H2_GT = PW_EL_toH2GT/10^6; %MW %Max power used in the EL to produce H2 to send to GT (no storage)
                            
                            if PW_max_H2_GT>max_size_el_peak
                                PW_max_H2_GT = max_size_el_peak;
                            end                            
                            if delta_PW<PW_max_H2_GT %All extra power can be used in EL to produce H2 to send to GT (no storage)                 
                                no_EL_on = delta_PW/P_ELY_module;
                                no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                                no_EL_on = ceil(no_EL_on); %Round to the next integer number
                                if no_EL_on==0
                                    PW_surplus_cell = 0;
                                else
                                     PW_surplus_cell = delta_PW/no_EL_on*10^6; %[W]
                                end
                                H2_in_GT(i,j) = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                                H2_in_GT(i,j) = H2_in_GT(i,j)*no_EL_on;
                                H2_produced_bis = H2_in_GT(i,j);
                            else %Some power must be dissipated
                                no_EL_on = PW_max_H2_GT/P_ELY_module;
                                no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                                no_EL_on = ceil(no_EL_on); %Round to the next integer number
                                if no_EL_on==0
                                    PW_surplus_cell = 0;
                                else
                                     PW_surplus_cell = PW_max_H2_GT/no_EL_on*10^6; %[W]
                                end
                                H2_in_GT(i,j) = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                                H2_in_GT(i,j) = H2_in_GT(i,j)*no_EL_on;
                                H2_produced_bis = H2_in_GT(i,j);
                            end
                            if delta_PW<PW_max_H2_GT 
                                PW_wind_h_dissipated(i,j)=PW_GT_d*GT_load_act(i,j)+PW_wind_h(i,1)-POWER(1,j)-delta_PW;
                            else
                                PW_wind_h_dissipated(i,j)=PW_GT_d*GT_load_act(i,j)+PW_wind_h(i,1)-POWER(1,j)-PW_max_H2_GT;
                            end                            
                        end
                        
                    end  %end of --> if (H2_stored_aux2_aux-min_level_storage) > limit_storage

                    H2_produced_h(i,j)=H2_produced_bis; %######
                    H2_stored_aux2 = H2_stored_aux2 + H2_produced_bis*3600-H2_in_GT(i,j)*3600;
                    H2_stored_aux = H2_stored_aux + H2_produced_bis*3600-H2_produced_already*3600-H2_in_GT(i,j)*3600;
                    H2_stored_h(i,j) = H2_stored_aux2;
                    PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j);
                                        
                    %Tracking max and min level of storage
                    if H2_stored_aux2 > max_level_storage
                        max_level_storage = H2_stored_aux2;
                    else if H2_stored_aux2 < min_level_storage
                            min_level_storage = H2_stored_aux2;
                        end                    
                    end

                end %end of --> if GT_load_act(i,j) == max_GT_l
             end %end of --> if H2_stored_aux < 0
            
           end
        
        H2_stored = H2_stored_aux;
        
%         if abs(H2_stored)>max_storage
%            here2 = 'we could not constrain the size of the H2 storage'
%            pause
%         end

    else  %H2_stored > 0
          %During the first iteration we registered a net positive balance of H2
          %Excess of power (typical of tail years). Need to use the H2 produced
          %The primary approach is to use the H2 in the FC stack rather than in the GT (better efficiency)

        for i=1:8760  
            if H2_stored_aux > 0  %This to check that we did not already used all H2 we produced

                if GT_load_act(i,j) == min_GT_l  %These are the cases when we cannot produce less from GT and we need to store some H2
                                                 %Here we can use H2 only in GT (which we do only if we are already over max storage size)
                    
                    H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_h(i,j)*3600; %-H2_in_GT(i,j)*3600;
                   
                    if H2_stored_aux2_aux > max_storage_surplus   %If we have the storage already full then we decrease the storage of H2 by using it in GT
                        %The H2 to the GT is the amount to keep the storage level within the limit set
                        H2_in_GT(i,j) = (H2_stored_aux2_aux-max_storage_surplus)/3600; %(H2_stored_aux2_aux-0)/3600; %
                        %Need to check we are not sending too much H2 to GT
                        if H2_in_GT(i,j)>max_H2_f;
                            H2_in_GT(i,j)=max_H2_f;
                            %If this happens we need to decrease H2_stored_aux by the H2 that cannot be
                            %produced as cannot be either stored nor used in GT
                            H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_in_GT(i,j)*3600; %
                            if H2_stored_aux2_aux > max_storage_surplus %If we cannot send enough H2 to GT to constrain storage we might need to dissipate wind power
                                H2_produced_old = H2_produced_h(i,j);
                                H2_produced_h(i,j) = max_H2_f; %(H2_stored_aux2_aux-max_storage_surplus)/3600;
                                % Less H2 produced than in 1st iteration --> H2_stored_aux to be adjusted
                                H2_stored_aux_adj = (H2_produced_old - H2_produced_h(i,j))*3600;
                                %Check how much less we use the electrolyser (thus how much more PW dissipation)         
                                H2_while = 0;
                                PW_EL_reduced = 0;                                
                                while H2_while < (H2_produced_old - H2_produced_h(i,j)) %Check how much less we use the electrolyser (thus more PW dissipation)
                                    PW_EL_reduced = PW_EL_reduced+1000; %W
                                    no_EL_on = PW_EL_reduced/10^6/P_ELY_module;
                                    no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                                    no_EL_on = ceil(no_EL_on); %Round to the next integer number
                                    if no_EL_on==0
                                        PW_surplus_cell = 0;
                                    else
                                        PW_surplus_cell = PW_EL_reduced/no_EL_on; %[W]
                                    end                    
                                    H2_while = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                                    H2_while = H2_while*no_EL_on;                                
                                end                                
                                PW_EL_reduced = PW_EL_reduced - 1000; %W
                                PW_EL_reduced = PW_EL_reduced/10^6; %MW
                                PW_wind_h_dissipated(i,j) = PW_wind_h_dissipated(i,j)+PW_EL_reduced;
                                
                            end                             
                            
                        end
                    end

    
                    H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_in_GT(i,j)*3600; %
                    H2_stored_aux = H2_stored_aux - H2_in_GT(i,j)*3600 - H2_stored_aux_adj;
                    H2_stored_h(i,j) = H2_stored_aux2;
                    H2_stored_aux_adj = 0;
                    
                    %Tracking max and min level of storage
                    if H2_stored_aux2 > max_level_storage
                        max_level_storage = H2_stored_aux2;
                    else if H2_stored_aux2 < min_level_storage
                            min_level_storage = H2_stored_aux2;
                        end                    
                    end

                 
                %...else... cases we can reduce power from GT and use the FC so we consume the H2 produced     
                else if H2_stored_aux2 > 0 %We have some produced H2 stored to consume
                                                
                        %Primary approach is to consume H2 by using FC (in some instances FC stack can already be exploited but not to its max capacity)                
                        GT_old = GT_load_act(i,j);  %Original GT load                                       
                        GT_load_act(i,j) = (PW_GT_load_act(i,j) - (size_FC-PW_fc_h(i,j)))/PW_GT_d;
                        if GT_load_act(i,j) < min_GT_l  %We cannot exploit entirely fuel cell if that means overpass the lower limit of GT load
                            GT_load_act(i,j) = min_GT_l;
                        end
 
                        %If we use polynomial interpolation from Marcin's model
                        no_FC_on = PW_fc_h(i,j)/P_FCS_module; %Number of cells of the FC stack to activate
                        no_FC_on = no_FC_on - 0.001;  %To avoid wrong rounding to next integer
                        no_FC_on = ceil(no_FC_on); %Round to the next integer number
                        if no_FC_on >0
                            PW_cell = PW_fc_h(i,j)/no_FC_on*10^6;
                        else
                            PW_cell = 0;
                        end 
                        old_H2_used_fc=A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC); %kg/s
                        old_H2_used_fc = old_H2_used_fc*no_FC_on;                        
                        
%                         %If we use constant FC efficiency
%                         old_H2_used_fc=PW_fc_h(i,j)/perf_fc;
                        
                        %PW produced by FC (could be increased now)      
                        PW_fc_h(i,j) = PW_GT_load_act(i,j) - GT_load_act(i,j)*PW_GT_d + PW_fc_h(i,j);
                        
                        %If we use polynomial interpolation from Marcin's model
                        no_FC_on = PW_fc_h(i,j)/P_FCS_module; %Number of cells of the FC stack to activate
                        no_FC_on = no_FC_on - 0.001;  %To avoid wrong rounding to next integer
                        no_FC_on = ceil(no_FC_on); %Round to the next integer number
                        if no_FC_on >0
                            PW_cell = PW_fc_h(i,j)/no_FC_on*10^6;
                        else
                            PW_cell = 0;
                        end                     
                        H2_used_fc_bis = A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC); %kg/s
                        H2_used_fc_bis = H2_used_fc_bis*no_FC_on;
%                         %If we use constant FC efficiency
%                         H2_used_fc_bis = PW_fc_h(i,j)/perf_fc;
                        
                        H2_stored_aux2_aux = H2_stored_aux2 - H2_used_fc_bis*3600; %(H2_used_fc_bis-old_H2_used_fc)*3600; %
                                                
                        %Check we are not overpassing maximum allowed storage size
                        if H2_stored_aux2_aux > max_storage_surplus %If we have the storage full we can consume some additional H2 in GT
                            H2_in_GT(i,j) = (H2_stored_aux2_aux-max_storage_surplus)/3600;
                            %Need to check we are not sending too much H2 to GT
                            if H2_in_GT(i,j)>max_H2_f;
                                H2_in_GT(i,j)=max_H2_f;
                            end
                        end

                  
                        H2_stored_aux2_aux = H2_stored_aux2 - H2_used_fc_bis*3600 - H2_in_GT(i,j)*3600; %(H2_used_fc_bis-old_H2_used_fc)*3600; %
                        H2_stored_aux = H2_stored_aux - (H2_used_fc_bis-old_H2_used_fc)*3600 - H2_in_GT(i,j)*3600;
                        H2_stored_aux2 = H2_stored_aux2 - H2_used_fc_bis*3600 - H2_in_GT(i,j)*3600;
                        PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j);
                                                
                        %H2_stored_aux = H2_stored_aux - (H2_used_fc_bis-old_H2_used_fc)*3600 - H2_in_GT(i,j)*3600;
                        %H2_stored_aux2 = H2_stored_aux2 - H2_used_fc_bis*3600 - H2_in_GT(i,j)*3600;
                        H2_stored_h(i,j) = H2_stored_aux2;
                        %PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j)                        
                        
                        %Tracking max and min level of storage
                        if H2_stored_aux2 > max_level_storage
                            max_level_storage = H2_stored_aux2;
                        else if H2_stored_aux2 < min_level_storage
                                min_level_storage = H2_stored_aux2;
                            end
                        end
                        
                    else %Not producing (because no excess wind here) and not using H2 (because storage is empty)
                        %But we might need fuel cell stack if the GT+wind is not enough to meet power demand
                        %Or we might want to produce H2 anyway (increasing GT load) if the storage level is too low                        
                        
                        if PW_fc_h(i,j)>0 %Cases where we have to use FC (GT+wind not enough)

                            %If we use polynomial interpolation from Marcin's model
                            no_FC_on = PW_fc_h(i,j)/P_FCS_module; %Number of cells of the FC stack to activate
                            no_FC_on = no_FC_on - 0.001;  %To avoid wrong rounding to next integer
                            no_FC_on = ceil(no_FC_on); %Round to the next integer number
                            if no_FC_on >0
                                PW_cell = PW_fc_h(i,j)/no_FC_on*10^6;
                            else
                                PW_cell = 0;
                            end                     
                            H2_used_fc_bis = A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC); %kg/s
                            H2_used_fc_bis = H2_used_fc_bis*no_FC_on;   

    %                         %If we use constant FC efficiency
    %                         H2_used_fc_bis = PW_fc_h(i,j)/perf_fc;

                            H2_stored_aux2 = H2_stored_aux2 - H2_used_fc_bis*3600 - H2_in_GT(i,j)*3600; 
                            H2_stored_h(i,j) = H2_stored_aux2;                            
                            
                        else if (-H2_stored_aux2_aux) > max_storage_lack %Cases where we want to produce H2 (storage level is too low)                                

                                PW_available=(max_GT_l-GT_load_act(i,j))*PW_GT_d;
                                if PW_available > max_size_el_peak
                                    %PW_wind_h_dissipated(i,j)= PW_available-max_size_el_peak;
                                    PW_available = max_size_el_peak;
                                    GT_load_act(i,j)= (max_size_el_peak+GT_load_act(i,j)*PW_GT_d)/PW_GT_d;                                
                                else
                                    GT_load_act(i,j)=max_GT_l;
                                end

                                %If we use polynomial interpolation from Marcin's model
                                no_EL_on = abs(PW_available/P_ELY_module); %Number of cells of the EL stack to activate
                                no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                                no_EL_on = ceil(no_EL_on); %Round to the next integer number
                                if no_EL_on > 0
                                    PW_surplus_cell = PW_available/no_EL_on*10^6;
                                else
                                    PW_surplus_cell = 0; %[W]
                                end
                                H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                                H2_produced = H2_produced*no_EL_on; %...the H2 produced would be this  
                                H2_stored_aux = H2_stored_aux + H2_produced*3600; %Some more H2 than expected is produced then we need to consume it again to close balance
                                H2_produced_h(i,j) = H2_produced; %+H2_produced_h(i,j)
                                PW_GT_load_act(i,j)=GT_load_act(i,j)*PW_GT_d;
                                H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;                                 
                                H2_stored_h(i,j) = H2_stored_aux2;    
                                
                            else %Not producing and not using H2
                                H2_stored_aux2 = H2_stored_aux2 - H2_in_GT(i,j)*3600; 
                                H2_stored_h(i,j) = H2_stored_aux2;                             
                            end
                        end
                        
                    end
                end
                                        
            else %H2_stored_aux <= 0  We already used all H2 we produced but we still have to keep under control storage
  
                H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600; %
                
                %H2_stored_h(i,j) = H2_stored_aux2;              
                
                if H2_stored_aux2_aux > max_storage_surplus %| (H2_stored_aux2-min_level_storage) > limit_storage %Storage full, then we decrease the storage of H2 by using it in GT
                    %The H2 to the GT is the amount to keep the storage level within the limit set
                    H2_in_GT_old = H2_in_GT(i,j);
                    H2_in_GT(i,j) = (H2_stored_aux2_aux-max_storage_surplus)/3600; %(H2_stored_aux2_aux-0)/3600; %
                    %H2_in_GT_old = H2_in_GT(i,j); 
                    %Need to check we are not sending too much H2 to GT
                    if H2_in_GT(i,j)>max_H2_f;
                        H2_in_GT(i,j)=max_H2_f;
                    end
                    %H2_produced_h(i,j) = H2_produced_h(i,j)-(H2_in_GT_old-H2_in_GT(i,j)); %if we cannot send all surplus of H2 to the GT then we need to decrease accordingly H2 production
                    %PW_wind_h_dissipated(i,j) = ?
                    H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;
                    H2_stored_aux = H2_stored_aux - (H2_in_GT(i,j)-H2_in_GT_old)*3600; %Some less H2 than expected is stored then we need to produce it again to close balance
                    
                else if (-H2_stored_aux2_aux) > max_storage_lack %| (max_level_storage-H2_stored_aux2) > limit_storage %Storage empty, we produce some hydrogen to fill it
                        if GT_load_act(i,j) < max_GT_l %if we are not using the GT at max load we can increase the load to allow EL producing some extra H2
                            PW_available=(max_GT_l-GT_load_act(i,j))*PW_GT_d;
                            if PW_available > max_size_el_peak
                                %PW_wind_h_dissipated(i,j)= PW_available-max_size_el_peak;
                                PW_available = max_size_el_peak;
                                GT_load_act(i,j)= (max_size_el_peak+GT_load_act(i,j)*PW_GT_d)/PW_GT_d;                                
                            else
                                GT_load_act(i,j)=max_GT_l;
                            end
                            
                            %If we use polynomial interpolation from Marcin's model
                            no_EL_on = abs(PW_available/P_ELY_module); %Number of cells of the EL stack to activate
                            no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                            no_EL_on = ceil(no_EL_on); %Round to the next integer number
                            if no_EL_on > 0
                                PW_surplus_cell = PW_available/no_EL_on*10^6;
                            else
                                PW_surplus_cell = 0; %[W]
                            end
                            H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                            H2_produced = H2_produced*no_EL_on; %...the H2 produced would be this  
                            H2_stored_aux = H2_stored_aux + H2_produced*3600; %Some more H2 than expected is produced then we need to consume it again to close balance
                            H2_produced_h(i,j) = H2_produced_h(i,j)+H2_produced;
                            PW_GT_load_act(i,j)=GT_load_act(i,j)*PW_GT_d;
                            H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;                            
                        else
                            H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;
                        end  

                    else
                        H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;                        
                        
                    end
              end
          end  

        end

    end
    %End of second run throughout the year %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
   H2_stored_aux
    
  %Now the storage balances should be closer to zero...
  H2_stored = 0;
  H2_GT_available=0;
  max_size_FC=0;
  %...but not necessarily zero yet
  if H2_stored_aux > 0  %We still have an excess of H2 to use. Here we check when we did not use H2 in GT and we then use it there
     hours_available = 0;
     H2_available_FC = 0;
     for i=1:8760
         if H2_in_GT(i,j)==0 %hours where we do not send H2 to GT
             hours_available = hours_available+1;
         end    
     end

     %If we prioritise use of H2 in GT
     H2_GT_available = H2_stored_aux/hours_available/3600; %kg/s %We equally split the H2 to the GT in all hours available to use some more H2
                                                           %This is the amount of H2 we would feed to GT
     %H2_available_FC = 0;
     if H2_GT_available > max_H2_f || hours_available==0 %Not possible to consume all H2 like this. Could be attempted to increase FC stack size
         statement2 = 'It will not be possible to use all H2'
         H2_GT_available
         size_FC
         max_size_FC
     end
     if H2_GT_available > max_H2_f
         H2_GT_available = max_H2_f; %H2 to GT cannot be higher than maximum limit
     end
  else if H2_stored_aux < 0%We still have a lack of H2. Here we check where this H2 can be produced
      for i=1:8760
          if GT_load_act(8760+1-i,j)<max_GT_l & H2_produced_h(8760+1-i,j)==0 %if we are not using the GT at max load we can increase the load to allow EL producing some extra H2
                           
                PW_available=(max_GT_l-GT_load_act(8760+1-i,j))*PW_GT_d;
                if PW_available > max_size_el_peak
                    PW_available = max_size_el_peak;
                    GT_load_act(8760+1-i,j)= (max_size_el_peak+GT_load_act(8760+1-i,j)*PW_GT_d)/PW_GT_d;                                
                else
                    GT_load_act(8760+1-i,j)=max_GT_l;
                end
                PW_GT_load_act(8760+1-i,j)=GT_load_act(8760+1-i,j)*PW_GT_d;

                %If we use polynomial interpolation from Marcin's model
                no_EL_on = abs(PW_available/P_ELY_module); %Number of cells of the EL stack to activate
                no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                no_EL_on = ceil(no_EL_on); %Round to the next integer number
                if no_EL_on > 0
                    PW_surplus_cell = PW_available/no_EL_on*10^6;
                else
                    PW_surplus_cell = 0; %[W]
                end              
                H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                H2_produced = H2_produced*no_EL_on; %...the H2 produced would be this  
                H2_stored_aux = H2_stored_aux + H2_produced*3600; %Some more H2 than expected is produced then we need to consume it again to close balance
                H2_produced_h(8760+1-i,j)=H2_produced;            
              
          end
          if H2_stored_aux>0
              break
          end
      end
      end %end of --> else if H2_stored_aux < 0%
  end %end of --> if H2_stored_aux > 0 


    %THIRD ITERATION throughout the year - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %In this iteration we implement the secondary strategies to produce H2 or consume H2 (starting point is conditions defined in second iteration)
    %Second attempt to set the storage balances to zero
    for i=1:8760   
        
        %Contribution of wind power
        PW_wind_h(i,1)=PW_norm_h(i,1)*WC;
        
        %Amount of H2 to feed to GT to set storage balance to zero
        if H2_in_GT(i,j) == 0
            H2_in_GT(i,j)=H2_GT_available;
        end
        
%         if PW_fc_h(i,j)<size_FC & GT_load_act(i,j)>min_GT_l & H2_available_FC>0    
%             GT_load_act(i,j)=GT_load_act(i,j)-extra_FC(i,j)/PW_GT_d;
%             PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j);
%             PW_fc_h(i,j)=PW_fc_h(i,j)+extra_FC(i,j);
%             H2_available_FC = H2_available_FC - extra_FC(i,j)/perf_fc*3600;
%         end

        %Remaining power to meet total demand
        PW_remaining = POWER(1,j) - PW_wind_h(i,1) - PW_GT_load_act(i,j) + PW_wind_h_dissipated(i,j); %PW_GT_load_act(i,j) was defined in second iteration

        if PW_remaining >= 0  %We need contribution of fuel cell

            %If we use polynomial interpolation from Marcin's model
            no_FC_on = PW_remaining/P_FCS_module; %Number of cells of the FC stack to activate
            no_FC_on = no_FC_on - 0.001;  %To avoid wrong rounding to next integer
            no_FC_on = ceil(no_FC_on); %Round to the next integer number
            if no_FC_on >0
                PW_remaining_cell = PW_remaining/no_FC_on*10^6; %[W]  
            else
                PW_remaining_cell = 0;
            end                          
            H2_used_fc = A*(C2_FC*PW_remaining_cell^2+C1_FC*PW_remaining_cell+C0_FC); %kg/s
            H2_used_fc = H2_used_fc*no_FC_on;
            
%             %If we use constant FC efficiency
%             H2_used_fc = PW_remaining/perf_fc; %H2 used to produce power with fuel cell
            H2_stored = H2_stored + H2_produced_h(i,j)*3600 - H2_used_fc*3600 - H2_in_GT(i,j)*3600;
            %H2_stored = H2_stored - H2_used_fc*3600 - H2_in_GT(i,j)*3600; %H2 stored (kg)
            H2_stored_h(i,j) = H2_stored;
            PW_fc_h(i,j) = PW_remaining;           
                      
        else  %PW_remaining < 0 --> We need to produce H2 not to waste wind power

%             %If we use polynomial interpolation from Marcin's model
%             no_EL_on = abs(PW_remaining/P_ELY_module); %Number of cells of the EL stack to activate
%             no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
%             no_EL_on = ceil(no_EL_on); %Round to the next integer number
%             if no_EL_on > 0
%                 PW_surplus_cell = -PW_remaining/no_EL_on*10^6;
%             else
%                 PW_surplus_cell = 0; %[W]
%             end
%             H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
%             H2_produced = H2_produced*no_EL_on; %...the H2 produced would be this 
%             
% %             %If we use constant EL efficiency
% %             H2_produced = -1*PW_remaining*perf_el;
% 
%             H2_stored = H2_stored + H2_produced*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;
%             %H2_stored = H2_stored + H2_produced*3600 - H2_in_GT(i,j)*3600;
%             H2_stored_h(i,j) = H2_stored;
%             H2_produced_h(i,j) = H2_produced;
            H2_stored = H2_stored + H2_produced_h(i,j)*3600 - H2_used_fc_h(i,j)*3600 - H2_in_GT(i,j)*3600;
            H2_stored_h(i,j) = H2_stored;
            
        end
    end
    
    %Re-check size FC
    size_FC_aux = max(PW_fc_h(:,j)); %Size of FC in MW
    if size_FC_aux > size_FC
       size_FC = size_FC_aux;
    end
    
    %Size of electrolyser
    size_electrolyser_aux = max(H2_produced_h(:,j)); %Size of FC in kg/s of H2
    
    H2_while = 0;
    PW_EL_size = size_electrolyser*10^6; %W
    while H2_while < size_electrolyser_aux %Check how much we can use the electrolyser to produce H2 that is used in the GT
        PW_EL_size = PW_EL_size+1000; %W
        no_EL_on = PW_EL_size/10^6/P_ELY_module;
        no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
        no_EL_on = ceil(no_EL_on); %Round to the next integer number
        if no_EL_on==0
            PW_surplus_cell = 0;
        else
            PW_surplus_cell = PW_EL_size/no_EL_on; %[W]
        end
        H2_while = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
        H2_while = H2_while*no_EL_on;                                
    end
    size_electrolyser_aux = (PW_EL_size-1000)/10^6; %MW
    
    if size_electrolyser_aux > size_electrolyser
       size_electrolyser = size_electrolyser_aux;
    end    
    
    %Size energy storage
    size_storage_max = max(H2_stored_h(:,j));
    size_storage_min = min(H2_stored_h(:,j));
    size_storage_aux = size_storage_max-size_storage_min; %Size energy storage in kg H2
    if size_storage_aux > size_storage
       size_storage = size_storage_aux;
    end       
    
    %OUTPUTS
    %H2 storage size (kg)
    size_storage
    %Size electrolyser stack (MW)
    size_electrolyser
    %Size fuel cell stack (MW)
    size_FC
    
    %Write results in Excel file
    filename = 'Outputs.xlsm';
    xlswrite(filename,PW_wind_h,j,'C4')
    xlswrite(filename,PW_wind_h_dissipated(:,j),j,'I4')
    xlswrite(filename,H2_stored_h(:,j),j,'E4')
    xlswrite(filename,PW_fc_h(:,j),j,'F4')
    xlswrite(filename,H2_produced_h(:,j),j,'G4')
    xlswrite(filename,max_size_el_peak,4,'D26')
    xlswrite(filename,PW_GT_load_act(:,j),j,'U4')
    xlswrite(filename,GT_load_act(:,j),j,'V4')
    xlswrite(filename,H2_in_GT(:,j),j,'X4')
    xlswrite(filename,POWER(j),j,'A1')
    xlswrite(filename,Num_GT(:,j),j,'T4')
      
end


%
%% CO2 emissions calculations
% % The precise CO2 emissions are given by the GT model in Excel given a energy input needed for each GT load
% % The approach here is to evaluate the energy input per each GT load instance...
% % ...then check what part of that energy input is covered by H2 (no CO2 emissions) or by natural gas
% % At this point calculate the related CO2 emissions
% 
% %     %Set the H2 to zero in the GT model (in excel) as initial approximation
% %     % Makes the 8th Excel sheet active.
% %     eSheet1 = eSheets.get('Item', 8); %get 8th sheet from eSheets. eSheets contain the Excel sheets
% %     eSheet1.Activate; %We activate the 8th sheet of the open excel file
% %     %H2 in GT
% %     e_loadH2 = get(e.Activesheet,'Range', 'E31');
% %     e_loadH2.Value = 0; %y_H2;
%     
% %  Check how many different GT load instances are there per year for which the CO2 emissions have to be calculated
% for j=1:3
%     
%     GT_load_instances = unique(GT_load_act(:,j));
%     instances = length(GT_load_instances);
%     CO2_instances = zeros(instances,1);
%     en_in_instances = zeros(instances,1);
%     
%     %For each insatnce identified, retrieve the energy input from GT model in excel
%     for i=1:instances
%         
% %          %Set the correct GT load in the GT model and read the outputs
% %          % Makes the 8th Excel sheet active.
% %          eSheet1 = eSheets.get('Item', 8); %get 8th sheet from eSheets. eSheets contain the Excel sheets
% %          eSheet1.Activate; %We activate the 8th sheet of the open excel file
% %          %Write the load of the GT
% %          e_loadGT = get(e.Activesheet,'Range', 'I13');
% %          e_loadGT.Value = GT_load_instances(i);
% %          %Read energy input
% %          e_en_in = get(e.Activesheet,'Range', 'E45');
% %          en_in_instances(i,1) = e_en_in.Value;
%          
%          % %Energy input from polynomial interpolation
%          if GT_type == 1
%              % Around GT load 50% the ploynomial interpolation overestimtes the EN in by about 2%
%              en_in_instances(i,1)=-8610.9*GT_load_instances(i)^6+34428*GT_load_instances(i)^5-55262*GT_load_instances(i)^4+45432*GT_load_instances(i)^3-20146*GT_load_instances(i)^2+4639.9*GT_load_instances(i)-394.54;
%          else
%              en_in_instances(i,1)=-369.45*GT_load_instances(i)^4+999.08*GT_load_instances(i)^3-953.2*GT_load_instances(i)^2+447.66*GT_load_instances(i)-17.295;
%          end
%         
%     end
%     
%     %Associate each of the 8760 instances in a year with an energy input 
%      for i=1:8760
%          for ii=1:instances
%             if GT_load_act(i,j)==GT_load_instances(ii)
%                  en_in(i,j) = en_in_instances(ii,1);
%             end             
%          end 
%      end
%     
%      %Now the contribution of H2 and natural gas to the energy input needs to be evaluated
%      for i=1:8760        
%          en_in_H2 = H2_in_GT(i,j)*LHV_H2; %depending on how much H2 is fed to GT, we estimate the related energy input
%          en_in_NG = en_in(i,j) - en_in_H2; %the energy input from NG will be the total energy input minus the energy input from H2
%          mm_NG = en_in_NG/LHV_NG; %we can then calculate the amount of NG needed and thus the CO2 emissions (complete combustion assumption)
%          CO2_act(i,j) = mm_NG*2.7433;  %2.7433 is the kgCO2 per kgCH4 (CO2 in air is neglected)
%          % If a natural gas is used instead of CH4 the factor needs to be changed!
%          if Num_GT(i,j) == 2
%              CO2_act(i,j) = CO2_act(i,j)*2;
%          end                  
%      end
% 
%     %OUTPUTS
%     %CO2 emissions (Mt)
%     CO2_year=mean(CO2_act(:,j))*365*24*3600/10^9
%     
%         %Write the results in excel
%         xlswrite(filename,CO2_act(:,j),j,'W4')
% end
% % 

%%
%COMMUNICATION WITH SHORT TERM MODEL

% %define how many GT are used 
% Num_GT_op = 1;
% for i=1:8760
%     if Num_GT(i,1)==2
%         Num_GT_op = 2;
%         break
%     end
% end
% 
% %Number of wind turbines
% if WT_type == 1
%     WT_n = WC/5;
% else
%    WT_n = WC/6;
% end
% 
% %Negative of electrolyser size
% size_electrolyser_bis = -size_electrolyser;
% 
% %Continuous load rated power - should become an array !?!?!??!
% Cont_load = POWER(1)/POWER(1);
% 
% %INPUTS:
% %GT_type --> GT_type (1=LM2500; 2=LM6000)
% %Num_GT_op --> number of Gas Turbines
% %min_GT_l --> minimum setpoint [pu]
% %max_GT_l --> maximum setpoint [pu]
% %WT_type --> Wind Turbine type (1 for NREL 5MW, 2 for Hywind 6MW)
% %WT_n --> number of Wind Turbines
% %size_electrolyser_bis --> size_electrolyser [MW] / electrolizer rated power (<= 0) [pu]
% %size_FC --> size_FC [MW] / fuel cell rated power (>= 0) [pu]
% %Cont_load --> continuous load rated power - should become an array !?!?!??!
% %POWER(1) --> base value apparent power [MVA]
% 
% %OUTPUTS
% %frequency variation within limits (Y/N) --> frequency_boolean
% %load ramp speed GT (Y/N) --> ramp_boolean
% %load ramp speed value --> ramp
% 
% [frequency_boolean,ramp_boolean,ramp] = short_term_analysis(GT_type,Num_GT_op,min_GT_l,max_GT_l,WT_type,WT_n,size_electrolyser_bis,size_FC,Cont_load,POWER(1));
% 
% %OUTPUTS
% %Frequency variation within limits
% frequency_boolean
% %Load ramp speed GT within limits
% ramp_boolean
% %Ramp speed GT within limits
% ramp

