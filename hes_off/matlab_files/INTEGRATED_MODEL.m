%%
%Script for GT model

global e workbook file eSheets eSheet1 GT_type

%Call "Excel_GT.m" file first
%That file while activate the excel file with the GT model
%The type of GT to be used is chosen there as well 
%GT_type = 1 --> LM2500
%GT_type = 2 --> LM6000

%%
%Script for wind power

%Wind turbine to use (power curve)
Power_curve = Hywind;
WT_type = 2;
% Power_curve = NREL;
% WT_type = 1;

%Data series of wind speed instances throughout one year (one minute resolution)
wdata = load('sleipnerwind.mat');
wspeed = wdata.wind(1).w;

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
for i=1:8760
   for j=1:60       
    PW_norm_h(i,1) = PW_norm_h(i,1)+PW_norm(j+60*(i-1));    
   end
   PW_norm_h(i,1)=PW_norm_h(i,1)/60;
end

%%
%List of assumptions

%POWER DEMAND
%Three cases considered: peak demand, mid-life demand and tail demand
%Power demand (MW)
POWER = [43.6, 35.2, 32.9]; %[45, 35, 30]; %

%WIND FARM
%Wind capacity installed (MW)
WC = 24;    %-----> Decision variable <-----%

%GAS TURBINE
%Design GT power (MW) (power provided by the GT at design - 100% load)
if GT_type == 1
    PW_GT_d = 33.31; %LM2500
else
    PW_GT_d = 44.70; %LM6000
end

%Possibility to start a second GT %1 yes // 2 no
if GT_type == 1
    GT2=1;  %LM2500
else
    GT2=2; %LM6000
end

%Variable that keeps tracks of #GT operating at any instance
Num_GT = ones(8760,3);
%Max GT load
max_GT_l = 0.95; %1.00;     %-----> Decision variable <-----%
%Min GT load (GT cannot be operated at a lower load)
min_GT_l = 0.40;     %-----> Decision variable <-----%
% %Fixed GT load (if we decided to have it constant)
% GT_load = [0.95, 0.70, 0.63]; %0.95; %0.70; 0.63;
%Variable that keeps tracks of actual GT loads
GT_load_act = zeros(8760,3);
%Variable that keeps tracks of actual GT power
PW_GT_load_act = zeros(8760,3);
% %Related GT power output at the fixed load
% PW_GT_load = GT_load * PW_GT_d;
% %CO2 emissions from GT with and without H2(kg/s)
% CO2_GT = zeros(1,3);
% CO2_H2inGT = zeros(1,3);
%Variable that keeps tracks of actual GT emissions
CO2_act = zeros(8760,3);
%Variable that keeps tracks of energy input
en_in = zeros(8760,3);
%Max fraction of H2 (volumetric) that can be used in GT  %-----> Decision variable <-----%
if GT_type == 1
    max_H2_f = 0.032; %0.015 kg/s is about 10%vol H2 inlet at GT load 40 -- 0.018 kg/s about 12% -- 0.032 kg/s about 20%
else
    max_H2_f = 0.037; %0.017 kg/s is about 10%vol H2 inlet at GT load 40 -- 0.02 kg/s about 12% -- 0.037 kg/s about 20%
end

%ELECTROLYSER
%Electrolyser performance (kgH2/MJ)
perf_el = 0.8*1/120; %0.005; %
%Max size electrolsyer in peak years
max_size_el_peak = 4; %MW  %-----> Decision variable <-----%
%Size elctrolyser
size_electrolyser = 0; %MW %Initialised
%Size single EL in a stack
P_ELY_module = 0.500;   %MW; 
%Coefficients for H2_produced
C3=-1.41175e-20;
C2 = 1.59861e-14;
C1 = - 8.29304e-09; 
C0 = 6.82299e-03; 
%Max number of module in stacks
max_modules_el = max_size_el_peak/P_ELY_module;

%FUEL CELL
%Fuel cell performance (MJ/kgH2)
perf_fc = 66; %53.4; %0.5*120;
%Size fuel cell (initialized)
size_FC = 3;
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
max_storage_half = 35000; %0; %99999999999999; %70800; %17750; %
max_storage = 5000;  %-----> Decision variable <-----%
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
%the hydrogen produced by electrolyser
H2_produced_h = zeros(8760,3);    
%the storage level
H2_stored_h = zeros(8760,3);
H2_stored_h_aux = zeros(8760,1);
%power used from fuel cell
PW_fc_h = zeros(8760,3);
%extra power to harvest from fuel cell
extra_FC= zeros(8760,3);

%CO2 EMISSIONS
%Cumulative CO2 emissions
CO2_tot = 0;

%%
%ITERATIONS throghout the year for the three power demands
for j=1:3 %Three cases considered: peak demand, mid-life demand and tail demand
    j
    
    %Define some auxiliary variables
    H2_stored = 0; %9999999;       
    H2_stored_aux = 0;
    H2_stored_aux2 = 0;
    H2_stored_aux2_aux=0;
    
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
                no_FC_on = ceil(no_FC_on); %Round to the next integer number
                if no_FC_on >0
                    PW_deficit_cell = PW_deficit/no_FC_on*10^6; %[W]  
                else
                    PW_deficit_cell = 0;
                end     
                H2_used_fc = A*(C2_FC*PW_deficit_cell^2+C1_FC*PW_deficit_cell+C0_FC); %kg/s
                H2_used_fc = H2_used_fc*no_FC_on;
                
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
                
%                 %Has to be zero here for the algorithm implemented later%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 H2_in_GT(i,j) = 0; %0.018; %kg/s  %0.018 kg/s is about 10%vol H2 inlet at GT load 40%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

    %Auxiliary variables for H2 storage
    H2_stored_aux = H2_stored;  %H2_stored_aux keeps tracks of the H2 still to produce or consume to reach zero balance
                                %The starting point is the net balance throughout a year calculated in first iteration (H2_stored)
    %H2_stored_h(:,j) = H2_stored_h_aux(:,1);
    H2_stored_aux2 = 0;  %H2_stored_aux2 keeps tracks of the H2 storage level with the strategy used in this second iteration
                         %The starting point is zero

    if H2_stored < 0  %During the first iteration we registered a net negative balance of H2
                      %Lack of power (typical of peak years). Need to produce the H2 to be used in the fuel cell
                      %The primary approach is to increase GT load when possible and to use extra power in electrolyser to produce H2
        
        for i=1:8760                 
             if H2_stored_aux < 0  %This to keep track that we haven't yet produced all the H2 needed
                if GT_load_act(i,j) == max_GT_l  %This are the cases when we do not have enough power and we need to use fuel cell
                                                 %Here we cannot produce H2
                    GT_load_act(i,j) = GT_load_act(i,j);
                    
                    %If we use polynomial interpolation from Marcin's model
                    no_FC_on = PW_fc_h(i,j)/P_FCS_module; %Number of cells of the FC stack to activate
                    no_FC_on = ceil(no_FC_on); %Round to the next integer number
                    if no_FC_on >0
                        PW_cell = PW_fc_h(i,j)/no_FC_on*10^6; %[W]  
                    else
                        PW_cell = 0;
                    end 
                    H2_stored_aux2 = H2_stored_aux2 - A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC)*3600*no_FC_on; %kg
                    
                else  %Here there is enough power and we can produce some H2 to be later used

                    GT_old = GT_load_act(i,j);  %Original GT load                   
                    GT_load_act(i,j) = max_GT_l;  %If we set GT at max load...
                    
                    %If we use polynomial interpolation from Marcin's model
                    no_EL_on = PW_GT_d*(GT_load_act(i,j)-GT_old)/P_ELY_module; %Number of cells of the EL stack to activate
                    no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                    no_EL_on = ceil(no_EL_on); %Round to the next integer number
                    if no_EL_on==0
                        PW_surplus_cell = 0; %[W]
                    else
                        PW_surplus_cell = PW_GT_d*(GT_load_act(i,j)-GT_old)/no_EL_on*10^6; %[W]
                    end                    
                    H2_produced_bis = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                    H2_produced_bis = H2_produced_bis*no_EL_on; %...the H2 produced would be this

                    H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_bis*3600; %the storage tank would fill until this level
%  
                    %Check we are not producing more H2 than what allowed by electrolyser stack
                    %If we use polynomial interpolation from Marcin's model
                    if no_EL_on > max_modules_el  %we would need bigger electrolyser stack...
                        GT_load_act(i,j)=(max_size_el_peak+GT_old*PW_GT_d)/PW_GT_d; %...thus we reduced H2 produced
                        no_EL_on = PW_GT_d*(GT_load_act(i,j)-GT_old)/P_ELY_module; %Number of cells of the EL stack to activate
                        no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                        no_EL_on = ceil(no_EL_on); %Round to the next integer number
                        if no_EL_on==0
                            PW_surplus_cell = 0;
                        else
                            PW_surplus_cell = PW_GT_d*(GT_load_act(i,j)-GT_old)/no_EL_on*10^6; %[W]
                        end                    
                        H2_produced_bis = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                        H2_produced_bis = H2_produced_bis*no_EL_on;
                        H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_bis*3600; %the storage tank would fill until this level
                    end
                  
%                                         
                    %This constrain can be possibly removed
                    %Check we are not overpassing maximum allowed storage size
                    %If we use polynomial interpolation from Marcin's model
                    if H2_stored_aux2_aux > max_storage  %we reached max storage level thus we need to decrease the production of H2
                        no_EL_on = PW_GT_d/P_ELY_module;
                        no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                        no_EL_on = ceil(no_EL_on); %Round to the next integer number
                        if no_EL_on==0
                            GT_load_act(i,j)=GT_old;
                        else
                            PW_surplus_cell = PW_GT_d/no_EL_on*10^6; %[W]
                            H2_produced_aux_d = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                            H2_produced_aux_d = H2_produced_aux_d*no_EL_on;
                            GT_load_act(i,j)=(H2_produced_bis-(abs(H2_stored_aux2_aux)-max_storage)/3600)/(H2_produced_aux_d)+GT_old;  %GT load set to the level that leads to max storage capacity

                        end
                            if GT_load_act(i,j) < GT_old
                                GT_load_act(i,j)=GT_old;
    %                              here1 = 'something wrong'
    %                             pause
                            end                        
                    end

                    %If we use polynomial interpolation from Marcin's model
                    no_EL_on = PW_GT_d*(GT_load_act(i,j)-GT_old)/P_ELY_module; %Number of cells of the EL stack to activate
                    no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
                    no_EL_on = ceil(no_EL_on); %Round to the next integer number
                    if no_EL_on==0
                       PW_surplus_cell = 0;
                    else
                        PW_surplus_cell = PW_GT_d*(GT_load_act(i,j)-GT_old)/no_EL_on*10^6; %[W]
                    end
                    H2_produced_bis = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
                    H2_produced_bis = H2_produced_bis*no_EL_on; %...the H2 produced would be this                    
                    
                    H2_stored_aux2 = H2_stored_aux2 + H2_produced_bis*3600;
                    H2_stored_aux = H2_stored_aux + H2_produced_bis*3600;
                    H2_stored_h(i,j) = H2_stored_aux2;
                    PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j);

                end
             end              
           end
        
        H2_stored = H2_stored_aux;
        
        if abs(H2_stored)>max_storage
           here2 = 'we could not constrain the size of the H2 storage'
           %pause
        end
        
    else  %During the first iteration we registered a net positive balance of H2
          %Excess of power (typical of tail years). Need to use the H2 produced
          %The primary approach is to use the H2 in the FC stack rather than in the GT (better efficiency)

        for i=1:8760        
            if H2_stored_aux > 0  %This to check that we did not already used all H2 we produced
                if GT_load_act(i,j) == min_GT_l  %These are the cases when we cannot produce less from GT and we need to store some H2
                                                 %Here we can use H2 only in GT (which we do only if we are already over max storage size)
                    
                    H2_stored_aux2_aux = H2_stored_aux2 + H2_produced_h(i,j)*3600; %-H2_in_GT(i,j)*3600;
                                        
                   
                    if H2_stored_aux2_aux > max_storage %> 0 % %  %If we have the storage already full then we decrease the production of H2 by using it in GT
                                                                  %Another possibility is to start already when the storage level is > 0
                        %The H2 to the GT is the amount to keep the storage level within the limit set
                        H2_in_GT(i,j) = (H2_stored_aux2_aux-max_storage)/3600; %(H2_stored_aux2_aux-0)/3600; %
                        %Need to check we are not sending too much H2 to GT
                        if H2_in_GT(i,j)>max_H2_f;
                            H2_in_GT(i,j)=max_H2_f;
                        end
                    end
                    
                    H2_stored_aux2 = H2_stored_aux2 + H2_produced_h(i,j)*3600 - H2_in_GT(i,j)*3600; %
                    H2_stored_aux = H2_stored_aux - H2_in_GT(i,j)*3600;
                    H2_stored_h(i,j) = H2_stored_aux2;
                 
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
                        no_FC_on = ceil(no_FC_on); %Round to the next integer number
                        if no_FC_on >0
                            PW_cell = PW_fc_h(i,j)/no_FC_on*10^6;
                        else
                            PW_cell = 0;
                        end 
                        old_H2_used_fc=A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC); %kg/s
                        old_H2_used_fc = old_H2_used_fc*no_FC_on;
                              
                        PW_fc_h(i,j) = PW_GT_load_act(i,j) - GT_load_act(i,j)*PW_GT_d + PW_fc_h(i,j);

                        %If we use polynomial interpolation from Marcin's model
                        no_FC_on = PW_fc_h(i,j)/P_FCS_module; %Number of cells of the FC stack to activate
                        no_FC_on = ceil(no_FC_on); %Round to the next integer number
                        if no_FC_on >0
                            PW_cell = PW_fc_h(i,j)/no_FC_on*10^6;
                        else
                            PW_cell = 0;
                        end                     
                        H2_used_fc_bis = A*(C2_FC*PW_cell^2+C1_FC*PW_cell+C0_FC); %kg/s
                        H2_used_fc_bis = H2_used_fc_bis*no_FC_on;
                        
                        H2_stored_aux2_aux = H2_stored_aux2 - H2_used_fc_bis*3600;
                                                
                        %Check we are not overpassing maximum allowed storage size
                        %Another possibility is to start already when the storage level is > 0
                        if H2_stored_aux2_aux > max_storage %> 0 %  %If we have the storage full we can consume some additional H2 in GT
                            H2_in_GT(i,j) = (H2_stored_aux2_aux-max_storage)/3600; %(H2_stored_aux2_aux-0)/3600; %
                            %Need to check we are not sending too much H2 to GT
                            if H2_in_GT(i,j)>max_H2_f;
                                H2_in_GT(i,j)=max_H2_f;
                            end
                        end
                        
                        H2_stored_aux = H2_stored_aux - (H2_used_fc_bis-old_H2_used_fc)*3600 - H2_in_GT(i,j)*3600;
                        H2_stored_aux2 = H2_stored_aux2 - H2_used_fc_bis*3600 - H2_in_GT(i,j)*3600;
                        H2_stored_h(i,j) = H2_stored_aux2;
                        PW_GT_load_act(i,j) = PW_GT_d*GT_load_act(i,j);
                        
                    else %Not producing (because no excess wind here) and not using H2 (because storage is empty)
                        %TO CHECK. MAYBE I CAN STILL CONSUME H2 EVEN THOUGH THE STORAGE IS TEMPORARILY EMPTY
                        %H2_stored_aux2 = H2_stored_aux2; 
                        H2_stored_h(i,j) = H2_stored_aux2;
                    end
                end
                
            end        
        end

    end
    %End of second run throughout the year %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
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
%          if PW_fc_h(i,j)>0 & GT_load_act(i,j)>min_GT_l %There is room to use more FC           size_FC
%              max_size_FC_aux = PW_fc_h(i,j)+(GT_load_act(i,j)-min_GT_l)*PW_GT_d;
%              if max_size_FC_aux > max_size_FC
%                  max_size_FC = max_size_FC_aux;
%              end
%          end         
%          
%          if PW_fc_h(i,j)<size_FC & GT_load_act(i,j)>min_GT_l %There is room to use more FC
%              extra_FC(i,j) = (GT_load_act(i,j)-min_GT_l)*PW_GT_d;  %PW_fc_h(i,j)+
%              if (extra_FC(i,j) + PW_fc_h(i,j)) > size_FC
%                  extra_FC(i,j) = size_FC - PW_fc_h(i,j);
%              end
%              H2_available_FC = H2_available_FC + extra_FC(i,j)/perf_fc*3600;
%          end
%         %hours_available = hours_available+(max_H2_f-H2_in_GT(i,j))/max_H2_f;     
     end

%      %If we prioritise use of H2 in FC
%      if H2_available_FC > H2_stored_aux %using FC we would consume more H2 than needed
%          H2_available_FC = H2_stored_aux;  %kg of H2 we can use in FC
%          H2_GT_available = 0;
%      else %using FC we would still need to use some H2 in GT
%          H2_GT_available = (H2_stored_aux - H2_available_FC)/hours_available/3600; %kg/s
%          if H2_GT_available > max_H2_f
%              statement2 = 'Will not be able to use all H2'
%              H2_GT_available = max_H2_f;
%          end        
%      end
         
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
          
 end

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
        PW_remaining = POWER(1,j) - PW_wind_h(i,1) - PW_GT_load_act(i,j); %PW_GT_load_act(i,j) was defined in second iteration
                
        if PW_remaining > 0  %We need contribution of fuel cell

            %If we use polynomial interpolation from Marcin's model
            no_FC_on = PW_remaining/P_FCS_module; %Number of cells of the FC stack to activate
            no_FC_on = ceil(no_FC_on); %Round to the next integer number
            if no_FC_on >0
                PW_remaining_cell = PW_remaining/no_FC_on*10^6; %[W]  
            else
                PW_remaining_cell = 0;
            end                          
            H2_used_fc = A*(C2_FC*PW_remaining_cell^2+C1_FC*PW_remaining_cell+C0_FC); %kg/s
            H2_used_fc = H2_used_fc*no_FC_on;
            
            H2_stored = H2_stored - H2_used_fc*3600 - H2_in_GT(i,j)*3600; %H2 stored (kg)
            H2_stored_h(i,j) = H2_stored;
            PW_fc_h(i,j) = PW_remaining;
                      
        else  %PW_remaining < 0 --> We need to produce H2 not to waste wind power
            
            %If we use polynomial interpolation from Marcin's model
            no_EL_on = abs(PW_remaining/P_ELY_module); %Number of cells of the EL stack to activate
            no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
            no_EL_on = ceil(no_EL_on); %Round to the next integer number
            if no_EL_on > 0
                PW_surplus_cell = -PW_remaining/no_EL_on*10^6;
            else
                PW_surplus_cell = 0; %[W]
            end
            H2_produced = PW_surplus_cell*(C3*PW_surplus_cell^3+C2*PW_surplus_cell^2+C1*PW_surplus_cell+C0)/10^6; %kg/s
            H2_produced = H2_produced*no_EL_on; %...the H2 produced would be this   
            
            H2_stored = H2_stored + H2_produced*3600 - H2_in_GT(i,j)*3600;
            H2_stored_h(i,j) = H2_stored;
            H2_produced_h(i,j) = H2_produced;
            
        end
    end
    
    %Re-check size FC
    size_FC_aux = max(PW_fc_h(:,j)); %Size of FC in MW
    if size_FC_aux > size_FC
       size_FC = size_FC_aux;
    end
    
    %Size of electrolyser
    size_electrolyser_aux = max(H2_produced_h(:,j)); %Size of FC in kg/s of H2
    %If we use polynomial interpolation from Marcin's model
    P_ELY_module_W=P_ELY_module*10^6; %[W]
    H2_produced_cell = P_ELY_module_W*(C3*P_ELY_module_W^3+C2*P_ELY_module_W^2+C1*P_ELY_module_W+C0)/10^6; %kg/s
    no_EL_on = size_electrolyser_aux/H2_produced_cell; %Number of cells of the EL stack to activate to meet max H2 production
    no_EL_on = no_EL_on - 0.001;  %To avoid wrong rounding to next integer
    no_EL_on = ceil(no_EL_on); %Round to the next integer number
    
    size_electrolyser_aux = P_ELY_module*no_EL_on;    
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



% %CO2 emissions calculations
% %The precise CO2 emissions are given by the GT model in Excel given a energy input needed for each GT load
% %The approach here is to evaluate the energy input per each GT load instance...
% %...then check what part of that energy input is covered by H2 (no CO2 emissions) or by natural gas
% %At this point calculate the related CO2 emissions
% 
%     %Set the H2 to zero in the GT model (in excel) as initial approximation
%     % Makes the 8th Excel sheet active.
%     eSheet1 = eSheets.get('Item', 8); %get 8th sheet from eSheets. eSheets contain the Excel sheets
%     eSheet1.Activate; %We activate the 8th sheet of the open excel file
%     %H2 in GT
%     e_loadH2 = get(e.Activesheet,'Range', 'E31');
%     e_loadH2.Value = 0; %y_H2;
%     
%  %Check how many different GT load instances are there per year for which the CO2 emissions have to be calculated
% for j=1:3
%     
%     GT_load_instances = unique(GT_load_act(:,j));
%     instances = length(GT_load_instances)
%     CO2_instances = zeros(instances,1);
%     en_in_instances = zeros(instances,1);
%     
%     %For each insatnce identified, retrieve the energy input from GT model in excel
%     for i=1:instances
%         
%          %Set the correct GT load in the GT model and read the outputs
%          % Makes the 8th Excel sheet active.
%          eSheet1 = eSheets.get('Item', 8); %get 8th sheet from eSheets. eSheets contain the Excel sheets
%          eSheet1.Activate; %We activate the 8th sheet of the open excel file
%          %Write the load of the GT
%          e_loadGT = get(e.Activesheet,'Range', 'I13');
%          e_loadGT.Value = GT_load_instances(i);
%          %Read energy input
%          e_en_in = get(e.Activesheet,'Range', 'E45');
%          en_in_instances(i,1) = e_en_in.Value;
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
%      %OUTPUTS
%      %Cumulative CO2 emissions (Mt)
%      if j==1
%          CO2_tot = CO2_tot+mean(CO2_act(:,j))*3600*24*365*2/10^9;
%      else if j==2
%              CO2_tot = CO2_tot+mean(CO2_act(:,j))*3600*24*365*4/10^9;
%          else
%              CO2_tot = CO2_tot+mean(CO2_act(:,j))*3600*24*365*12/10^9;
%          end         
%      end
%              
% %         %Write the results in excel
% %         xlswrite(filename,CO2_act(:,j),j,'W4')
% 
% end
% 
% %OUTPUTS
% %Cumulative CO2 emissions (Mt)
% CO2_tot


%%
%COMMUNICATION WITH SHORT TERM MODEL

%define how many GT are used 
Num_GT_op = 1;
for i=1:8760
    if Num_GT(i,1)==2
        Num_GT_op = 2;
        break
    end
end

%Number of wind turbines
if WT_type == 1
    WT_n = WC/5;
else
   WT_n = WC/6;
end

%Negative of electrolyser size
size_electrolyser_bis = -size_electrolyser;

%Continuous load rated power - should become an array !?!?!??!
Cont_load = POWER(1)/POWER(1);

%INPUTS:
%GT_type --> GT_type (1=LM2500; 2=LM6000)
%Num_GT_op --> number of Gas Turbines
%min_GT_l --> minimum setpoint [pu]
%max_GT_l --> maximum setpoint [pu]
%WT_type --> Wind Turbine type (1 for NREL 5MW, 2 for Hywind 6MW)
%WT_n --> number of Wind Turbines
%size_electrolyser_bis --> size_electrolyser [MW] / electrolizer rated power (<= 0) [pu]
%size_FC --> size_FC [MW] / fuel cell rated power (>= 0) [pu]
%Cont_load --> continuous load rated power - should become an array !?!?!??!
%POWER(1) --> base value apparent power [MVA]

%OUTPUTS
%frequency variation within limits (Y/N) --> frequency_boolean
%load ramp speed GT (Y/N) --> ramp_boolean
%load ramp speed value --> ramp

[frequency_boolean,ramp_boolean,ramp] = short_term_analysis(GT_type,Num_GT_op,min_GT_l,max_GT_l,WT_type,WT_n,size_electrolyser_bis,size_FC,Cont_load,POWER(1));

%OUTPUTS
%Frequency variation within limits
frequency_boolean
%Load ramp speed GT within limits
ramp_boolean
%Ramp speed GT within limits
ramp


