
function [frequency_boolean,ramp_boolean,ramp] = short_term_analysis(GT_type,Num_GT_op,min_GT_l,max_GT_l,WT_type,WT_n,size_electrolyser_bis,size_FC,Cont_load,POWER)



%% Wind data - from field measurements or a database
wraw = load('sleipnerwind.mat');
Fs = 1/60;              % Sample frequency of the data [Hz]

% Correction of wind speed for height
href = 14;              % Reference height
hhub = 90;              % Hub height
h0 = 0.0002;            % Roughness lenght - Class 0 (Water surface), according to the European Wind Atlas
wcorr = wraw.wind(1).w * log(hhub/h0) / log(href/h0);

% Calculate the derivative of the wind speed
wdiff = diff(wcorr);

% Find maximum derivative
[diffMax,idxMax] = max(wdiff);
wdata.max = timeseries(wcorr(idxMax-1:idxMax+2),0:1/Fs:3/Fs);

% Find minimum derivative
[diffMin,idxMin] = min(wdiff);
wdata.min = timeseries(wcorr(idxMin-1:idxMin+2),0:1/Fs:3/Fs);

%% Load data - SCADA data provided by Lundin
GTdata = load('GTLoadMin.mat');
GTLoad = GTdata.GTLoad; %fillmissing(GTdata.GTLoad,'linear');

% Calculate the derivative of the load
ldiff = diff(GTLoad);

% Find maximum derivative
[diffMax,idxMax] = max(ldiff);
ldata.On = diffMax;

% Find minimum derivative
[diffMin,idxMin] = min(ldiff);
ldata.Off = diffMin;

%% Simulation parameters

% Basic data
par.elGrid.model = 'average_model_WT_2015a.slx';  % File name of the model to be evaluated
par.elGrid.Sbase = POWER; % 44.7; %               % Base value apparent power [MVA]    
par.elGrid.fbase = 60;                  % Base value frequency [Hz]
par.elGrid.fmax = 1.02;                 % Maximum allowed frequency in the electrical grid [pu]
par.elGrid.fmin = 0.98;                 % Minimum allowed frequency in the electrical grid [pu]

%% Input data

% Gas turbine
in.GT.type = GT_type; % 2; %                       % Gas Turbine type: 1 for LM2500, 2 for LM6000
in.GT.n = Num_GT_op; %  1; %                          % Number of Gas Turbines
in.GT.Pmin = min_GT_l; % 0.3; %                      % Minimum setpoint [pu]
in.GT.Pmax = max_GT_l; % 0.95; %                     % Maximum setpoint [pu]

% Wind turbine
in.WT.type = WT_type; % 1; %                       % Wind Turbine type: 1 for NREL 5MW, 2 for Hywind 6MW
in.WT.n = WT_n; %5; %                            % Number of Wind Turbines
in.WT.wdata = wdata;                    % Timeseries with wind speeds [m/s]

% Energy storage
in.ES.Pmin = size_electrolyser_bis / par.elGrid.Sbase; %-2 / par.elGrid.Sbase; %     % Minimum setpoint / electrolizer rated power (<= 0) [pu] 
in.ES.Pmax = size_FC / par.elGrid.Sbase; %2 / par.elGrid.Sbase; %      % Maximum sepoint / fuel cell rated power (>= 0) [pu] 
in.ES.Erated = 300;                     % Storage system rated capacity [MWh]
in.ES.Eff.dis = 63;                     % Discharge / fuel-cell efficiency [%]
in.ES.Eff.ch = 63;                      % Charge / electrolyzer efficiency [%]
in.ES.SOC.min = 10;                     % Minimum state-of-charge [%]
in.ES.SOC.max = 90;                     % Maximum state-of-charge [%]
in.ES.SOC.initial = 50;                 % Initial state-of-charge [%]

in.LD.Cont = Cont_load; %(mean(GTLoad) - ldata.On) / par.elGrid.Sbase; %     % Continuous load rated power [pu]
in.LD.TOn.Load = ldata.On / par.elGrid.Sbase; % On transient load rated power [pu]
in.LD.TOn.Time = 70;                    % Time to switch on transient load [s]
in.LD.TOff.Load = -ldata.On / par.elGrid.Sbase; % Off transient load rated power [pu]
in.LD.TOff.Time = 130;                  % Time to switch off transient load [s]
in.LD.Flex = 7.6 / par.elGrid.Sbase;    % Flexible load rated power [pu]
in.LD.Td = 0.1;                         % Flexible load time delay [s]

[out,rv] = elgridEval(in, par);

if rv == 1
    save('test_W24_LM6000','out','-v7.3');
end

frequency_boolean = out.f.OK;
ramp_boolean = out.dPdt.OK;
ramp = out.dPdt.max;

end
