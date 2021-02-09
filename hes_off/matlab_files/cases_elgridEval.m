clear 

%% Wind data - from field measurements or a database
wraw = load('sleipnerwind.mat');
Fs = 1/60;              % Sample frequency of the data [Hz]

% Correction of wind speed for height
href = 14;              % Reference height
hhub = 90;              % Hub height
%h0 = 0.0002;            % Roughness lenght - Class 0 (Water surface), according to the European Wind Atlas
%wcorr = fillmissing(wraw.wind(1).w,'linear') * log(hhub/h0) / log(href/h0);
wcorr = fillmissing(wraw.wind(1).w,'linear') * (hhub/href)^(1/7);

% Calculate the derivative of the wind power
wdiff = diff(wcorr);

% Find maximum derivative
[diffMax,idxMax] = max(wdiff);
wdata.max = timeseries(wcorr(idxMax-1:idxMax+2),0:1/Fs:3/Fs);

% Find minimum derivative
[diffMin,idxMin] = min(wdiff);
wdata.min = timeseries(wcorr(idxMin-1:idxMin+2),0:1/Fs:3/Fs);

%% Load data - SCADA data provided by Lundin
GTdata = load('GTLoadMin.mat');
GTLoad = fillmissing(GTdata.GTLoad,'linear');
Fs = 1/60;              % Sample frequency of the data [Hz]

% Calculate the derivative of the load
ldiff = diff(GTLoad);

% Find maximum derivative
[diffMax,idxMax] = max(ldiff);
ldata.max = timeseries([43.6 - diffMax 43.6 43.6 43.6],0:1/Fs:3/Fs);
% ldata.max = timeseries([36.8 40.6 40.6 40.6],0:1/Fs:3/Fs); %% Data from SCADA used for validation

% Find minimum derivative
[diffMin,idxMin] = min(ldiff);
ldata.min = timeseries(GTLoad(idxMin-1:idxMin+2),0:1/Fs:3/Fs);

%% Simulation parameters

% Basic data
par.elGrid.model = 'average_model_WT.slx';  % File name of the model to be evaluated
par.elGrid.Sbase = 43.6;                % Base value apparent power [MVA]    
par.elGrid.fbase = 60;                  % Base value frequency [Hz]
par.elGrid.fmax = 1.02;                 % Maximum allowed frequency in the electrical grid [pu]
par.elGrid.fmin = 0.98;                 % Minimum allowed frequency in the electrical grid [pu]

%% Input data

% Gas turbine
in.GT.type = 2;                         % Gas Turbine type: 1 for LM2500, 2 for LM6000
in.GT.n = 1;                            % Number of Gas Turbines
in.GT.Pmin = 0.4;                       % Minimum setpoint [pu]
in.GT.Pmax = 0.95;                       % Maximum setpoint [pu]

% Wind turbine
in.WT.type = 2;                         % Wind Turbine type: 1 for NREL 5MW, 2 for Hywind 6MW
in.WT.n = 5;                            % Number of Wind Turbines
in.WT.wdata = wdata;                    % Timeseries with wind speeds [m/s]

% Energy storage
in.ES.Pmin = -2 / par.elGrid.Sbase;      % Minimum setpoint / electrolizer rated power (<= 0) [pu] 
in.ES.Pmax = 2 / par.elGrid.Sbase;      % Maximum sepoint / fuel cell rated power (>= 0) [pu] 
in.ES.Erated = 300;                     % Storage system rated capacity [MWh]
in.ES.Eff.dis = 63;                     % Discharge / fuel-cell efficiency [%]
in.ES.Eff.ch = 63;                      % Charge / electrolyzer efficiency [%]
in.ES.SOC.min = 10;                     % Minimum state-of-charge [%]
in.ES.SOC.max = 90;                     % Maximum state-of-charge [%]
in.ES.SOC.initial = 50;                 % Initial state-of-charge [%]

in.LD.ldata = ldata.max / par.elGrid.Sbase;     % Timeseries with load values [pu]

%% Simulation cases

i = 0;

%% LM 2500

% Reference: 2 x LM2500, no WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 1;
simCase(i).GT.n = 2;
simCase(i).GT.Pmin = 0.2;
simCase(i).WT.n = 0;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM2500 + 2 x WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 1;
simCase(i).GT.n = 1;
simCase(i).WT.n = 2;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM2500 + 2 x WT, with ES
i = i + 1;
simCase(i) = simCase(i-1);
simCase(i).ES.Pmin = -4 / par.elGrid.Sbase;
simCase(i).ES.Pmax = 3.6 / par.elGrid.Sbase;

% 1 x LM2500 + 3 x WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 1;
simCase(i).GT.n = 1;
simCase(i).WT.n = 3;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM2500 + 3 x WT, with ES
i = i + 1;
simCase(i) = simCase(i-1);
simCase(i).ES.Pmin = -4 / par.elGrid.Sbase;
simCase(i).ES.Pmax = 3.6 / par.elGrid.Sbase;

% 1 x LM2500 + 4 x WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 1;
simCase(i).GT.n = 1;
simCase(i).WT.n = 4;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM2500 + 4 x WT, with ES
i = i + 1;
simCase(i) = simCase(i-1);
simCase(i).ES.Pmin = -6.1 / par.elGrid.Sbase;
simCase(i).ES.Pmax = 3.6 / par.elGrid.Sbase;


%% LM 6000

% Reference: 1 x LM6000, no WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 2;
simCase(i).GT.n = 1;
simCase(i).GT.Pmin = 0.4;
simCase(i).WT.n = 0;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM6000 + 2 x WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 2;
simCase(i).GT.n = 1;
simCase(i).WT.n = 2;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM6000 + 2 x WT, with ES
i = i + 1;
simCase(i) = simCase(i-1);
simCase(i).ES.Pmin = -4 / par.elGrid.Sbase;
simCase(i).ES.Pmax = 1.1 / par.elGrid.Sbase;

% 1 x LM6000 + 3 x WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 2;
simCase(i).GT.n = 1;
simCase(i).WT.n = 3;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM6000 + 3 x WT, with ES
i = i + 1;
simCase(i) = simCase(i-1);
simCase(i).ES.Pmin = -4 / par.elGrid.Sbase;
simCase(i).ES.Pmax = 1.1 / par.elGrid.Sbase;

% 1 x LM6000 + 4 x WT, no ES
i = i + 1;
simCase(i) = in;
simCase(i).GT.type = 2;
simCase(i).GT.n = 1;
simCase(i).WT.n = 4;
simCase(i).ES.Pmin = 0;
simCase(i).ES.Pmax = 0;

% 1 x LM6000 + 4 x WT, with ES
i = i + 1;
simCase(i) = simCase(i-1);
simCase(i).ES.Pmin = -6.6 / par.elGrid.Sbase;
simCase(i).ES.Pmax = 1.6 / par.elGrid.Sbase;


%% Simulation
for i = 1:length(simCase)
    [out,rv] = elgridEval(simCase(i), par);
    simResults(i) = out;
end

save('test','simResults', 'simCase', '-v7.3');
