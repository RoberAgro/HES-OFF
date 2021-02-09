function Pw = windpower(Power_curve,wspeed,wspeed_factor,Pw_rat)
% WINDPOWER Convert wind to power given speed and power curve
% 
%   Pw = windpower(Power_curver,wspeed,Nw,nshift,wspeed_factor,Pw_rat)

% nargin
% pause

% Simulation parameters and data
if nargin<4
  Pw_rat = 1;          % Rated wind power [MW]
end

if nargin<3||isempty(wspeed_factor)
  wspeed_factor = 1.25; %1; % Multiplication factor for wind speed
end

%% Create wind power output for each time step
Pw = interp1(Power_curve(:,1),Power_curve(:,2),wspeed_factor*wspeed).*Pw_rat;
Pw(Pw<0)=0;

return
