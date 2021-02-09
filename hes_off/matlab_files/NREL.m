function Power_Curve = NREL
% NREL generic wind turbine Power curve 
%
%    Power_curve = NREL
%

Power_Curve = [0 0;1 0;2 0;3 170;4 391;5 731;6 1173;7 1752;8 2534;9 3452;10 ...
               4558;11 5000;12 5000;13 5000;14 5000;15 5000;16 5000;25 ...
               5000;30 0;100 0];
Power_Curve(:,2)=Power_Curve(:,2)./max(Power_Curve(:,2)); % Normalised power curve

