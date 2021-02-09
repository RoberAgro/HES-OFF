function Power_Curve = Hywind
% NREL generic wind turbine Power curve 
%
%    Power_curve = Hywind
%

Power_Curve = [0 0;1 0;2 0;3 0;4 0.029;5 0.0725;6 0.1304;7 0.2101;8 0.3261;9 0.4638;10 ...
               0.6232;11 0.7754;12 0.8913;13 0.9565;14 0.9855;15 1;16 1;25 1;26 0;50 0];         
Power_Curve(:,2)=Power_Curve(:,2)./max(Power_Curve(:,2)); % Normalised power curve




