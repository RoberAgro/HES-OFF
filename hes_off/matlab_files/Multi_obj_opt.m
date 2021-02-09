
    A = []; b = [];
    Aeq = []; beq = [];
%     % Wind capacity & Electrolyser stack
%     lb = [12 2];
%     ub = [24 8];
%     % Wind capacity & Electrolyser stack & GT type 
%     lb = [12 2 1];
%     ub = [24 8 2];
    % Wind capacity & Electrolyser stack & GT type & Fuel cell stack & H2 in GT factor
    lb = [12 2 1 1 0];
    ub = [24 8 2 5 3];
%     % Wind capacity & Electrolyser stack & GT type & Fuel cell stack  & H2 in GT factor & Max storage level surplus & Max storage level lack & Max GT load & min GT load 
%     lb = [12 2 1 1 0 3000 3000 0.85 0.40];
%     ub = [24 8 2 5 3 10000 10000 0.95 0.50];   
    
    nvars = length(lb);  %nvars equal to the number of values in the array lb
   
   % options = gaoptimset('populationsize',30,'generations',10,'display','iter','PlotFcn',@gaplotpareto);
    options = gaoptimset('TolFun',1e-3,'populationsize',150,'generations',2,'display','iter','PlotFcn',@gaplotpareto);
    %options = gaoptimset('populationsize',25,'StallGenLimit',2,'TolFun',1e-4,'display','iter','PlotFcns',@gaplotbestf); %'UseParallel',
    %[x,fval] = gamultiobj(@INTEGRATED_MODEL_v2_constantFCandEL_GA,nvars,A,b,Aeq,beq,lb,ub,options);
    %[x,fval] = gamultiobj(@INTEGRATED_MODEL_v2_modelFCandEL_GA,nvars,A,b,Aeq,beq,lb,ub,options);
    [x,fval] = gamultiobj(@INTEGRATED_MODEL_v2_constantFCandEL_storageCONSTRAINT_GA,nvars,A,b,Aeq,beq,lb,ub,options);