
lb = [100000 100000 100000 100000]; % Lower bound
ub = []; % No Upper bound
A = []; % No linear inequality constraints
b = []; % No linear inequality constraints
Aeq = [1 1 1 1]; % linear equality constraints
beq = [1000000]; % linear equality constraints
FitnessFunction = @subsystemmultiobjective;
numberOfVariables = 4;
[x,fval] = gamultiobj(FitnessFunction,numberOfVariables,A,b,Aeq,beq,lb,ub);

objectivesToPlot = [1,2,3]; % 3D pareto plot or whatever two objectives you want
plotfn = @(options,state,flag)gaplotpareto(options,state,flag,objectivesToPlot);
options = gaoptimset('PlotFcns',plotfn);
gamultiobj(FitnessFunction,numberOfVariables,A,b,Aeq,beq,lb,ub,options);

% X=array2table(x)
% FVAL=array2table(fval)

%plotting pareto solutions table
J = [x fval]
fig = uifigure;
uit = uitable(fig,'Data',J);
uit.ColumnName = {'x1','x2','x3','x4','Objfun1','objfun2','objfun3'};
