%% Minimising Difference of Key Nutrients using fmincon
obj1 = @(x) lsq(x); %call f1a

T = readmatrix('Matlab_data.xlsx'); %read excel table
DF1 = T(1,2:9); %create array of defficiency factors for first 8 nutrients
D1 = T(2,2:9); %create array of dairy values for first 8 nutrients 
NF = T(3:6,2:21); %create matrix of nutrient factors

x0 = [1,1,1,1]; %starting values

A =[1,1,1,1;
    -1,0,0,0;
    0,-1,0,0;
    0,0,-1,0;
    0,0,0,-1;
    1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    0.0131,0.033,0.0121,0.045;
    0.0304,0.015,0.005,0.043;
    0.0000101,0.00005,0.00002,0;
    1.73,5.53,3.89,6.28
    ]; %linear constraints

b=[160;
   0;
   0;
   0;
   0;
   128;
   47;
   143;
   46;
   5.467; 
   53;
   0.136;
   400
   ]; %linear constrain values

Aeq = [];
beq = []; %no equality constraints

lb = []; %lower bound of soy, hemp, oat, hazelnut
ub = []; %upper bound of soy, hemp, oat, hazelnut

nonlcon = []; %no nonlinear constraints

options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); %find active set of solutions
[y,gval] = fmincon(obj1,x0,A,b,Aeq,beq,lb,ub,nonlcon,options) %minimise sum of squares of difference between dairy and sum of plant crops

Nut1 = y'.*NF; %calculate nutrient values of each crop using function solution
%Profile1 = total nutrient values for final product solution 
Profile1 = [sum(Nut1(1:4,1)),sum(Nut1(1:4,2)),sum(Nut1(1:4,3)),sum(Nut1(1:4,4)),sum(Nut1(1:4,5)),sum(Nut1(1:4,6)),sum(Nut1(1:4,7)),sum(Nut1(1:4,8)),sum(Nut1(1:4,9)),sum(Nut1(1:4,10)),sum(Nut1(1:4,11)),sum(Nut1(1:4,12)),sum(Nut1(1:4,13)),sum(Nut1(1:4,14)),sum(Nut1(1:4,15)),sum(Nut1(1:4,16)),sum(Nut1(1:4,17)),sum(Nut1(1:4,18)),sum(Nut1(1:4,19)),sum(Nut1(1:4,20))]

X = [0:25:600]; %X = calorie constraint values from 0 to 600 (dairy milk value)
%Y = fmincon output values (gval) for different calorie constraint
%function aims to reach 0 (i.e. 0 difference between dairy and plant milk nutrient values
Y = [191,163,138,114.8,93.7,74.8,58,43.4,30.9,21,15.8,11.31,7.6,4.64,2.4,0.99,0.3,0.23,0.23,0.23,0.23,0.25,0.25,0.24,0.24];

figure %create figure for calorie parametric study
scatter(X,Y) %scatter graph of calorie constrain limit vs f1a function output
title('Effect of Calorie Constraint (g13) on f1a')
xlabel('Calories (Kcal)') %axis labels 
ylabel('f1')

%% Minimising Difference of Key Nutrients Individually using fminimax
obj4 = @(x) goal(x); %call f1b

T = readmatrix('Matlab_data.xlsx'); %read excel table
T(1:6,1:20); 
DF1 = T(1,2:9); %create array of defficiency factors for first 8 nutrients
D1 = T(2,2:9); %create array of dairy values for first 8 nutrients
NF = T(3:6,2:21); %create matrix of nutrient factors 

x0 = [1,1,1,1]; %create starting values

A = [1,1,1,1;
    -1,0,0,0;
    0,-1,0,0;
    0,0,-1,0;
    0,0,0,-1;
    1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    0.0131,0.033,0.0121,0.045;
    0.0304,0.015,0.005,0.043;
    0.0000101,0.00005,0.00002,0;
    1.73,5.53,3.89,6.28]; %create constraints

b =[160;
    0;
    0;
    0;
    0;
    128;
    47;
    143;
    46;
    5.467; 
    53;
    0.136;
    400]; %create constraints

Aeq = []; %no equality constraints
beq = [];

lb = [0,0,0,0]; %lower bound of grams of each crop
ub = []; %upper bound of soy, hemp, oat, hazelnut

nonlcon = []; %no nonlinear constraints

   
[x,fval] = fminimax(obj4,x0,A,b,Aeq,beq,lb,ub) %minimise the difference between dairy and solution for each individual nutrient 

Nut2 = x'.*NF; %calculate values of each nutrient for each crop in function solution
%calculate total nutrient values for final solution of f1b
Profile2 = [sum(Nut2(1:4,1)),sum(Nut2(1:4,2)),sum(Nut2(1:4,3)),sum(Nut2(1:4,4)),sum(Nut2(1:4,5)),sum(Nut2(1:4,6)),sum(Nut2(1:4,7)),sum(Nut2(1:4,8)),sum(Nut2(1:4,9)),sum(Nut2(1:4,10)),sum(Nut2(1:4,11)),sum(Nut2(1:4,12)),sum(Nut2(1:4,13)),sum(Nut2(1:4,14)),sum(Nut2(1:4,15)),sum(Nut2(1:4,16)),sum(Nut2(1:4,17)),sum(Nut2(1:4,18)),sum(Nut2(1:4,19)),sum(Nut2(1:4,20))]


%% Maximising Healthy Nutrients using linprog
T = readmatrix('Matlab_data.xlsx'); %read table from excel
T(1:6,1:20);
NF = T(3:6,2:21); %create matrix of nutrient factors 
DF2 = T(1,10:17); %create array of defficiency factors for nutrients 10 to 17 
NF2 = [T(3,10:17); %create matrix of nutrient factors for relevant healthy nutrients 
    T(4,10:17);
    T(5,10:17);
    T(6,10:17)]; 

NFX2 = DF2.*NF2; %multiply nutrient proportions by deficiency factor 
A2 = [sum(NFX2(1,1:8)),sum(NFX2(2,1:8)),sum(NFX2(3,1:8)),sum(NFX2(4,1:8))]; %sum coefficients for x1 (soy), x2 (hemp), x3 (oat), x4 (hazelnut)
obj3 = -A2; %create f1b coefficients (in negative null form)

A =[1,1,1,1;
    -1,0,0,0;
    0,-1,0,0;
    0,0,-1,0;
    0,0,0,-1;
    1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    0.0131,0.033,0.0121,0.045;
    0.0304,0.015,0.005,0.043;
    0.0000101,0.00005,0.00002,0;
    1.73,5.53,3.89,6.28
    ]; %inequality constraint coefficients

b =[160;
    0;
    0;
    0;
    0;
    128;
    47;
    143;
    46;
    5.467; 
    53;
    0.136;
    400
    ]; %inequality constrain bounds

Aeq = [];
beq = []; %no equality constraints

lb = []; %lower bound of soy, hemp, oat, hazelnut
ub = []; %upper bound of soy, hemp, oat, hazelnut

[w,eval] = linprog(obj3,A,b,Aeq,beq,lb,ub) %use linear programming optimisation to solve f2

Nut3 = w.*NF; %calculate values of each nutrient for each crop in function solution
%calculate total nutrient values for final solution of f2
Profile3 = [sum(Nut3(1:4,1)), sum(Nut3(1:4,2)), sum(Nut3(1:4,3)), sum(Nut3(1:4,4)), sum(Nut3(1:4,5)), sum(Nut3(1:4,6)), sum(Nut3(1:4,7)), sum(Nut3(1:4,8)), sum(Nut3(1:4,9)), sum(Nut3(1:4,10)), sum(Nut3(1:4,11)), sum(Nut3(1:4,12)), sum(Nut3(1:4,13)), sum(Nut3(1:4,14)), sum(Nut3(1:4,15)), sum(Nut3(1:4,16)), sum(Nut3(1:4,17)), sum(Nut3(1:4,18)), sum(Nut3(1:4,19)), sum(Nut3(1:4,20))]

compval = [Profile1(1:19);Profile2(1:19);Profile3(1:19)]; %compare nutrient profile of f1a, f1b, and f2 solutions 
%create nutrient labels
compvar = categorical({'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium'});
compvar = reordercats(compvar,{'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium'});

bar(compvar,compval) %bar chart of f1a, f1b, f2 nutrient profiles (calories not included)
title('f1a, f1b, f2 Solution Nutrient Profiles')
ylabel('Nutrient Value (g)')
lgd = legend('f1a','f1b','f2');

%% Multiobjective - obtaining multiobjective solutions to minimise difference and maximise beneficial nutrients 

fun = @(x) big_obj(x); %call multi-objective function (f1a and f2)

T = readmatrix('Matlab_data.xlsx'); %read excel table
T(1:6,1:20);
NF = T(3:6,2:21); %create matrix of nutrient factors 
D = [T(2,2:21)]; %nutrient values for dairy

x0 = [1,1,1,1]; %starting values 

A = [1,1,1,1;
    -1,0,0,0;
    0,-1,0,0;
    0,0,-1,0;
    0,0,0,-1;
    1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    0.0131,0.033,0.0121,0.045;
    0.0304,0.015,0.005,0.043;
    0.0000101,0.00005,0.00002,0;
    1.73,5.53,3.89,6.28]; %inequality constraint coefficients

b=[160;
     0;
    0;
    0;
    0;
    128;
    47;
    143;
    46;
    5.467; 
    53;
    0.136;
    400]; %inequality constraint bounds

Aeq = [];
beq = []; %no equality constraints

lb = [0,0,0,0]; %lower bound of soy, hemp, oat, hazelnut
ub = []; %upper bound of soy, hemp, oat, hazelnut

nonlcon = []; %no nonlinear constraints
options = optimoptions('gamultiobj','MaxGenerations',2,'PlotFcn','gaplotpareto'); %plot pareto front as algorithm runs, limit generations to prevent long run time 
% find pareto front and show generation of possible solutions as algorithm runs 
[t,bval,exitflag,output,population,scores] = gamultiobj(fun,4,A,b,Aeq,beq,lb,ub,options); %use multiobjective algorithm 

[numRows,numCols] = size(t); %size of array containing solutions
I = numRows; %create variable 'I' which counts number of solutions

for i = 1:I %create for loop to extract solutions individually from array
    newt = strcat('x',num2str(i)); %create structure that contains each solution and names each solution using string
    sol.(newt) = t(i,1:4); %assign to solutions values 't'
    sol.(newt) = (sol.(newt)).*NF' %multiply solutions by nutrient proportions to get actual amount of each nutrient in solutions
    sol.(newt) = [sum(sol.(newt)(1,1:4)),sum(sol.(newt)(2,1:4)),sum(sol.(newt)(3,1:4)),sum(sol.(newt)(4,1:4)),sum(sol.(newt)(5,1:4)),sum(sol.(newt)(6,1:4)),sum(sol.(newt)(7,1:4)),sum(sol.(newt)(8,1:4)),sum(sol.(newt)(9,1:4)),sum(sol.(newt)(10,1:4)),sum(sol.(newt)(11,1:4)),sum(sol.(newt)(12,1:4)),sum(sol.(newt)(13,1:4)),sum(sol.(newt)(14,1:4)),sum(sol.(newt)(15,1:4)),sum(sol.(newt)(16,1:4)),sum(sol.(newt)(17,1:4)),sum(sol.(newt)(18,1:4)),sum(sol.(newt)(19,1:4)),sum(sol.(newt)(20,1:4))];
    sol.(newt) = sol.(newt)%sum the amount of each nutrient from each crop to create a total nutrient value for that solution
   
end

Y = D; %create array Y containing dairy nutrient profile 

for j = [1:length(fieldnames(sol))] %loop to add each solution to comparison array Y
    test = sol.(strcat('x',num2str(j)))
    Y(j+1,:) = test; %index solutions 
end

%solution visualisation

%create nutrient names for graphs
X = categorical({'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium','calories'});
X = reordercats(X,{'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium','calories'});
J = normalize(Y); %normalise solution and dairy nutrient values for initial comparison 
figure(); 
bar(X,J) %create bar graph of all solutions and dairy nutrient values
title('Solutions vs Dairy')
xlabel('Nutrients')
ylabel('Normalised Nutrient Score')

Ymacro = [Y(1:19,1),Y(1:19,9),Y(1:19,8),Y(1:19,17),Y(1:19,18)]; %matrix of macro nutrient solution values only
Xmacro = categorical({'protein','fibre','Unsaturated fat','saturated fat','sugar'}); %create nutrient names for graphs
Xmacro = reordercats(Xmacro,{'protein','fibre','Unsaturated fat','saturated fat','sugar'}) 
lgd = legend('Dairy');
figure();
bar(Xmacro,Ymacro) %bar graph of macronutrient value comparison
title('Solutions vs Dairy')
xlabel('Macronutrients')
ylabel('Nutrient Value (g)')
lgd = legend('Dairy');

%create matrix of micronutrient values only 
Ymicro = [Y(1:19,2),Y(1:19,3),Y(1:19,4),Y(1:19,5),Y(1:19,6),Y(1:19,7),Y(1:19,10),Y(1:19,11),Y(1:19,12),Y(1:19,13),Y(1:19,14),Y(1:19,15),Y(1:19,16),Y(1:19,19)] %solutions micronutrient values
Xmicro = categorical({'calcium','potassium','vitamin a','vitamin e','vitamin k','B2','iron','magnesium','vitamin c','zinc','B1','B3','B6','sodium'});%create nutrient names for graphs
Xmicro = reordercats(Xmicro,{'calcium','potassium','vitamin a','vitamin e','vitamin k','B2','iron','magnesium','vitamin c','zinc','B1','B3','B6','sodium'})
figure();
bar(Xmicro,Ymicro) %bar graph of micronutrient value comparison
title('Solutions vs Dairy')
xlabel('Micronutrients')
ylabel('Nutrient Value (g)')
lgd = legend('Dairy');

Ycal = [Y(1:19,20)]; %solution and dairy calorie values
Xcal = categorical({'calories'}); %create nutrient names for graphs
figure();
bar(Xcal,Ycal) %bar graph of calorie value comparison
title('Solutions vs Dairy')
xlabel('-')
ylabel('Calories (Kcal)')
lgd = legend('Dairy');

%Comparison of individual random solution from pareto front vs dairy

Ymaccomp = [Y(1,1),Y(1,9),Y(1,8),Y(1,17),Y(1,18);Y(10,1),Y(10,9),Y(10,8),Y(10,17),Y(10,18)]; %random solution macro nutrient values
Xmaccomp = categorical({'protein','fibre','Unsaturated fat','saturated fat','sugar'}); %create nutrient names for graphs
Xmaccomp = reordercats(Xmaccomp,{'protein','fibre','Unsaturated fat','saturated fat','sugar'}) 
figure();
macbar = bar(Xmaccomp,Ymaccomp) %bar graph of macronutrient value comparison of possible solution
title('Possible Solution vs Dairy')
xlabel('Macronutrients')
ylabel('Nutrient Value (g)')
lgd = legend('Dairy','Possible Solution');

%random solution micronutrient comparison
Ymiccomp = [Y(1,2),Y(1,3),Y(1,4),Y(1,5),Y(1,6),Y(1,7),Y(1,10),Y(1,11),Y(1,12),Y(1,13),Y(1,14),Y(1,15),Y(1,16),Y(1,19);Y(10,2),Y(10,3),Y(10,4),Y(10,5),Y(10,6),Y(10,7),Y(10,10),Y(10,11),Y(10,12),Y(10,13),Y(10,14),Y(10,15),Y(10,16),Y(10,19)]; %solutions macro nutrient values
Xmiccomp = categorical({'calcium','potassium','vitamin a','vitamin e','vitamin k','B2','iron','magnesium','vitamin c','zinc','B1','B3','B6','sodium'});%create nutrient names for graphs
Xmiccomp = reordercats(Xmiccomp,{'calcium','potassium','vitamin a','vitamin e','vitamin k','B2','iron','magnesium','vitamin c','zinc','B1','B3','B6','sodium'})
figure();
micbar = bar(Xmiccomp,Ymiccomp) %bar graph of micronutrient value comparison of possible solution
title('Possible Solution vs Dairy')
xlabel('Micronutrients')
ylabel('Nutrient Value (g)')
lgd = legend('Dairy','Possible Solution');

%Random solution calorie comparison with dairy 
Ycalcomp = [Y(1,20);Y(10,20)]; %solution and dairy calorie values
Xcalcomp = categorical({'calories'}); %create nutrient names for graphs
figure();
calbar = bar(Xcalcomp,Ycalcomp)
title('Possible Solution vs Dairy')
ylabel('Calories (Kcal)')
lgd = legend('Dairy','Possible Solution');


%% Final Solution Optimisation. 
%Using Healthy Nutrients as a Constraint, Minimising Difference with fminimax

obj4 = @(x) goal(x); %call f1b
maxcon = @(x) max(x); %call f2 as a function 

T = readmatrix('Matlab_data.xlsx'); %read excel table
T(1:6,1:20); 
DF1 = T(1,2:9); %create array of defficiency factors for first 8 nutrients
D1 = T(2,2:9); %create array of dairy values for first 8 nutrients
NF = T(3:6,2:21); %create matrix of nutrient factors 

x0 = [1,1,1,1]; %create starting values

A = [1,1,1,1;
    -1,0,0,0;
    0,-1,0,0;
    0,0,-1,0;
    0,0,0,-1;
    1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    0.0131,0.033,0.0121,0.045;
    0.0304,0.015,0.005,0.043;
    0.0000101,0.00005,0.00002,0;
    1.73,5.53,3.89,6.28;
    ]; %create constraints

b =[160;
    0;
    0;
    0;
    0;
    128;
    47;
    143;
    46;
    5.467; 
    53;
    0.136;
    400]; %create constraints

Aeq = []; %no equality constraints
beq = [];

lb = [0,0,0,0]; %lower bound of grams of each crop
ub = []; %upper bound of soy, hemp, oat, hazelnut

nonlcon = [maxcon]; %set non linear constraints as f2 (maximise healthy nutrients function)

  
[fin,finval] = fminimax(obj4,x0,A,b,Aeq,beq,lb,ub) %minimise the difference of each nutrient for dairy and the sum of each plant crop

Nut4 = fin'.*NF; %calculate amounts of each nutrient in each crop in final solution
%calculate final solution nutrient values 
Profile4 = [sum(Nut4(1:4,1)),sum(Nut4(1:4,2)),sum(Nut4(1:4,3)),sum(Nut4(1:4,4)),sum(Nut4(1:4,5)),sum(Nut4(1:4,6)),sum(Nut4(1:4,7)),sum(Nut4(1:4,8)),sum(Nut4(1:4,9)),sum(Nut4(1:4,10)),sum(Nut4(1:4,11)),sum(Nut4(1:4,12)),sum(Nut4(1:4,13)),sum(Nut4(1:4,14)),sum(Nut4(1:4,15)),sum(Nut4(1:4,16)),sum(Nut4(1:4,17)),sum(Nut4(1:4,18)),sum(Nut4(1:4,19)),sum(Nut4(1:4,20))]


FinComp = [D(1:19);Profile4(1:19)];
compvar = categorical({'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium'});
compvar = reordercats(compvar,{'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium'});
figure
bar(compvar,FinComp) %bar chart of final solution and dairy nutrient profiles (calories not included)
title('Final Solution Nutrient Profile vs Dairy')
ylabel('Nutrient Value (g)')
lgd = legend('Dairy','Final Solution');

Table = [D;Profile4]; %create table of final solution vs dairy 
Table_of_Solution = array2table(Table, 'VariableNames',{'protein','calcium','potassium','vitamin a','vitamin e','vitamin k','B2','Unsaturated fat','fibre','iron','magnesium','vitamin c','zinc','B1','B3','B6','saturated fat','sugar','sodium','calories'}, 'RowNames', {'Dairy','Final Solution'})