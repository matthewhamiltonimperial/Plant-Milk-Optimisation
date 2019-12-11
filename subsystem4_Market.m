clear
clc

%% Read raw data from .csv file

Fat = csvread('Fat_Protein_Type.csv',0,0,[0,0,3,0]);                        % [0g/8oz,2.5g/8oz,5g/8oz,8g/8oz]
fat_utility = csvread('Fat_Protein_Type.csv',0,1,[0,1,3,1]);                % Utility score of fat contents
Fat_Utility = normalize(fat_utility,'range');                               % Normalised so that 0 < Utility < 1
Protein = csvread('Fat_Protein_Type.csv',0,2,[0,2,4,2]);                    % [1g/8oz,2g/8oz,5g/8oz,8g/8oz,9g/oz]
protein_utility = csvread('Fat_Protein_Type.csv',0,3,[0,3,4,3]);            % Utility score of protein contents
Protein_Utility = normalize(protein_utility,'range');                       % Normalised so that 0 < Utility < 1
type_utility = csvread('Fat_Protein_Type.csv',0,4,[0,4,3,4]);               % [Hazelnut, Oat, Soy, Hemp]
Type_Utility = normalize(type_utility,'range');                             % Normalised so that 0 < Utility < 1
Fat_Content = csvread('Fat_Protein_Type.csv',0,5,[0,5,3,5]);                % [Hazelnut fat, Oat fat, Soy fat, Hemp fat]
Protein_Content = csvread('Fat_Protein_Type.csv',0,6,[0,6,3,6]);            % [Hazelnut protein, Oat protein, Soy protein, Hemp protein]
input_costs = csvread('Optimisation_Main_2.0.csv',0,2,[0,2,0,8]);           % [Hazelnut,Oat,Soy,Hemp,Packaging,Production,Margin]
quantity_sold = csvread('Optimisation_Main_2.0.csv',0,1,[0,1,8,1]);         % Market data of plant milk - quantity sold
Quantity_Sold = normalize(quantity_sold,'range');                           % Normalised so that 0 < Quantity Sold < 1
Price_Per_Litre = csvread('Optimisation_Main_2.0.csv',0,0,[0,0,8,0]);       % Market data of plant milk - price
coeff_fat = polyfit(Fat,Fat_Utility,2);                                     % Fitting a curve of fat content vs utility and returning coefficients
coeff_protein = polyfit(Protein,Protein_Utility,2);                         % Fitting a curve of protein content vs utility and returning coefficients
coeff_demand = polyfit(Price_Per_Litre,Quantity_Sold,2);                    % Fitting a curve of type vs utility and returning coefficients

global c1           % Cost of 10g (1%) of Hazelnuts
global c2           % Cost of 10g (1%) of Oats
global c3           % Cost of 10g (1%) of Soy
global c4           % Cost of 10g (1%) of Hemp
global c5           % Cost to process 10g (1%) of plant
global c6           % Packaging costs per litre plant milk
global k1           % Quadratic coefficient of demand curve
global k2           % Linear coefficient of demand curve
global k3           % Constant coefficient of demand curve
global k4           % Quadratic coefficient of fat curve
global k5           % Linear coefficient of fat curve
global k6           % Constant coefficient of fat curve
global k7           % Quadratic coefficient of protein curve
global k8           % Linear coefficient of protein curve
global k9           % Constant coefficient of protein curve
global k10          % Utility of Hazelnut
global k11          % Utility of Oat
global k12          % Utility of Soy
global k13          % Utility of Hemp
global p1           % Protein content in 10g of Hazelnut
global p2           % Protein content in 10g of Oat
global p3           % Protein content in 10g of Soy
global p4           % Protein content in 10g of Hemp
global f1           % Fat content in 10g of Hazelnut
global f2           % Fat content in 10g of Oat
global f3           % Fat content in 10g of Soy
global f4           % Fat content in 10g of Hemp

c6 = input_costs(1,5);          % ASSIGNING GLOBAL VARIABLES (DEFINED ABOVE)                                                      
c5 = input_costs(1,6);                                                      
c1 = input_costs(1,1);                                                       
c2 = input_costs(1,2);                                                      
c3 = input_costs(1,3);                                                      
c4 = input_costs(1,4);          
k1 = coeff_demand(1,1);
k2 = coeff_demand(1,2);
k3 = coeff_demand(1,3);
k4 = coeff_fat(1);
k5 = coeff_fat(2);
k6 = coeff_fat(3);
k7 = coeff_protein(1);
k8 = coeff_protein(2);
k9 = coeff_protein(3);
k10 = Type_Utility(1);
k11 = Type_Utility(2);
k12 = Type_Utility(3);
k13 = Type_Utility(4);
p1 = Protein_Content(1);
p2 = Protein_Content(2);
p3 = Protein_Content(3);
p4 = Protein_Content(4);
f1 = Fat_Content(1);
f2 = Fat_Content(2);
f3 = Fat_Content(3);
f4 = Fat_Content(4);

%% Visualisation 

figure('Name','Relationships between protein, fat and type on Utility Coefficient')
subplot(2,2,1)
plot(Fat,Fat_Utility, 'r.');                                % Plotting survey data points, fat vs utility
hold on
X_fat = linspace(0,40);
Y_fat = polyval(coeff_fat,X_fat);
plot(X_fat,Y_fat)                                           % Plotting 2nd order regression curve of survey data
hold off
ylim([-0.1 1.1])
title('Effect of Fat content on Utility')
xlabel('Fat content (grams/litre)')
ylabel('Utility coefficient')
subplot(2,2,2)
plot(Protein,Protein_Utility, 'r.');                        % Plotting survey data points, protein vs utility
hold on
X_protein = linspace(0,40);
Y_protein = polyval(coeff_protein,X_protein);
plot(X_protein,Y_protein)                                   % Plotting 2nd order regression curve of survey data
hold off
ylim([-0.1 1.1])
title('Effect of Protein content on Utility')
xlabel('Protein content (grams/litre)')
ylabel('Utility coefficient')
subplot(2,2,[3 4])
bar(type_utility)                                           % Plotting survey data points of type utility
set(gca,'xticklabel',{'Hazelnut','Oat','Soy','Hemp'})       % (Not using normalised data for the purpose of plotting)
title('Effect of Type on Utility')
ylabel('Utility coefficient')


figure('Name','Plant Milk demand schedule')
plot(Price_Per_Litre,Quantity_Sold, 'b.');                  % Plotting market data, price per litre vs volume plant milk sold
hold on
X_demand = linspace(0,5);
Y_demand = polyval(coeff_demand,X_demand);
plot(X_demand,Y_demand, 'r')                                % Plotting 2nd order regression curve of market data                                
hold off
ylim([-0.1 1.1])
title('Plant Milk Demand Schedule')
xlabel('Price ($/litre)')
ylabel('Normalised Quantity Sold')


figure('Name','Surface Utility Plots')
colormap winter                                             % Changing colormap so that higher values are Green not red
subplot(2,3,1)                                              % Hazelnuts vs Oat 
[X,Y] = meshgrid(0:.5:15);                                  % Generating 41x41 meshgrid from 0% to 20%
Z = (1/3).*((((abs(k4.*(X.*f1 + Y.*f2).^2 + k5.*(X.*f1 + Y.*f2) + k6)+(k4.*(X.*f1 + Y.*f2).^2 + k5.*(X.*f1 + Y.*f2) + k6))./2))...
     + (((abs(k7.*(X.*p1+Y.*p2).^2 + k8.*(X.*p1+Y.*p2) + k9)+(k7.*(X.*p1+Y.*p2).^2 + k8.*(X.*p1+Y.*p2) + k9))./2))...
     + ((X./(X+Y)).*k10 + (Y./(X+Y)).*k11));                % Calculating utility given X% Haz and Y% Oat (with 0% Soy/Hemp)
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Oat % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+12;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hazelnut % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+3;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Hazelnut vs Oat')
xlabel('% Hazelnut')
ylabel('% Oat')
zlabel('Utility')
subplot(2,3,2)                                              % Hazelnuts vs Soy
[X,Y] = meshgrid(0:.5:15);                                  % Generating 41x41 meshgrid from 0% to 20%
Z = (1/3).*((((abs(k4.*(X.*f1 + Y.*f3).^2 + k5.*(X.*f1 + Y.*f3) + k6)+(k4.*(X.*f1 + Y.*f3).^2 + k5.*(X.*f1 + Y.*f3) + k6))./2))...
     + (((abs(k7.*(X.*p1+Y.*p3).^2 + k8.*(X.*p1+Y.*p3) + k9)+(k7.*(X.*p1+Y.*p3).^2 + k8.*(X.*p1+Y.*p3) + k9))./2))...
     + ((X./(X+Y)).*k10 + (Y./(X+Y)).*k12));                % Calculating utility given X% Haz and Y% Soy (with 0% Oat/Hemp)
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Soy % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+10;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hazelnut % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+3;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Hazelnut vs Soy')
xlabel('% Hazelnut')
ylabel('% Soy')
zlabel('Utility')
subplot(2,3,3)                                              % Hazelnuts vs Hemp
[X,Y] = meshgrid(0:.5:15);                                  % Generating 41x41 meshgrid from 0% to 20%
Z = (1/3).*((((abs(k4.*(X.*f1 + Y.*f4).^2 + k5.*(X.*f1 + Y.*f4) + k6)+(k4.*(X.*f1 + Y.*f4).^2 + k5.*(X.*f1 + Y.*f4) + k6))./2))...
     + (((abs(k7.*(X.*p1+Y.*p4).^2 + k8.*(X.*p1+Y.*p4) + k9)+(k7.*(X.*p1+Y.*p4).^2 + k8.*(X.*p1+Y.*p4) + k9))./2))...
     + ((X./(X+Y)).*k10 + (Y./(X+Y)).*k13));                % Calculating utility given X% Haz and Y% Hemp (with 0% Soy/Oat)
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hemp % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+8;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hazelnut % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+3;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Hazelnut vs Hemp')
xlabel('% Hazelnut')
ylabel('% Hemp')
zlabel('Utility')
subplot(2,3,4)                                              % Oat vs Soy
[X,Y] = meshgrid(0:.5:15);                                  % Generating 41x41 meshgrid from 0% to 20%
Z = (1/3).*((((abs(k4.*(X.*f2 + Y.*f3).^2 + k5.*(X.*f2 + Y.*f3) + k6)+(k4.*(X.*f2 + Y.*f3).^2 + k5.*(X.*f2 + Y.*f3) + k6))./2))...
     + (((abs(k7.*(X.*p2+Y.*p3).^2 + k8.*(X.*p2+Y.*p3) + k9)+(k7.*(X.*p2+Y.*p3).^2 + k8.*(X.*p2+Y.*p3) + k9))./2))...
     + ((X./(X+Y)).*k11 + (Y./(X+Y)).*k12));                % Calculating utility given X% Oat and Y% Soy (with 0% Haz/Hemp)
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Soy % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+10;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Oat % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+12;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Oat vs Soy')
xlabel('% Oat')
ylabel('% Soy')
zlabel('Utility')
subplot(2,3,5)                                              % Oat vs Hemp
[X,Y] = meshgrid(0:.5:15);                                  % Generating 41x41 meshgrid from 0% to 20%
Z = (1/3).*((((abs(k4.*(X.*f2 + Y.*f4).^2 + k5.*(X.*f2 + Y.*f4) + k6)+(k4.*(X.*f2 + Y.*f4).^2 + k5.*(X.*f2 + Y.*f4) + k6))./2))...
     + (((abs(k7.*(X.*p2+Y.*p4).^2 + k8.*(X.*p2+Y.*p4) + k9)+(k7.*(X.*p2+Y.*p4).^2 + k8.*(X.*p2+Y.*p4) + k9))./2))...
     + ((X./(X+Y)).*k11 + (Y./(X+Y)).*k13));                % Calculating utility given X% Oat and Y% Hemp (with 0% Haz/Soy)
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hemp % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+8;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Oat % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+12;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Oat vs Hemp')
xlabel('% Oat')
ylabel('% Hemp')
zlabel('Utility')
subplot(2,3,6)                                              % Soy vs Hemp
[X,Y] = meshgrid(0:.5:15);                                  % Generating 41x41 meshgrid from 0% to 20%
Z = (1/3).*((((abs(k4.*(X.*f3 + Y.*f4).^2 + k5.*(X.*f3 + Y.*f4) + k6)+(k4.*(X.*f3 + Y.*f4).^2 + k5.*(X.*f3 + Y.*f4) + k6))./2))...
     + (((abs(k7.*(X.*p3+Y.*p4).^2 + k8.*(X.*p3+Y.*p4) + k9)+(k7.*(X.*p3+Y.*p4).^2 + k8.*(X.*p3+Y.*p4) + k9))./2))...
     + ((X./(X+Y)).*k12 + (Y./(X+Y)).*k13));                % Calculating utility given X% Soy and Y% Hemp (with 0% Hemp/Haz)
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hemp % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+8;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Soy % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+10;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Soy vs Hemp')
xlabel('% Soy')
ylabel('% Hemp')
zlabel('Utility')


figure('Name','Surface Demand Plots')
colormap hsv                                                % Changing colormap for more contrast
subplot(2,3,1)                                              % Hazelnuts vs Oat 
[X,Y] = meshgrid(0:.5:15);                                  %Calculating price - and therefore demand (Z) - given X% Haz and Y% Oat
Z = k1*((X.*(c1+c5)+Y.*(c2+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c1+c5)+Y.*(c2+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Oat % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+12;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hazelnut % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+3;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Hazelnut vs Oat')
xlabel('% Hazelnut')
ylabel('% Oat')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Demanded')
subplot(2,3,2)                                              % Hazelnuts vs Soy 
[X,Y] = meshgrid(0:.5:15);                                  %Calculating price - and therefore demand (Z) - given X% Haz and Y% Soy
Z = k1*((X.*(c1+c5)+Y.*(c3+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c1+c5)+Y.*(c3+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Soy % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+10;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hazelnut % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+3;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Hazelnut vs Soy')
xlabel('% Hazelnut')
ylabel('% Soy')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Demanded')
subplot(2,3,3)                                              % Hazelnuts vs Hemp 
[X,Y] = meshgrid(0:.5:15);                                  %Calculating price - and therefore demand (Z) - given X% Haz and Y% Hemp
Z = k1*((X.*(c1+c5)+Y.*(c4+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c1+c5)+Y.*(c4+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hemp % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+8;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hazelnut % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+3;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Hazelnut vs Hemp')
xlabel('% Hazelnut')
ylabel('% Hemp')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Demanded')
subplot(2,3,4)                                              % Oat vs Soy 
[X,Y] = meshgrid(0:.5:15);                                  %Calculating price - and therefore demand (Z) - given X% Oat and Y% Soy
Z = k1*((X.*(c2+c5)+Y.*(c3+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c2+c5)+Y.*(c3+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Soy % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+10;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Oat % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+12;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Oat vs Soy')
xlabel('% Oat')
ylabel('% Soy')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Demanded')
subplot(2,3,5)                                              % Oat vs Hemp 
[X,Y] = meshgrid(0:.5:15);                                  %Calculating price - and therefore demand (Z) - given X% Oat and Y% Hemp
Z = k1*((X.*(c2+c5)+Y.*(c4+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c2+c5)+Y.*(c4+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hemp % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+8;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Oat % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+12;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Oat vs Hemp')
xlabel('% Oat')
ylabel('% Hemp')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Demanded')
subplot(2,3,6)                                              % Soy vs Hemp 
[X,Y] = meshgrid(0:.5:15);                                  %Calculating price - and therefore demand (Z) - given X% Soy and Y% Hemp
Z = k1*((X.*(c3+c5)+Y.*(c4+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c3+c5)+Y.*(c4+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);                                            % Surface plot
hold on
[Xx,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Hemp % Constraint
Yy = (Xx-Xx)+(Zz-Zz)+8;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
[Yy,Zz] = meshgrid(0:15:15,0:1:1);                          % Plot plane at Soy % Constraint
Xx = (Yy-Yy)+(Zz-Zz)+10;
S = surf(Xx,Yy,Zz);
alpha(S,.1)                                                 % Make plane translucent 
hold off
title('Soy vs Hemp')
xlabel('% Soy')
ylabel('% Hemp')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Demanded')

%% SQP Initial Implimentation and Multi-Start 

options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1000);    % Using SQP algorithm 
A = [1 1 1 1];
b = [100];
Aeq = [];
beq = [];
x0 = [1 1 1 1];
lb = [0 0 0 0];                                                             % Lower bounds
ub = [3 12 10 8];                                                           % Upper bounds
nonlcon = [];
[xopt,fval] = fmincon(@objective_initial,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
fprintf('An initial formulation using SQP yielded a minima of %.3f with optimisers: %.3f, %.3f, %.3f, and %.3f',fval,xopt(1),xopt(2),xopt(3),xopt(4)) 

rng default;                                                                % For reproducibility
ms = MultiStart;
problem = createOptimProblem('fmincon','x0',[0 0 0 0],'objective',@objective_initial,'lb',[0 0 0 0],'ub',[3 12 10 8]);
[X,fval,exitflag,output,manymins] = run(ms,problem,50);
fprintf('Multi-start converged from 50 locations to %.3f with optimisers = %.3f, %.3f, %.3f, and %.3f',fval,X(1),X(2),X(3),X(4)) 


%% Parametric Analysis

x_p5 = linspace(0,0.1);
x_p6 = linspace(0,1);
x_p7 = linspace(0,0.8);
n = 0;
y_p5 = [];
y_p6 = [];
y_p7 = [];

for i = 0:0.001:0.1
    n = n+1;
    y = objective_parametric([0.31 12 4.3 0],0.0616,0.0028,0.0035,0.0097,i,0.2,0.57,1,0.5437,0.3786,0);
    y_p5(n) = y;
end
n = 0;
for i = 0:0.01:1
    n = n+1;
    y = objective_parametric([0.31 12 4.3 0],0.0616,0.0028,0.0035,0.0097,0.04,i,0.57,1,0.5437,0.3786,0);
    y_p6(n) = y;
end
n = 0;
for i = 0:0.008:0.8
    n = n+1;
    y = objective_parametric([0.31 12 4.3 0],0.0616,0.0028,0.0035,0.0097,0.04,0.2,i,1,0.5437,0.3786,0);
    y_p7(n) = y;
end

figure('Name','Parametric Study')
hold on
plot(x_p5,y_p5(1:100))
plot(x_p6,y_p6(1:100))
plot(x_p7,y_p7(1:100))
xlabel('Parameter Value')
ylabel('f(x)')
title('Parametric Analysis')
hold off
legend('Processing cost','Packaging cost','Profit margin')

%% Pareto Set

figure('Name','Attainable and Pareto Set')
hold on
for i = 0:0.3:3                                                      % From 0% to 20% Hazelnut,
    for j = 0:1.2:12                                                  % 0% to 20% Oat,
        for k = 0:1:10                                              % 0% to 20% Soy,
            for l = 0:0.8:8
                plot(f([i j k l]),u([i j k l]),'r.')    % Plot attainable set in red                                     
            end
        end
    end
end
xlabel('Quantity Sold')
ylabel('Utility')
title('Attainable Set and Pareto Set')
xlim([-0.805 -0.75])                                                   % Focus graph in the vicinity of the pareto set
ylim([-0.825 -0.76])                                                   % Focus graph in the vicinity of the pareto set                                       

fun = @(x)[f(x);u(x)];
lb = [0,0,0,0];
ub = [3,12,10,8];
rng default 
[x,fval] = paretosearch(fun,4,[],[],[],[],lb,ub,[]);                % Find the pareto set (fval) and the paramaters (x) to create it 
plot(fval(:,1),fval(:,2),'m*')                                      % Plot pareto set on top of attainable set
hold off

%% Re-formulation SQP Optimisation

options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1000);
A = [1 1 1 1];
b = [100];
Aeq = [];
beq = [];
x0 = [1 1 1 1];
lb = [0 0 0 0];
ub = [3 12 10 8];
nonlcon = @qd;
[xopt,fval] = fmincon(@objective,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
fprintf('The reformulated problem solved with SQP yielded a minima of %.3f with optimisers: %.3f, %.3f, %.3f, and %.3f',fval,xopt(1),xopt(2),xopt(3),xopt(4)) 

%% Sensitivity Analysis

x_x1 = linspace(0,3);
x_x2 = linspace(0,12);
x_x3 = linspace(0,10);
x_x4 = linspace(0,8);
n = 0;
y_x1 = [];
y_x2 = [];
y_x3 = [];
y_x4 = [];

for i = 0:0.03:3
    n = n+1;
    y_x1(n) = objective([i 12 4.3107 0]);
end
n = 0;
for i = 0:0.12:12
    n = n+1;
    y_x2(n) = objective([0.3131 i 4.3107 0]);  
end 
n = 0;
for i = 0:0.1:10
    n = n+1;
    y_x3(n) = objective([0.3131 12 i 0]);
end
n = 0;
for i = 0:0.08:8
    n = n+1;
    y_x4(n) = objective([0.3131 12 4.3107 i]);  
end 

figure('Name','Sensitivity Analysis')
hold on
plot(x_x1,y_x1(1:100))
plot(x_x2,y_x2(1:100))
plot(x_x3,y_x3(1:100))
plot(x_x4,y_x4(1:100))
hold off
xlabel('% Ingredient')
ylabel('f(x)')
title('Sensitivity Analysis')
legend('Hazelnut','Oat','Soy','Hemp')

%% Function definitions

function U = u(x)                                   % Utility function
global k4                                           % Call global variables 
global k5
global k6
global k7
global k8
global k9
global k10
global k11
global k12
global k13
global p1
global p2
global p3
global p4
global f1
global f2
global f3
global f4
pfc = x(1)*f1 + x(2)*f2 + x(3)*f3 + x(4)*f4;        % Calculate the fat content of a milk with x(1)% Haz, ...
ppc = x(1)*p1 + x(2)*p2 + x(3)*p3 + x(4)*p4;        % Calculate the protein content of a milk with x(1)% Haz, ...

uf = k4*pfc^2 + k5*pfc + k6;                        % Input fat content into utility fat curve
UF = ((abs(uf)+uf)/2);                              % Make -ve values = 0
up = k7*ppc^2 + k8*ppc + k9;                        % Input protein content into utility protein curve
UP = ((abs(up)+up)/2);                              % Make -ve values = 0
ut = (x(1)/(x(1)+x(2)+x(3)+x(4)))*k10 + (x(2)/(x(1)+x(2)+x(3)+x(4)))*k11...
        + (x(3)/(x(1)+x(2)+x(3)+x(4)))*k12...
        + (x(4)/(x(1)+x(2)+x(3)+x(4)))*k13;         % Calculate type utility: %(of product)Haz * Haz utility + ...
UT = ((abs(ut)+ut)/2);                              % Make -ve values = 0
U = -(1/3)*(UF + UP + UT);                          % Combine by weighting 1/3 each and put into -ve null form
end

function Qs = f(x)                                  % Quantity demanded function
global k1                                           % Call global variable
global k2
global k3
global c1
global c2
global c3
global c4
global c5
global c6                                           % Calculate price based on input costs and profit margin               
P = (x(1)*(c1+c5)+x(2)*(c2+c5)+x(3)*(c3+c5)+x(4)*(c4+c5)+c6)/(1-0.57);
qs = (k1*P^2 + k2*P + k3);                          % Input price into plant milk demand schedule curve
Qs = -(abs(qs)+qs)/2;                               % Make -ve values 0 THEN put into -ve null form
end

function obj = objective_initial(x)                 % Objective function
obj = -u(x)*f(x);                                         % Need to put in -ve null form again (-ve * -ve = +ve)
end

function U = U(x,P8,P9,P10,P11)      % FUNCTION REDEFINED FOR PARAMETRIC ANALYSIS                            
global k4                                          
global k5
global k6
global k7
global k8
global k9
global p1
global p2
global p3
global p4
global f1
global f2
global f3
global f4
pfc = x(1)*f1 + x(2)*f2 + x(3)*f3 + x(4)*f4;       
ppc = x(1)*p1 + x(2)*p2 + x(3)*p3 + x(4)*p4;        

uf = k4*pfc^2 + k5*pfc + k6;                       
UF = ((abs(uf)+uf)/2);                              
up = k7*ppc^2 + k8*ppc + k9;                        
UP = ((abs(up)+up)/2);                              
ut = (x(1)/(x(1)+x(2)+x(3)+x(4)))*P8 + (x(2)/(x(1)+x(2)+x(3)+x(4)))*P9...
        + (x(3)/(x(1)+x(2)+x(3)+x(4)))*P10...
        + (x(4)/(x(1)+x(2)+x(3)+x(4)))*P11;         
UT = ((abs(ut)+ut)/2);                              
U = -(1/3)*(UF + UP + UT);                          
end

function Qs = F(x,P1,P2,P3,P4,P5,P6,P7)  % FUNCTION REDEFINED FOR PARAMETRIC ANALYSIS                                    
global k1                                           
global k2
global k3                                                          
P = (x(1)*(P1+P5)+x(2)*(P2+P5)+x(3)*(P3+P5)+x(4)*(P4+P5)+P6)/(1-P7);
qs = (k1*P^2 + k2*P + k3);                          
Qs = -(abs(qs)+qs)/2;                               
end

function obj = objective_parametric(x,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11) % FUNCTION REDEFINED FOR PARAMETRIC ANALYSIS                      
obj = -U(x,P8,P9,P10,P11)*F(x,P1,P2,P3,P4,P5,P6,P7);                                  
end

function [c,ceq] = qd(x,k)
k = 0.76;                                           % Set -0.795 < k < 0
c = f(x)+k;                                      
ceq = [];
end

function obj = objective(x)                         % FINAL OBJECTIVE FUNCTION
obj = u(x);                               
end
