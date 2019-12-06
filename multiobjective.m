clear
clc

%% Read raw data from .csv file

Fat = csvread('Fat_Protein_Type.csv',0,0,[0,0,3,0]);                        % [0g/8oz,2.5g/8oz,5g/8oz,8g/8oz]
fat_utility = csvread('Fat_Protein_Type.csv',0,1,[0,1,3,1]);                % Utility score of fat contents
Fat_Utility = normalize(fat_utility,'range');
Protein = csvread('Fat_Protein_Type.csv',0,2,[0,2,4,2]);                    % [1g/8oz,2g/8oz,5g/8oz,8g/8oz,9g/oz]
protein_utility = csvread('Fat_Protein_Type.csv',0,3,[0,3,4,3]);            % Utility score of protein contents
Protein_Utility = normalize(protein_utility,'range');
type_utility = csvread('Fat_Protein_Type.csv',0,4,[0,4,3,4]);               % [Hazelnut, Oat, Soy, Hemp]
Type_Utility = normalize(type_utility,'range');
Fat_Content = csvread('Fat_Protein_Type.csv',0,5,[0,5,3,5]);                % [Hazelnut fat, Oat fat, Soy fat, Hemp fat]
Protein_Content = csvread('Fat_Protein_Type.csv',0,6,[0,6,3,6]);            % [Hazelnut protein, Oat protein, Soy protein, Hemp protein]

input_costs = csvread('Optimisation_Main_2.0.csv',0,2,[0,2,0,8]);           % [Hazelnut,Oat,Soy,Hemp,Packaging,Production,Margin]
quantity_sold = csvread('Optimisation_Main_2.0.csv',0,1,[0,1,6,1]);
Quantity_Sold = normalize(quantity_sold,'range'); 
Price_Per_Litre = csvread('Optimisation_Main_2.0.csv',0,0,[0,0,6,0]);                     %Price per Litre
coeff_fat = polyfit(Fat,Fat_Utility,2); 
coeff_protein = polyfit(Protein,Protein_Utility,2); 
coeff_demand = polyfit(Price_Per_Litre,Quantity_Sold,2); 

global c1
global c2
global c3
global c4
global c5
global c6
global k1
global k2
global k3
global k4
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

c6 = input_costs(1,5);                                                      %Packaging Costs
c5 = input_costs(1,6);                                                      %Plant processing costs
c1 = input_costs(1,1);                                                      %Hazelnut cost 
c2 = input_costs(1,2);                                                      %Oat cost
c3 = input_costs(1,3);                                                      %Soy cost
c4 = input_costs(1,4);          %Hemp cost 
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


%% Data visualisation 

figure('Name','Relationships between protein, fat and type on Utility Coefficient')
subplot(3,1,1)
plot(Fat,Fat_Utility, 'r.');
hold on
X_fat = linspace(0,40);
Y_fat = polyval(coeff_fat,X_fat);
plot(X_fat,Y_fat)
hold off
ylim([-0.1 1.1])
title('Effect of Fat content on Utility')
xlabel('Fat content (grams/litre)')
ylabel('Utility coefficient')
subplot(3,1,2)
plot(Protein,Protein_Utility, 'r.');
hold on
X_protein = linspace(0,40);
Y_protein = polyval(coeff_protein,X_protein);
plot(X_protein,Y_protein)
hold off
ylim([-0.1 1.1])
title('Effect of Protein content on Utility')
xlabel('Protein content (grams/litre)')
ylabel('Utility coefficient')
subplot(3,1,3)
bar(type_utility)
set(gca,'xticklabel',{'Hazelnut','Oat','Soy','Hemp'})
title('Effect of Type on Utility')
ylabel('Utility coefficient')


figure('Name','Plant Milk demand schedule')
plot(Price_Per_Litre,Quantity_Sold, 'b.');
hold on
X_demand = linspace(0,5);
Y_demand = polyval(coeff_demand,X_demand);
plot(X_demand,Y_demand, 'r')
hold off
ylim([-0.1 1.1])
title('Plant Milk Demand Schedule')
xlabel('Price ($/litre)')
ylabel('Normalised Quantity Sold')


figure('Name','Surface Utility Plots')
subplot(3,2,1) % Hazelnuts vs Oat 
[X,Y] = meshgrid(0:.5:20);
Z = (1/3).*((((abs(k4.*(X.*f1 + Y.*f2).^2 + k5.*(X.*f1 + Y.*f2) + k6)+(k4.*(X.*f1 + Y.*f2).^2 + k5.*(X.*f1 + Y.*f2) + k6))./2))...
     + (((abs(k7.*(X.*p1+Y.*p2).^2 + k8.*(X.*p1+Y.*p2) + k9)+(k7.*(X.*p1+Y.*p2).^2 + k8.*(X.*p1+Y.*p2) + k9))./2))...
     + ((X./(X+Y)).*k10 + (Y./(X+Y)).*k11));
s = surf(X,Y,Z);
title('Hazelnut vs Oat')
xlabel('% Hazelnut')
ylabel('% Oat')
zlabel('Utility')
subplot(3,2,2) % Hazelnuts vs Soy
[X,Y] = meshgrid(0:.5:20);
Z = (1/3).*((((abs(k4.*(X.*f1 + Y.*f3).^2 + k5.*(X.*f1 + Y.*f3) + k6)+(k4.*(X.*f1 + Y.*f3).^2 + k5.*(X.*f1 + Y.*f3) + k6))./2))...
     + (((abs(k7.*(X.*p1+Y.*p3).^2 + k8.*(X.*p1+Y.*p3) + k9)+(k7.*(X.*p1+Y.*p3).^2 + k8.*(X.*p1+Y.*p3) + k9))./2))...
     + ((X./(X+Y)).*k10 + (Y./(X+Y)).*k12));
s = surf(X,Y,Z);
title('Hazelnut vs Soy')
xlabel('% Hazelnut')
ylabel('% Soy')
zlabel('Utility')
subplot(3,2,3) % Hazelnuts vs Hemp
[X,Y] = meshgrid(0:.5:20);
Z = (1/3).*((((abs(k4.*(X.*f1 + Y.*f4).^2 + k5.*(X.*f1 + Y.*f4) + k6)+(k4.*(X.*f1 + Y.*f4).^2 + k5.*(X.*f1 + Y.*f4) + k6))./2))...
     + (((abs(k7.*(X.*p1+Y.*p4).^2 + k8.*(X.*p1+Y.*p4) + k9)+(k7.*(X.*p1+Y.*p4).^2 + k8.*(X.*p1+Y.*p4) + k9))./2))...
     + ((X./(X+Y)).*k10 + (Y./(X+Y)).*k13));
s = surf(X,Y,Z);
title('Hazelnut vs Hemp')
xlabel('% Hazelnut')
ylabel('% Hemp')
zlabel('Utility')
subplot(3,2,4) % Oat vs Soy
[X,Y] = meshgrid(0:.5:20);
Z = (1/3).*((((abs(k4.*(X.*f2 + Y.*f3).^2 + k5.*(X.*f2 + Y.*f3) + k6)+(k4.*(X.*f2 + Y.*f3).^2 + k5.*(X.*f2 + Y.*f3) + k6))./2))...
     + (((abs(k7.*(X.*p2+Y.*p3).^2 + k8.*(X.*p2+Y.*p3) + k9)+(k7.*(X.*p2+Y.*p3).^2 + k8.*(X.*p2+Y.*p3) + k9))./2))...
     + ((X./(X+Y)).*k11 + (Y./(X+Y)).*k12));
s = surf(X,Y,Z);
title('Oat vs Soy')
xlabel('% Oat')
ylabel('% Soy')
zlabel('Utility')
subplot(3,2,5) % Oat vs Hemp
[X,Y] = meshgrid(0:.5:20);
Z = (1/3).*((((abs(k4.*(X.*f2 + Y.*f4).^2 + k5.*(X.*f2 + Y.*f4) + k6)+(k4.*(X.*f2 + Y.*f4).^2 + k5.*(X.*f2 + Y.*f4) + k6))./2))...
     + (((abs(k7.*(X.*p2+Y.*p4).^2 + k8.*(X.*p2+Y.*p4) + k9)+(k7.*(X.*p2+Y.*p4).^2 + k8.*(X.*p2+Y.*p4) + k9))./2))...
     + ((X./(X+Y)).*k11 + (Y./(X+Y)).*k13));
s = surf(X,Y,Z);
title('Oat vs Hemp')
xlabel('% Oat')
ylabel('% Hemp')
zlabel('Utility')
subplot(3,2,6) % Soy vs Hemp
[X,Y] = meshgrid(0:.5:20);
Z = (1/3).*((((abs(k4.*(X.*f3 + Y.*f4).^2 + k5.*(X.*f3 + Y.*f4) + k6)+(k4.*(X.*f3 + Y.*f4).^2 + k5.*(X.*f3 + Y.*f4) + k6))./2))...
     + (((abs(k7.*(X.*p3+Y.*p4).^2 + k8.*(X.*p3+Y.*p4) + k9)+(k7.*(X.*p3+Y.*p4).^2 + k8.*(X.*p3+Y.*p4) + k9))./2))...
     + ((X./(X+Y)).*k12 + (Y./(X+Y)).*k13));
s = surf(X,Y,Z);
title('Soy vs Hemp')
xlabel('% Soy')
ylabel('% Hemp')
zlabel('Utility')


figure('Name','Surface Demand Plots')
subplot(3,2,1) % Hazelnuts vs Oat 
[X,Y] = meshgrid(0:.5:20);
Z = k1*((X.*(c1+c5)+Y.*(c2+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c1+c5)+Y.*(c2+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);
title('Hazelnut vs Oat')
xlabel('% Hazelnut')
ylabel('% Oat')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Sold')
subplot(3,2,2) % Hazelnuts vs Soy 
[X,Y] = meshgrid(0:.5:20);
Z = k1*((X.*(c1+c5)+Y.*(c3+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c1+c5)+Y.*(c3+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);
title('Hazelnut vs Soy')
xlabel('% Hazelnut')
ylabel('% Soy')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Sold')
subplot(3,2,3) % Hazelnuts vs Hemp 
[X,Y] = meshgrid(0:.5:20);
Z = k1*((X.*(c1+c5)+Y.*(c4+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c1+c5)+Y.*(c4+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);
title('Hazelnut vs Hemp')
xlabel('% Hazelnut')
ylabel('% Hemp')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Sold')
subplot(3,2,4) % Oat vs Soy 
[X,Y] = meshgrid(0:.5:20);
Z = k1*((X.*(c2+c5)+Y.*(c3+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c2+c5)+Y.*(c3+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);
title('Oat vs Soy')
xlabel('% Oat')
ylabel('% Soy')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Sold')
subplot(3,2,5) % Oat vs Hemp 
[X,Y] = meshgrid(0:.5:20);
Z = k1*((X.*(c2+c5)+Y.*(c4+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c2+c5)+Y.*(c4+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);
title('Oat vs Hemp')
xlabel('% Oat')
ylabel('% Hemp')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Sold')
subplot(3,2,6) % Soy vs Hemp 
[X,Y] = meshgrid(0:.5:20);
Z = k1*((X.*(c3+c5)+Y.*(c4+c5)+c6)/(1-0.57)).^2 + k2*((X.*(c3+c5)+Y.*(c4+c5)+c6)/(1-0.57)) + k3;
s = surf(X,Y,Z);
title('Soy vs Hemp')
xlabel('% Soy')
ylabel('% Hemp')
zlim([-0.1 1.1])
zlabel('Normalised Quantity Sold')

%% Pareto Set Visualisation

figure('Name','Attainable Set')
hold on
for i = 0:4:20
    for j = 0:4:20
        for k = 0:4:20
            for l = 0:4:20
                if f([i j k l]) < -0.2 && u([i j k l]) < -0.2
                    plot(f([i j k l]),u([i j k l]),'r-')
                end                                
            end
        end
    end
end
for i = 1:0.4:2
    for j = 10:2:20
        for k = 0:0.8:4
            for l = 0:0.4:2
                if f([i j k l]) < -0.2 && u([i j k l]) < -0.2
                    plot(f([i j k l]),u([i j k l]),'b.')
                end                                
            end
        end
    end
end
xlabel('Quantity Sold')
ylabel('Utility')
title('Attainable Set and Pareto Set')
xlim([-0.8 -0.5])
ylim([-0.9 -0.5])
xlabel('Quantity Demanded')
ylabel('Utility')

fun = @(x)[f(x);u(x)];
lb = [0,0,0,0];
ub = [20,20,20,20];
rng default 
[x,fval] = paretosearch(fun,4,[],[],[],[],lb,ub,[]);
plot(fval(:,1),fval(:,2),'m*')
hold off

%% fmincon Solver

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0 0 0];
ub = [17 17 17 17];
x0 = [1 1 1 1];
xopt = fmincon(@objective,x0,A,b,Aeq,beq,lb,ub)

% figure('Name','fmincon solution')
% bar(xopt)
% set(gca,'xticklabel',{'Hazelnut','Oat','Soy','Hemp'})
% title('fmincon solution')
% ylabel('% Ingredient')

%% Sensitivity Analysis



%% Functions

function U = u(x)

global k4
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

pfc = x(1)*f1 + x(2)*f2 + x(3)*f3 + x(4)*f4;
ppc = x(1)*p1 + x(2)*p2 + x(3)*p3 + x(4)*p4; 

uf = k4*pfc^2 + k5*pfc + k6;
UF = ((abs(uf)+uf)/2);
up = k7*ppc^2 + k8*ppc + k9;
UP = ((abs(up)+up)/2);
ut = (x(1)/(x(1)+x(2)+x(3)+x(4)))*k10 + (x(2)/(x(1)+x(2)+x(3)+x(4)))*k11...
        + (x(3)/(x(1)+x(2)+x(3)+x(4)))*k12...
        + (x(4)/(x(1)+x(2)+x(3)+x(4)))*k13; 
UT = ((abs(ut)+ut)/2);
U = -(1/3)*(UF + UP + UT);
end

function Qs = f(x) 
global k1
global k2
global k3
global c1
global c2
global c3
global c4
global c5
global c6
P = (x(1)*(c1+c5)+x(2)*(c2+c5)+x(3)*(c3+c5)+x(4)*(c4+c5)+c6)/(1-0.57);
qs = (k1*P^2 + k2*P + k3); 
Qs = -(abs(qs)+qs)/2;
end

function obj = objective(x)
obj = -f(x)*u(x);
end
