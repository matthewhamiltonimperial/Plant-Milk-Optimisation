function[obj4] = goal(x) %create function for minimising individual nutrient difference with dairy (f1b)

T = readmatrix('Matlab_data.xlsx'); %read excel table 
T(1:6,1:20); 
DF1 = T(1,2:9); %create array of relevant deficiency factors
D1 = T(2,2:9); %create array of relevant dairy nutrient values

NF1 = [T(3,2:9) %create matrix of relevant nutrient factors
    T(4,2:9)
    T(5,2:9)
    T(6,2:9)];

%put table into format for fminimax optimisation function

a = T(3:6,2); %create individual nutrient equations using relevant nutrient factors. i.e. a is for protien
b = T(3:6,3);
c = T(3:6,4);
d = T(3:6,5);
e = T(3:6,6);
f = T(3:6,7);
g = T(3:6,8);
h = T(3:6,9);

a0 = D1(1,1); %create individual nutrient goals using dairy values 
b0 = D1(1,2);
c0 = D1(1,3);
d0 = D1(1,4);
e0 = D1(1,5);
f0 = D1(1,6);
g0 = D1(1,7);
h0 = D1(1,8);

diff = [(a0-x*a)^2,(b0-x*b)^2,(c0-x*c)^2,(d0-x*d)^2,(e0-x*e)^2,(f0-x*f)^2,(g0-x*g)^2,(h0-x*h)^2]; %create array of squared differences of individual nutrients with individual nutrient goals 
obj4 = [DF1.*diff]; %add deficiency factor weighting 

end
