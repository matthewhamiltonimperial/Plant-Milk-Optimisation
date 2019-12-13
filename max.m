function[obj2] = max(x) %create maximising healthy nutrients function for f2

T = readmatrix('Matlab_data.xlsx'); %read excel table 
T(1:6,1:20);

DF2 = T(1,10:17); %create deficiency factor for nutrients (i) 10 to 17
x = x'; 
NF2 = -[T(3,10:17); %create matrix of null form nutrient coefficients
    T(4,10:17);
    T(5,10:17);
    T(6,10:17)];
NFX2 = DF2.*(x.*NF2); %multiply by relevant deficiency factors 
A2 = [sum(NFX2(1,1:8)),sum(NFX2(2,1:8)),sum(NFX2(3,1:8)),sum(NFX2(4,1:8))]; %sum each variable coefficients 
obj2 = sum(A2(1,1:4)) %create array objective function 
end