function[obj1] = lsq(x) %create minimise sum of least squares function (f1a)

T = readmatrix('Matlab_data.xlsx'); %read excel table 
T(1:6,1:20);
DF1 = T(1,2:9); %create array of deficiency factors for relevant nutrients 
x = x'
NF1 = [T(3,2:9) %create matrix of nutrient factors
    T(4,2:9)
    T(5,2:9)
    T(6,2:9)];


NFX1 = x.*NF1; %multiply each variable by relevant nutrient factors 
D1 = T(2,2:9); %create array of dairy values 
A1 = [sum(NFX1(1:4,1)),sum(NFX1(1:4,2)),sum(NFX1(1:4,3)),sum(NFX1(1:4,4)),sum(NFX1(1:4,5)),sum(NFX1(1:4,6)),sum(NFX1(1:4,7)),sum(NFX1(1:4,8))]; %create array of the sum of each nutrient as a proportion of each crop
ob1 = DF1.*(D1-A1).^2; %find difference between sum of plant-based nutrient value and dairy, square it to remove negative values, and multiply by relevant deficiency factor
obj1 = sum(ob1(1,1:4)); %sum the total squared differences 
end