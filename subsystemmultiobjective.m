
function z = subsystemmultiobjective(A)
    %1.hazelnut,2.hemp,3.oat,4.soy
    Y = [1.5 1.5 4.9 2.81];  %yield for each crop given in kg/m^2
    P = [2770 1010.33 181.5 264.31];  %selling price given in GBP/kg
    IC = 45; %irrigation cost GBP/m^3
    CWR = [366 350 329.35 575]; %constants for crop water requirement
    CS = [8 5 1 38]; %constants for carbon sequestered by plant biomass


    TY = Y(1)*P(1)*A(1) - CWR(1)*A(1)*IC/10000000000+ Y(2)*P(2)*A(2) - CWR(2)*A(2)*IC/10000000000 + Y(3)*P(3)*A(3) - CWR(3)*A(3)*IC/10000000000 + Y(4)*P(4)*A(4) - CWR(4)*A(4)*IC/10000000000;
    %TY = Y(1)*P(1)*A(1) + Y(2)*P(2)*A(2) + Y(3)*P(3)*A(3) + Y(4)*P(4)*A(4);

    WU = CWR(1)*A(1) + CWR(2)*A(2) +CWR(3)*A(3) + CWR(3)*A(3);

    TCS = CS(1)*A(1) + CS(2)*A(2) + CS(3)*A(3) + CS(4)*A(4);


    z(1) = -TY;
    z(2) = WU;
    z(3) = -TCS;

end



