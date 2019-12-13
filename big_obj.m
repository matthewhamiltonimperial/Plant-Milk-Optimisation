function[z] = big_obj(x) %create multi-objective function (f)
z = zeros(1,2); % empty array
z(1) = lsq(x); %assign first objective function as f1 (saved as lsq)
z(2) = max(x); %assign second ojective function as f2 (saved as max)
end 
