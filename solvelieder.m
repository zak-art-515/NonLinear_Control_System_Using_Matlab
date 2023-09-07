function [lhf lhg]=solvelieder(Lh,fx,g)
% The solvelieder MATLAB function is used to find
% the lie derivatives of the functions f(x) and g(x)
% along the vector field h(x)
% Lh : The jacobian vector of h along x
% fx : The function f(x)
% g : The input function g(x)
LHg=[];
lhf=Lh*fx;
[n,b]=size(g);

for ii=1:b
    Lgh=Lh*g(:,ii);
    LHg=[LHg,Lgh];
end

lhg=LHg;