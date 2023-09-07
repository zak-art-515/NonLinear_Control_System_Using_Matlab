function df=Lie_Derivative(h,x)
% The LieDerivative MATLAB function is used
% to find the jacobian vector of a given output
% h : Is the output function h(x)=[x1;x2;...;xn]
% x : The state vector x=[x1,x2,.....,xn]
% df : The jacobian of h along x
if nargin<2 & nargin==0
error('not enough input argument');
end
df=[];
n=length(x);
for ii=1:n
xx=x(ii);
dff=diff(h,xx);
df=[df,dff];
end
df;
end
function [lhf lhg]=solve_lie_der(df,fx,G)
% This equation solves the Lie derivativies that is
% multiplied by
% The output vectors Lhg=Lh*g*u where
% LH : The Lie derivativies of h(x) along the vector
% field fx
% fx : The vector field f(x) that describes the system
% G : The vector field of the input
LHg=[];
lhf=df*fx;
[n,b]=size(G);
for ii=1:b
Lgh=df*G(:,ii);
LHg=[LHg,Lgh];
end
lhg=LHg;
end
