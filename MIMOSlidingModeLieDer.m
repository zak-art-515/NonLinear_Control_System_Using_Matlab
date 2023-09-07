function [Lhf Lhg L r]=MIMOSlidingModeLieDer(fx,g,h,x)
% This programme finds a nonlinear controller which is the
% input output feedback linearization controller where the
% system should be given in symbolic format this function
% works well with SISO and MIMO systems which is written in
% compact form like dx=f(x)+gu;
% y=h(x);
% x : The state vector x=[x1,x2,...,xn]
% h : The output vector h(x)=[x1;x2;...;xn]
% fx : The f(x) the describes the system
% g : The g.u vector of the output system
% u : The output vector
% Lhf : The lie derivative of the vector file
% h along the function f(x)
% Lhf=[Lˆ{r1−1}fh1;Lˆ{r2−1}fh2;...;
% Lˆ{rn−1}fhn]
% Lhg : The lie derivative of the vector
% field h along the function g
% Lhg=[Lg1Lfh1, Lg2Lfh1,.....,LgnLfh1;
.........;
% Lg1Lfhn,Lg2Lfhn,......, LgnLfhn]
% r : The degree relative of the system
% r=[r1;r2;........;rn]
% L : The lie derivative vector
% : L=[Lhf1 Lhf2....... Lhfn;
% Lˆ2hf1 Lˆ2hf2....Lˆ2hfn;
%
% Lˆ(r1−1)hf1 Lˆ(r2−1)hf2.....
% Lˆ(rn−1)hfn]

if nargin <4
    error('Not enough input argument');
end
k=1;L=[];kk1=0;
Lhg=[];Lhf=[];vc=1;
nb=length(h);
r=zeros(1,nb);
while k<length(h)+1
    h1=h(k);
    for i=1:nb+1
        df=Lie_Derivative(h1,x); % this Lie derivative
        % function
        [lhf lhg]=solvelieder(df,fx,g); % solve for the g
        L=[L,lhf];
        [n b]=size(lhg);
        for ii=1:n
            d=any(lhg(ii,:) ~= 0);
        end
        if d==1;
            disp(['The relative degree of h',num2str(k)]);
            disp(['equal:=',num2str(i)]);
            r(1,k)=i;
            break;
        else
            h1=lhf;
        end
        if i==nb+1 && d==0
            disp(['The system dose not admit an input output feedback linearization']);
            return;
    end
end
Lhg=[Lhg;lhg];
Lhf=[Lhf;lhf];
k=k+1;
end