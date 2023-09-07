function [Lhf,Lhg,dh,L,u,r]=NonContSidFed(f,g,Hc,x,v)
% The NonContSidFed MATLAB function finds the lie
% derivatives and the decoupling matrix that has to be used
% to find the Sliding mode controller, by providing this
% function with:
% f : The system function f(x)
% g : The system output function g(x)
% x : The state vector x=[x1,x2,....,xn]
% v : The new input vector
% Hc : The output of the system hn(x)=xn
% The program has to give the user the following functions
% Lhf : The lie derivative of the vector h(x) along
% f(x) Lhf=[Lˆ{r−1}hfn]
% Lhg : The lie derivative of the vector field h
% along the function g; Lhg=[Lg1Lfhn]
% L : The lie derivatives L=[Lhf1 Lˆ{2}hf1.......
% Lˆ{r−1}hf1]
% x=f(x)+gu
L=[];
Lhf=[];Lhg=[];
if(nargin<4)
    error('Not enough input arguments!');
end
n1=length(Hc);
n=length(x);
k=1;kk=1;
while(k<n+1)
    dh=LieDer(Hc,x,n);
    Lg=dh*g;
    L=[L,dh*f];
    if(Lg~=0)
        r=k;
        Lhf=[Lhf;dh*f];
        Lhg=[Lhg;Lg];
        break;
    end
    Hc=dh*f;
    k=k+1;
end
u=inv(Lhg)*(Lhf-v);

function dh=LieDer(Hc,x,n)
dh=[];
for ii=1:n
    d=diff(Hc,x(ii));
    dh=[dh,d];
end
end
r=k;
end