function [Surf,dSurf,dd,K,Uc]=SlidingModeTerms(Hc,L,r,Lhg)
% The SlidingModeTerms MATLAB function finds the sliding mode
% surface and its derivatives for SISO nonlinear systems
% That takes the following form:
% xp=f(x)+gu
% y=h(x)
% Hc : The output of the system Hc=x
% L : The lie derivatives L=[Lhf1 Lˆ2hf1.......
% Lˆ{r−1}hf1]
% r : The relative degree of the system r=r1
% Surf : Sliding Mode surface for SISO systems
% Surf=(d/dt+k)ˆ(r−1)e
% dSurf : The derivative of the sliding surface
% Surf=d((d/dt+k)ˆ(r−1)e)/dt
% K : The controller parameter
% vector K=[1 k1 k2,....kn]
k=[];
syms kp Uc U sgnS
if (r==1)
    K=1;
else
    k=sym(zeros(1,r-1));
    for jj=1:r-1
        eval(sprintf('syms k%d',jj));
        k(:,jj)=sprintf('k%d',jj);
    end
    K=[1 k];
end
dd=sym(zeros(1,r));
for ii=1:r+1
    eval(sprintf('syms d%dyr',ii));
    dd(:,ii)=sprintf('d%dyr',ii);
end
e=[dd(1)-Hc];
%s=[dd(r)−Hc];
s=[];
k=1;
for kk=2:r
    e=[e;(dd(kk)-L(k))];
    k=k+1;
end
LL=fliplr(L);
dl=fliplr(dd);
for kk=1:1:r
    s=[s;(dl(kk)-LL(kk))];
end
S=flipud(e);
Sp=s;
Uc=Lhg*U;
disp(['The sliding mode control law for SISO systems=:']);
Uc=K*Sp+kp*sgnS
disp(['The Sliding mode surface:=']);
Surf=K*S
disp(['The derivative of Sliding mode surface:=']);
dSurf=K*Sp