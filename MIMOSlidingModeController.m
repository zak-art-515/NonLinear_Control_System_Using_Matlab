function [e,der,Surf,dSurf,Uc,S]=MIMOSlidingModeController(h,L,r,Lhg)
% The mimo nonlinear systems described by the following
% xp=f(x)+gu
% e : Is the error vector
% e=[e1 e2 ....en;de1 de2 ....den;.....;dne1
% dNe.... dNen]
% der : The derivative of the output vector
% der=[h1 h2 ...hn;d1h1 d1h2....dNhn;......;dNh1
% dNe....dNen]
% Surf : The Sliding Mode Surface for MIMO systems
% Surf=(d/dt+k)ˆ(r−1)e
% dSurf : The derivative of the sliding mode surface
% dSurf=d((d/dt+k)ˆ(r−1)e)/dt
% Uc : The Sliding mode control law
% h : The output vector y=[h1;h2;....;hN] for MIMO
% systems
% L : The lie derivative vector of the MIMO system
% L=[Lhf1 Lhf2....... Lhfn;
% Lˆ2hf1 Lˆ2hf2....Lˆ2hfn;
%
% Lˆrhf1 Lˆrhf2..... Lˆrhfn]
% r : The relative degree vector r=[r1,r2,r3,....,rn]
% Lhg : The lie derivative of h along the vector field g
% Lhg=[Lg1Lfh1, Lg2Lfh1,.....,LgnLfh1;
.........;
% Lg1Lfhn,Lg2Lfhn,......, LgnLfhn]
e=[];kk=1;
rr=max(r);
der=[];S=[];Surf=[];Ss=[];
nb=length(h);d=[];qq=1;q=1;
if (r==1)& (nb==1)
error('The output vector h should be of length >=2');
end
if nargin <4
error('Not enough input argument');
end
LL=sym(zeros(rr+1,nb));
der=sym(zeros(rr+1,nb));
syms Uc U
LL(1,:)=h;
k=[];
K=sym(zeros(nb,rr));
KK=sym(zeros(nb,rr+1));
KK(:,1)=sym('1');
%K(:,1)=1;
for ll=1:nb
R=r(ll);
k=sym(zeros(1,R-1));
for jj=1:R-1
eval(sprintf('syms k%d',jj));
k(:,jj)=sprintf('k%d',jj);
end
K(ll,:)=[1 k];
end
KK(:,nb:rr+1)=K;
F=subs(KK(:,nb+1:rr+1),{1},{0});
KK(:,nb+1:rr+1)=F;
%% The functions sgnS1, Sgn(S2),....sgn(Sn)
% The parameters kp1, kp2,....kp3
sgns=sym(zeros(1,nb));
kp=sym(zeros(1,nb));
for jj=1:nb
eval(sprintf('syms sgnS%d',jj));
sgns(:,jj)=sprintf('sgnS%d',jj);
end
for jj=1:nb
eval(sprintf('syms kp%d',jj));
kp(:,jj)=sprintf('kp%d',jj);
end
%%
for jj=1:nb
R=r(jj);
dd=sym(zeros(1,R));
for ii=1:R+1
eval(sprintf('syms d%dyr%d',ii,kk));
dd(:,ii)=sprintf('d%dyr%d',ii,kk);
end
d=[d,dd];
kk=kk+1;
end
for jj=1:nb
R=r(jj);
for mm=2:R+1
LL(mm,jj)=L(qq);
qq=qq+1;
end
for nn=1:R+1
der(nn,jj)=d(q);
q=q+1;
end
end
e=flipud(der-LL);
for bb=1:nb
Kd=fliplr(KK(bb,nb:rr+1));
ee=e(1:rr,bb);
nm=kp(bb)*sgns(bb);
sc=Kd*ee+kp(bb)*sgns(bb);
S=[S;sc];
Ss=[Ss;Kd*ee];
end
es=flipud(e);
for cc=1:nb
R=r(cc);
dK=fliplr(KK(cc,nb:R+1));
esf=(es(1:R,cc));
surface=dK*(esf);
Surf=[Surf;surface];
end
dSurf=Ss;
Surf=Surf;
[nh bh]=size(Lhg);
disp(['The sliding mode control law for MIMO systems']);
disp(['−−−−−−Is given by Uc=:inv(Lhg)*(S)−−−−−−−−−−−−−']);
disp(['−−−−−−The function S:=−−−−−−−−−−']);
S
disp(['−−−−−−The Matrix Lhg:=−−−−−−−−−−']);
Lhg
for ii=1:nh
invg=any(Lhg(ii,:)~=0);
end
if(invg~=0)
disp(['The sliding mode control law for MIMO systems=:']);
Uc=inv(Lhg)*(S)
else
disp(['The sliding mode controller is not possible']);
disp(['for this output vector try to choose another one!']);
end

