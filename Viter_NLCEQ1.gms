*value functiuon iteration + NLCEQ for credit card puzzle paper
scalars
beta discount rate /0.99833/
rho1 wage auto correlation /0.9/
sigma wage standard deviation/0.01/
rho2 CRRA /2/
age household age /30/
rs saving interest rate /0.00033/
avy permanent average monthly income /1/
rb credit card debt interest rate /0.013/
Blim credit limit of all credit cards /10/
epsilon A small numbner /0.001/
myear mortgage length /120/

;

* Step 1: Define the transformed deterministic problem for renters


Sets t2 time month/1*480/;

parameters
y(t2) monthly income
alpha parameters of preference to housing service
M     mortgage payment
betas(t2) discount factor;
y(t2) = 1;
betas(t2) = beta**(ord(t2)-1);

variables
objR objective creterion of renters
uR(t2) utility of renting

s(t2) saving
b(t2) credit debt
;

positive variables
c(t2) non-durable goods consumption
p(t2) cash payment of non-durable goods
rent(t2) cash payment of rent
;

Equations
objR1 Objective function for renters
utilR(t2) utility function of renting
savR(t2) process of saving growth of house renters
bor(t2) process of credit card debt grwoth
debt0(t2) lower bound of credit card debt
debt1(t2) upper bound of credit card debt
cashR(t2)  saving is nonnegative--renters
debtT(t2)    credit card debt of the last month of one's life
;

savR(t2)..  s(t2+1) =e= (1+rs)*(s(t2)+y(t2)-rent(t2)-p(t2));
bor(t2) .. b(t2+1) =e= (1+rb)*(b(t2)+c(t2)-p(t2));
debt0(t2) .. b(t2)+c(t2)-p(t2) =g= 0;
debt1(t2) .. b(t2)+c(t2)-p(t2) =l= Blim;
cashR(t2) .. s(t2)+y(t2)-rent(t2)-p(t2) =g= 0;
debtT(t2)$(ord(t2)=card(t2)) .. b(t2) =e= epsilon;
utilR(t2) .. uR(t2) =e= (rent(t2)**alpha * c(t2)**(1-alpha))**(1-rho2)/(1-rho2);
objR1 .. objR =e= sum(t2, betas(t2)*uR(t2))

* Step 2: Define the transformed deterministic problem for house owners without mortgage
* because the borrowing process are the same for owners and renters, we dont have to model the process twice.

variables
uO(t2) utility of owning house
objO objective creterion of house owner
;

Equations
objO1 Objective function for house owners without mortgage
utilO(t2) utility function of owing
savO(t2) process of saving growth of house owners
cashO(t2)  saving is nonnegative--house owners without mortgage
;

savO(t2) .. s(t2+1) =e= (1+rs)*(s(t2)+y(t2)-p(t2));
cashO(t2).. s(t2)+y(t2)-p(t2) =g= 0;
utilO(t2) .. uO(t2) =e= (M**alpha * c(t2)**(1-alpha))**(1-rho2)/(1-rho2);
objO1 .. objO =e= sum(t2, betas(t2)*uO(t2));

* bound constraints
s.lo(t2) = epsilon;
s.up(t2) = 20;
b.lo(t2) = epsilon;
b.up(t2) = Blim;
p.lo(t2) = epsilon;
c.lo(t2) = epsilon;
c.up(t2) = 10;
rent.lo(t2) = epsilon;
rent.up(t2) = 0.9;


* initial guess
s.l(t2) = epsilon;
b.l(t2) = epsilon;
rent.l(t2) = 0.2*y(t2);
c.l(t2) = 0.5*y(t2);
p.l(t2) = c.l(t2);
alpha = 0.5;
M = 0.5;



uR.l(t2) = (rent.l(t2)**alpha * c.l(t2)**(1-alpha))**(1-rho2)/(1-rho2);
objR.l = sum(t2, betas(t2)*uR.l(t2));

uO.l(t2) = (M**alpha * c.l(t2)**(1-alpha))**(1-rho2)/(1-rho2);
objO.l = sum(t2,betas(t2)*uO.l(t2));

Option limrow=0,limcol=0,solprint=off;
OPTION ITERLIM = 500000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION NLP= conopt;
option decimals = 7;

model dpRenter "value function of renters at t=myear+1" /objR1,utilR,savR,bor,debt0,debt1,cashR,debtT/
      dpOwner "value function of houseowner at t=myear+1"/objO1,utilO,savO,bor,debt0,debt1,cashO,debtT/;

*******************************************************
* Step 2: Optimization step

set v1 grid space/1*6/
alias(v1,v2,v3,v4,v5,d1,d2,d3,d4,d5);

*approximation domain

scalars
svmin lower bound of s /0/
svmax upper bound of s /20/
bmin lower bound of b /0/
bmax upper bound of b /10/
amin lower bound of alpha /0.1/
amax upper bound of alpha /0.9/
mpmin lower bound of mortgage payment /0.05/
mpmax upper bound of mortgage payment /0.9/
ymin lower bound of income /0.3/
ymax upper bound of income /2/
;

*Chebyshev parameters

parameters zsv(v1) saving
         zb(v2)
         za(v3)
         zmp(v4)
         zy(v5)
         extsv
         extb
         exta
         extmp
         exty
         extsvmin
         extsvmax
         extbmin
         extbmax
         extamin
         extamax
         extmpmin
         extmpmax
         extymin
         extymax
         dzsv
         dzb
         dza
         dzmp
         dzy
;

zsv(v1) = - cos((2*ord(v1)-1)/(2*card(v1))*pi);
zb(v2) = - cos((2*ord(v2)-1)/(2*card(v2))*pi);
za(v3) = - cos((2*ord(v3)-1)/(2*card(v3))*pi);
zmp(v4) = - cos((2*ord(v4)-1)/(2*card(v4))*pi);
zy(v5) =  - cos((2*ord(v5)-1)/(2*card(v5))*pi);

extsv = -(zsv('1')+1)*(svmax-svmin)/2/zsv('1');
extb = -(zb('1')+1)*(bmax-bmin)/2/zb('1');
exta = -(za('1')+1)*(amax-amin)/2/za('1');
extmp = -(zmp('1')+1)*(mpmax-mpmin)/2/zmp('1');
exty = -(zy('1')+1)*(ymax-ymin)/2/zy('1');

extsvmin = svmin - extsv;
extbmin = bmin - extb;
extamin = amin - exta;
extmpmin = mpmin - extmp;
extymin = ymin -exty;

extsvmax = svmax + extsv;
extbmax = bmax + extb;
extamax = amax + exta;
extmpmax = mpmax - extmp;
extymax = ymax - exty;

dzsv = 2/(extsvmax-extsvmin);
dzb = 2/(extbmax-extbmin);
dza = 2/(extamax - extamin);
dzmp = 2/(extmpmax - extmpmin);
dzy = 2/ (extymax - extymin);

parameters
sinit(v1) Grid Saving
binit(v2) Grid Credit Card Debt
ainit(v3) Grid preference coefficient alpha
mpinit(v4) Grid of mortgage payment
yinit(v5) grid of income
;

sinit(v1) = extsvmin + (1+zsv(v1))/dzsv;
binit(v2) = extbmin + (1+zb(v2))/dzb;
ainit(v3) = extamin + (1+za(v3))/dza;
mpinit(v4) = extmpmin + (1+zmp(v4))/dzmp;
yinit(v5) = extymin + (1+zy(v5))/dzy;

Parameters
vosol(v1,v2,v3,v4,v5) value of the house owners value function
vrsol(v1,v2,v3,v5) value of the renters value function
;

loop((v1,v2,v3,v4,v5),
  s.fx('1') = sinit(v1);
  b.fx('1') = binit(v2);
  alpha = ainit(v3);
  M = mpinit(v4);
  y('1') = yinit(v5);
  loop (t2$(ord(t2) lt card(t2)),
    y(t2+1) = y(t2)**rho1;
  );

  solve dpOwner using nlp maximizing objO;
  vosol(v1,v2,v3,v4,v5) = objO.l;
  solve dpRenter using nlp maximizing objR;
  vrsol(v1,v2,v3,v5) = objR.l;
);

************************************************
* Step 3: Approximation step

parameters TTs(v1,d1)
           TTb(v2,d2)
           TTa(v3,d3)
           TTmp(v4,d4)
           TTy(v5,d5)
;
TTs(v1,'1') = 1;
TTs(v1,'2') = zsv(v1);
loop(d1$(ord(d1)>=2 and ord(d1)<card(d1)),
  TTs(v1,d1+1) = 2*zsv(v1)*TTs(v1,d1) - TTs(v1,d1-1);
);

TTb(v2,'1') = 1;
TTb(v2,'2') = zb(v2);
loop(d2$(ord(d2)>=2 and ord(d2)<card(d2)),
  TTb(v2,d2+1) = 2*zb(v2)*TTb(v2,d2) - TTb(v2,d2-1);
);

TTa(v3,'1') = 1;
TTa(v3,'2') = za(v3);
loop(d3$(ord(d3)>=2 and ord(d3)<card(d3)),
  TTa(v3,d3+1) = 2*za(v3)*TTa(v3,d3) - TTa(v3,d3-1);
);

TTmp(v4,'1') = 1;
TTmp(v4,'2') = zmp(v4);
loop(d4$(ord(d4)>=2 and ord(d4)<card(d4)),
  TTmp(v4,d4+1) = 2*zmp(v4)*TTmp(v4,d4) - TTmp(v4,d4-1);
);

TTy(v5,'1') = 1;
TTy(v5,'2') = zy(v5);
loop(d5$(ord(d5)>=2 and ord(d5)<card(d5)),
  TTy(v5,d5+1) = 2*zy(v5)*TTy(v5,d5) - TTy(v5,d5-1);
);

parameter coefsVo(d1,d2,d3,d4,d5) approximation coefficients of house owners value function;
parameter coefsVr(d1,d2,d3,d5) approximation coefficients of renters value function;

coefsVo(d1,d2,d3,d4,d5) = 0;
coefsVo(d1,d2,d3,d4,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)+ord(d5)<=card(d1)+4) =
    32*sum((v1,v2,v3,v4,v5), TTs(v1,d1)*TTb(v2,d2)*TTa(v3,d3)*TTmp(v4,d4)*TTy(v5,d5)*vosol(v1,v2,v3,v4,v5)) / (card(v1)*card(v2)*card(v3)*card(v4)*card(v5));

coefsVo('1',d2,d3,d4,d5) = coefsVo('1',d2,d3,d4,d5)/2;
coefsVo(d1,'1',d3,d4,d5) = coefsVo(d1,'1',d3,d4,d5)/2;
coefsVo(d1,d2,'1',d4,d5) = coefsVo(d1,d2,'1',d4,d5)/2;
coefsVo(d1,d2,d3,'1',d5) = coefsVo(d1,d2,d3,'1',d5)/2;
coefsVo(d1,d2,d3,d4,'1') = coefsVo(d1,d2,d3,d4,'1')/2;

coefsVr(d1,d2,d3,d5) = 0;
coefsVr(d1,d2,d3,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d5)<=card(d1)+3) =
    16*sum((v1,v2,v3,v5), TTs(v1,d1)*TTb(v2,d2)*TTa(v3,d3)*TTy(v5,d5)*vrsol(v1,v2,v3,v5)) / (card(v1)*card(v2)*card(v3)*card(v5));
coefsVr('1',d2,d3,d5) = coefsVr('1',d2,d3,d5)/2;
coefsVr(d1,'1',d3,d5) = coefsVr(d1,'1',d3,d5)/2;
coefsVr(d1,d2,'1',d5) = coefsVr(d1,d2,'1',d5)/2;
coefsVr(d1,d2,d3,'1') = coefsVr(d1,d2,d3,'1')/2;

* Value fucntion iteration

set k1 future income nodes /1*7/;

parameter GHNodes(k1) Gaussian Hermite quadrature nodes;
GHNodes('1') =  -2.651961356835233*sqrt(2.0);
GHNodes('2') =  -1.673551628767471*sqrt(2.0);
GHNodes('3') =  -0.8162878828589647*sqrt(2.0);
GHNodes('4') =  0;
GHNodes('5') =  0.8162878828589647*sqrt(2.0);
GHNodes('6') =  1.673551628767471*sqrt(2.0);
GHNodes('7') =  2.651961356835233*sqrt(2.0);

parameter GHWeights(k1) Gaussian Hermite quadrature weights;
GHWeights('1') =  0.0009717812450995192/sqrt(pi);
GHWeights('2') =  0.05451558281912703/sqrt(pi);
GHWeights('3') =  0.4256072526101278/sqrt(pi);
GHWeights('4') =  0.8102646175568073/sqrt(pi);
GHWeights('5') =  0.4256072526101278/sqrt(pi);
GHWeights('6') =  0.05451558281912703/sqrt(pi);
GHWeights('7') =  0.0009717812450995192/sqrt(pi);

Parameters
alpha1
M1
y1 income at t1
y1plus(k1) income at t1+1
;
alpha1 = 0.5;
M1 = 0.09;
y1 = 1;
y1plus(k1) = y1**rho1 * exp(sigma*GHNodes(k1));

Variables
Vond value function of no default house owner
Vondp(k1) Vond at t+1
EVondp  expected Vond at t+1

Vod  value function of default house owner
Vodp(k1) Vod at t1+1
EVodp    expected Vod at t1+1

Vrent value function of renters
Vrentp(k1)
EVrentp   expected Vrent at t1+1

rent1 rent
c1 consumption
s1 saving
s1plus saving at t1+1
b1 credit card balance
b1plus balance at t1+1
p1 cash payment
uo1 utility of house owner
ur1 utility of renter

arccossond
arccosbond

arccossod
arccosbod

arccossr
arccosbr

;

Equations
savond
bond
bond0
bond1
cashond
cashondup
utilond
arccossavond
arccosborond
Vondplus(k1)
Evondplus
objvond

savod
bod
bod0
bod1
cashod
cashodup
utilod
arccossavod
arccosborod
Vodplus(k1)
Evodplus
objvod

savr1
br1
br10
br11
cashr1
cashr1up
utilr1
arccossavr
arccosborr
Vrentplus(k1)
Evrentplus
objvrent

;
* value funciton if no default
savond.. s1plus =e= (1+rs)*(s1 + y1 - M1 -p1);
bond.. b1plus =e= (1+rb)*(b1 + c1 - p1);
bond1.. b1+c1-p1 =g= 0;
bond0.. b1+c1-p1 =l= Blim/(1+rb);
cashond.. s1 + y1 - M1 -p1 =g= 0;
cashondup.. (1+rs)*(s1 + y1 - M1 -p1) =l= 19.9;
utilond.. uo1 =e= (M1**alpha1 * c1**(1-alpha1))**(1-rho2)/(1-rho2);
arccossavond.. arccossond =e= arccos(dzsv*(s1plus - extsvmin)-1);
arccosborond.. arccosbond =e= arccos(dzb*(b1plus - extbmin)-1);

Vondplus(k1)..  Vondp(k1) =e=  sum( (d1,d2,d3,d4,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)+ord(d5)<=card(d1)+4), coefsVo(d1,d2,d3,d4,d5)*
                 cos((ord(d1)-1)*arccossond) *
                 cos((ord(d2)-1)*arccosbond) *
                 cos((ord(d3)-1)*arccos(dza*(alpha1 - extamin)-1))*
                 cos((ord(d4)-1)*arccos(dzmp*(M1 - extmpmin)-1))*
                 cos((ord(d5)-1)*arccos(dzy*(y1plus(k1) - extymin)-1)) );

Evondplus.. Evondp =e= sum(k1, GHWeights(k1)*Vondp(k1));
objvond.. Vond =e= uo1 + beta*Evondp;

*value function if default
savod.. s1plus =e= (1+rs)*(s1 + y1 - p1);
bod.. b1plus =e= (1+rb)*(b1 + c1 -p1);
bod0.. b1 + c1 -p1 =g= 0;
bod1.. (1+rb)*(b1 + c1 -p1) =l= Blim-0.01;
cashod.. s1 + y1 -p1 =g= 0;
cashodup.. (1+rs)*( s1 + y1 -p1 ) =l= 19.9;
utilod.. uo1 =e= (M1**alpha1 * c1**(1-alpha1))**(1-rho2)/(1-rho2);
arccossavod.. arccossod =e= arccos(dzsv*(s1plus - extsvmin)-1);
arccosborod.. arccosbod =e= arccos(dzb*(b1plus - extbmin)-1);

Vodplus(k1).. Vodp(k1) =e=  sum( (d1,d2,d3,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d5)<=card(d1)+3), coefsVr(d1,d2,d3,d5)*
                 cos((ord(d1)-1)*arccossod) *
                 cos((ord(d2)-1)*arccosbod) *
                 cos((ord(d3)-1)*arccos(dza*(alpha1 - extamin)-1))*
                 cos((ord(d5)-1)*arccos(dzy*(y1plus(k1)-extymin)-1)) );

Evodplus.. Evodp =e= sum(k1, GHWeights(k1)*Vodp(k1));
objvod..  Vod =e= uo1 + beta*Evodp;

*value function if rent
savr1..  s1plus =e= (1+rs)*(s1 + y1 - rent1 - p1);
br1.. b1plus =e= (1+rb)*(b1 + c1 - p1);
br10.. b1 + c1 -p1 =g= 0;
br11.. b1 + c1 -p1 =l= Blim/(1+rb)-0.01;
cashr1.. s1 + y1 -rent1 -p1 =g= 0;
cashr1up.. (1+rs)*( s1 + y1 -rent1 -p1 ) =l= 19.9;
utilr1.. ur1 =e= (rent1**alpha1 * c1**(1-alpha1))**(1-rho2)/(1-rho2);
arccossavr.. arccossr =e= arccos(dzsv*(s1plus - extsvmin)-1);
arccosborr.. arccosbr =e= arccos(dzb*(b1plus - extbmin)-1);

Vrentplus(k1).. Vrentp(k1) =e=  sum( (d1,d2,d3,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d5)<=card(d1)+3), coefsVr(d1,d2,d3,d5)*
                 cos((ord(d1)-1)*arccossr) *
                 cos((ord(d2)-1)*arccosbr) *
                 cos((ord(d3)-1)*arccos(dza*(alpha1  - extamin)-1))*
                 cos((ord(d5)-1)*arccos(dzy*(y1plus(k1)-extymin)-1)) );

Evrentplus.. EVrentp =e= sum(k1, GHWeights(k1)*Vrentp(k1));
objvrent.. Vrent =e= ur1 + beta*Evrentp;

* bound constraints
s1.lo = epsilon;
s1.up = 10;
b1.lo = epsilon;
b1.up = Blim;
s1plus.lo = epsilon;
s1plus.up = 19.9;
b1plus.lo = epsilon;
b1plus.up = Blim;
c1.lo = epsilon;
c1.up = 5;
rent1.lo = epsilon;
rent1.up = 2;
p1.lo = epsilon;
p1.up = 10;

model dpRenter1 "value function of renters at t1" /objvrent,Evrentplus,Vrentplus,arccossavr,arccosborr,utilr1,cashr1,cashr1up,br11,br10,br1,savr1/
      dpOnd "value function of no default houseowner at t1"/objvond,Evondplus,Vondplus,arccossavond,arccosborond,utilond,cashond,cashondup,bond0,bond1,bond,savond/
      dpOd "value function of default houseowner at t1"/objvod,Evodplus,Vodplus,arccossavod,arccosborod,utilod,cashod,cashodup,bod0,bod1,bod,savod/;

parameters
vondsol(v1,v2,v3,v4,v5)
vodsol(v1,v2,v3,v4,v5)
Vowner1(v1,v2,v3,v4,v5)
vr1sol(v1,v2,v3,v5)
vondsolstatus(v1,v2,v3,v4,v5)
D(v1,v2,v3,v4,v5)
c1ndsol(v1,v2,v3,v4,v5)
p1ndsol(v1,v2,v3,v4,v5)
bndpsol(v1,v2,v3,v4,v5)

c1dsol(v1,v2,v3,v4,v5)
p1dsol(v1,v2,v3,v4,v5)
bdpsol(v1,v2,v3,v4,v5)
;

D(v1,v2,v3,v4,v5) = 0;


Sets t1 period 1 /1/;

loop(t1,
         loop((v1,v2,v3,v4,v5),
                 s1.fx = sinit(v1);
                 b1.fx = binit(v2);
                 alpha1 = ainit(v3);
                 M1 = mpinit(v4);
                 y1 = yinit(v5);
                 y1plus(k1) = y1**rho1 * exp(sigma*GHNodes(k1));
         if (sinit(v1) + yinit(v5) < mpinit(v4) ,
                 vondsol(v1,v2,v3,v4,v5) = -999999;
                 vondsolstatus(v1,v2,v3,v4,v5) = 4;
                 D(v1,v2,v3,v4,v5) = 1;
         else
*initial guess for Vond
                 p1.l = 0.1;
                 c1.l = 0.1;
                 s1plus.l = (s1.l + y1 - M1 - p1.l);
                 b1plus.l = (b1.l + c1.l - p1.l);
                 arccossond.l = arccos(dzsv*(s1plus.l - extsvmin)-1);
                 arccosbond.l = arccos(dzb*(b1plus.l - extbmin)-1);
                 uo1.l = (M1**alpha1 * c1.l**(1-alpha1))**(1-rho2)/(1-rho2);
                 Vondp.l(k1) =  sum( (d1,d2,d3,d4,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)+ord(d5)<=card(d1)+4), coefsVo(d1,d2,d3,d4,d5)*
                                  cos((ord(d1)-1)*arccossond.l) *
                                  cos((ord(d2)-1)*arccosbond.l) *
                                  cos((ord(d3)-1)*arccos(dza*(alpha1 - extamin)-1))*
                                  cos((ord(d4)-1)*arccos(dzmp*(M1 - extmpmin)-1))*
                                  cos((ord(d5)-1)*arccos(dzy*(y1plus(k1)-extymin)-1)) );

                 Evondp.l = sum(k1, GHWeights(k1)*Vondp.l(k1));
                 Vond.l = uo1.l + beta*Evondp.l;

                 solve dpOnd using nlp maximizing Vond;
                 vondsol(v1,v2,v3,v4,v5) = objvond.l;
                 vondsolstatus(v1,v2,v3,v4,v5) = dpOnd.ModelStat;
                 D(v1,v2,v3,v4,v5) $ (vondsolstatus(v1,v2,v3,v4,v5) > 2 ) = 1;
                 c1ndsol(v1,v2,v3,v4,v5) = c1.l;
                 p1ndsol(v1,v2,v3,v4,v5) = p1.l;
                 bndpsol(v1,v2,v3,v4,v5) = b1plus.l;
         );
*initial guess for Vod
         p1.l = 0.1;
         c1.l = 0.1;
         s1plus.l = (s1.l + y1 - p1.l);
         b1plus.l = (b1.l + c1.l - p1.l);
         arccossod.l = arccos(dzsv*(s1plus.l - extsvmin)-1);
         arccosbod.l = arccos(dzb*(b1plus.l - extbmin)-1);
         uo1.l = (M1**alpha1 * c1.l**(1-alpha1))**(1-rho2)/(1-rho2);
         Vodp.l(k1) =  sum( (d1,d2,d3,d4,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)+ord(d5)<=card(d1)+3), coefsVr(d1,d2,d3,d5)*
                 cos((ord(d1)-1)*arccossod.l) *
                 cos((ord(d2)-1)*arccosbod.l) *
                 cos((ord(d3)-1)*arccos(dza*(alpha1 - extamin)-1))*
                 cos((ord(d5)-1)*arccos(dzy*(y1plus(k1)-extymin)-1)) );
         Evodp.l = sum(k1, GHWeights(k1)*Vodp.l(k1));
         Vond.l = uo1.l + beta*Evodp.l;

         solve dpOd using nlp maximizing Vod;
         vodsol(v1,v2,v3,v4,v5) = objvod.l;
         c1dsol(v1,v2,v3,v4,v5) = c1.l;
         p1dsol(v1,v2,v3,v4,v5) = p1.l;
         bdpsol(v1,v2,v3,v4,v5) = b1plus.l;

         Vowner1(v1,v2,v3,v4,v5) = max(vondsol(v1,v2,v3,v4,v5),vodsol(v1,v2,v3,v4,v5));
         Vowner1(v1,v2,v3,v4,v5)$(vondsolstatus(v1,v2,v3,v4,v5) > 2) = vodsol(v1,v2,v3,v4,v5);
*initial guess for Vrent
         rent1.l = 0.1;
         c1.l = 0.1;
         p1.l = 0.1;
         s1plus.l = (s1.l + y1 - rent1.l - p1.l);
         b1plus.l = (b1.l + c1.l - p1.l);
         arccossr.l = arccos(dzsv*(s1plus.l - extsvmin)-1);
         arccosbr.l = arccos(dzb*(b1plus.l - extbmin)-1);
         ur1.l = (rent1.l**alpha1 * c1.l**(1-alpha1))**(1-rho2)/(1-rho2);
         Vrentp.l(k1) =  sum( (d1,d2,d3,d4,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)+ord(d5)<=card(d1)+3), coefsVr(d1,d2,d3,d5)*
                 cos((ord(d1)-1)*arccossr.l) *
                 cos((ord(d2)-1)*arccosbr.l) *
                 cos((ord(d3)-1)*arccos(dza*( alpha1  - extamin)-1))*
                 cos((ord(d5)-1)*arccos(dzy*(y1plus(k1)-extymin)-1)) );
         EVrentp.l  = sum(k1, GHWeights(k1)*Vrentp.l(k1));
         Vrent.l = ur1.l + beta*Evrentp.l;

         solve dpRenter1 using nlp maximizing Vrent;
         vr1sol(v1,v2,v3,v5) = objvrent.l;
         );

         coefsVo(d1,d2,d3,d4,d5) = 0;
         coefsVo(d1,d2,d3,d4,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)+ord(d5)<=card(d1)+4) =
    32*sum((v1,v2,v3,v4,v5), TTs(v1,d1)*TTb(v2,d2)*TTa(v3,d3)*TTmp(v4,d4)*TTy(v5,d5)*Vowner1(v1,v2,v3,v4,v5)) / (card(v1)*card(v2)*card(v3)*card(v4)*card(v5));

         coefsVo('1',d2,d3,d4,d5) = coefsVo('1',d2,d3,d4,d5)/2;
         coefsVo(d1,'1',d3,d4,d5) = coefsVo(d1,'1',d3,d4,d5)/2;
         coefsVo(d1,d2,'1',d4,d5) = coefsVo(d1,d2,'1',d4,d5)/2;
         coefsVo(d1,d2,d3,'1',d5) = coefsVo(d1,d2,d3,'1',d5)/2;
         coefsVo(d1,d2,d3,d4,'1') = coefsVo(d1,d2,d3,d4,'1')/2;

         coefsVr(d1,d2,d3,d5) = 0;
         coefsVr(d1,d2,d3,d5)$(ord(d1)+ord(d2)+ord(d3)+ord(d5)<=card(d1)+3) =
    16*sum((v1,v2,v3,v5), TTs(v1,d1)*TTb(v2,d2)*TTa(v3,d3)*TTy(v5,d5)*vr1sol(v1,v2,v3,v5)) / (card(v1)*card(v2)*card(v3)*card(v5));
         coefsVr('1',d2,d3,d5) = coefsVr('1',d2,d3,d5)/2;
         coefsVr(d1,'1',d3,d5) = coefsVr(d1,'1',d3,d5)/2;
         coefsVr(d1,d2,'1',d5) = coefsVr(d1,d2,'1',d5)/2;
         coefsVr(d1,d2,d3,'1') = coefsVr(d1,d2,d3,'1')/2;
);


File modeloutput2 /sol_path.csv/;
modeloutput2.nw=18;
modeloutput2.nr=2;
modeloutput2.nz=1e-15;

Put modeloutput2;
Put "s, b, alpha, mp, y,Value_fun, defualt, cnd, pnd, bndplus, cd, pd, bdplus" /;

modeloutput2.pc=5;
*modeloutput3.pw=4000;

loop((v1,v2,v3,v4,v5),
  put  sinit(v1):14:6;
  put  binit(v2):14:6;
  put  ainit(v3):14:6;
  put  mpinit(v4):14:6;
  put  yinit(v5):14:6;
  put  Vowner1(v1,v2,v3,v4,v5);
  put  D(v1,v2,v3,v4,v5):14:6;
  put  c1ndsol(v1,v2,v3,v4,v5):14:6;
  put  p1ndsol(v1,v2,v3,v4,v5):14:6;
  put  bndpsol(v1,v2,v3,v4,v5):14:6;
  put  c1dsol(v1,v2,v3,v4,v5):14:6;
  put  p1dsol(v1,v2,v3,v4,v5):14:6;
  put  bdpsol(v1,v2,v3,v4,v5):14:6;
  put /;
);







