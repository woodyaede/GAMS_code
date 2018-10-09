scalars
beta discount rate /0.99833/
rho1 persistence parameter for income growth HSV /0.975/
sigma standard deviation of income shock /0.01/
gamma intratemporal substitutability between housing and non-durable consumption goods /-3/
rho2 the degree of curvature of the utility funciton /3/
age household age /35/
cbss consumption of non-durable goods in the steady state
cpss cash payment of non-durable goods in the steady state
rs saving interest rate /0.00033/
avy permanent average monthly income /1/
rb credit card debt interest rate /0.05/
Blim credit limit of all credit cards /9/
epsilon A small numbner /0.001/
myear mortgage length /140/
;

* Step 1: Define the transformed deterministic problem

Sets t time /1*301/;

parameters
y(t)            income path
betas(t)        discounts
ch(t)           housing expenditure
alpha           housing cost share
h               housing service
;

betas(t) = beta**(ord(t)-1);
y(t) = 1;



variables
obj objective creterion
u(t) utility
w(t) saving
b(t) credit debt
;

positive variables
cb(t) non-durable goods consumption
cp(t) cash payment of non-durable goods
;

Equations
obj1 Objective function
util(t) utility function
sav(t) process of saving growth
bor(t) process of credit card debt grwoth
debt0(t) lower bound of credit card debt
debt1(t) upper bound of credit card debt
cash(t)  saving is nonnegative
;

sav(t)$(ord(t) lt card(t)) .. w(t+1) =e= (1+rs)*(w(t)+y(t)-ch(t)-cp(t));
bor(t)$(ord(t) lt card(t)) .. b(t+1) =e= (1+rb)*(b(t)+cb(t)-cp(t));
debt0(t) .. b(t)+cb(t)-cp(t) =g= 0;
debt1(t) .. b(t)+cb(t)-cp(t) =l= Blim;
cash(t) .. w(t)+y(t)-ch(t)-cp(t) =g= 0;
util(t)$(ord(t) lt card(t)) .. u(t) =e= (h**alpha * cb(t)**(1-alpha))**rho2/(1-rho2);
obj1 .. obj =e= sum(t$(ord(t) lt card(t)), betas(t)*u(t)) +
  sum(t$(ord(t)=card(t)), betas(t)/(1-beta)*(h**alpha * (y(t)+ rs*w(t)/(1+rs)- rb*b(t)/(1+rb))**(1-alpha))**rho2/(1-rho2));

* bound constraints
w.lo(t) = epsilon;
w.up(t) = 20;
b.lo(t) = epsilon;
b.up(t) = Blim;
cp.lo(t) = epsilon;
cb.lo(t) = epsilon;
cb.up(t) = 10;


* initial guess
w.l(t) = epsilon;
b.l(t) = epsilon;
cb.l(t) = y(t)+ rs*w.l(t)/(1+rs)- rb*b.l(t)/(1+rb);
cp.l(t) = cb.l(t);
alpha = 0.5;
h = 0.5;


u.l(t) = (h**alpha * cb.l(t)**(1-alpha))**rho2/(1-rho2);
obj.l = sum(t$(ord(t) lt card(t)), betas(t)*u.l(t)) +
  sum(t$(ord(t)=card(t)), betas(t)/(1-beta)*(h**alpha * cb.l(t)**(1-alpha))**rho2/(1-rho2));

Option limrow=0,limcol=0,solprint=off;
OPTION ITERLIM = 500000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION NLP= conopt;
option decimals = 7;

model CCDP /all/;

*******************************************************
* Step 2: Optimization step

set v1 grid space/1*10/
set v3 grid space/1*3/
alias(v1,v2,d1,d2)
alias(v3,v4,d3,d4);

*approximation domain

scalars
wmin lower bound of w /0/
wmax upper bound of w /20/
bmin lower bound of b /0/
bmax upper bound of b /9/
amin lower bound of alpha /0.1/
amax upper bound of alpha /0.9/
mpmin lower bound of mortgage payment /0.05/
mpmax upper bound of mortgage payment /0.9/
;

*Chebyshev parameters

parameters zw(v1)
         zb(v2)
         za(v3)
         zmp(v4)
         extw
         extb
         exta
         extmp
         extwmin
         extwmax
         extbmin
         extbmax
         extamin
         extamax
         extmpmin
         extmpmax
         dzw
         dzb
         dza
         dzmp
;

zw(v1) = - cos((2*ord(v1)-1)/(2*card(v1))*pi);
zb(v2) = - cos((2*ord(v2)-1)/(2*card(v2))*pi);
za(v3) = - cos((2*ord(v3)-1)/(2*card(v3))*pi);
zmp(v4) = - cos((2*ord(v4)-1)/(2*card(v4))*pi);

extw = -(zw('1')+1)*(wmax-wmin)/2/zw('1');
extb = -(zb('1')+1)*(bmax-bmin)/2/zb('1');
exta = -(za('1')+1)*(amax-amin)/2/za('1');
extmp = -(zmp('1')+1)*(mpmax-mpmin)/2/zmp('1');

extwmin = wmin - extw;
extbmin = bmin - extb;
extamin = amin - exta;
extmpmin = mpmin - extmp;

extwmax = wmax + extw;
extbmax = bmax + extb;
extamax = amax + exta;
extmpmax = mpmax - extmp;

dzw = 2/(extwmax-extwmin);
dzb = 2/(extbmax-extbmin);
dza = 2/(extamax - extamin);
dzmp = 2/(extmpmax - extmpmin);

parameters
winit(v1) Grid Saving
binit(v2) Grid Credit Card Debt
ainit(v3) Grid preference coefficient alpha
mpinit(v4) Grid of mortgage payment
;

winit(v1) = extwmin + (1+zw(v1))/dzw;
binit(v2) = extbmin + (1+zb(v2))/dzb;
ainit(v3) = extamin + (1+za(v3))/dza;
mpinit(v4) = extmpmin + (1+zmp(v4))/dzmp;

Parameters
cbsol(v1,v2,v3,v4) solution for consumption of non-durable good
cpsol(v1,v2,v3,v4) solution for cash payment of non-durable good
lambdas(v1,v2,v3,v4) multiplier for lower bound of saving
lambdab0(v1,v2,v3,v4) multiplier for lower bound of credit card debt
lambdab1(v1,v2,v3,v4) multiplier for upper bound of credit card debt

w2sol(v1,v2,v3,v4)
b2sol(v1,v2,v3,v4)
solstatus(v1,v2,v3,v4) model termination message
;

loop((v1,v2,v3,v4),
  w.fx('1') = winit(v1);
  b.fx('1') = binit(v2);
  alpha = ainit(v3);
  h = mpinit(v4);
  loop (t$(ord(t) le myear ),
    ch(t) = h;
  );
  loop (t$(ord(t) ge myear ),
    ch(t) = 0;
  );

  solve CCDP using nlp maximizing obj;

  cbsol(v1,v2,v3,v4) = cb.l('1');
  cpsol(v1,v2,v3,v4) = cp.l('1');
  lambdas(v1,v2,v3,v4) = abs(cash.m('1'));
  lambdab0(v1,v2,v3,v4) = abs(debt0.m('1'));
  lambdab1(v1,v2,v3,v4) = abs(debt1.m('1'));
  w2sol(v1,v2,v3,v4) = w.l('2');
  b2sol(v1,v2,v3,v4) = b.l('2');
  solstatus(v1,v2,v3,v4) = CCDP.ModelStat;
);

File modeloutput2 /sol_cash_w2b2_path.csv/;
modeloutput2.nw=18;
modeloutput2.nr=2;
modeloutput2.nz=1e-15;

Put modeloutput2;
Put "Saving, Borrowing, alpha, ch, s1, b1, solver_status" /;

modeloutput2.pc=5;
*modeloutput3.pw=4000;

loop((v1,v2,v3,v4),
  put  winit(v1):14:6;
  put  binit(v2):14:6;
  put  ainit(v3):14:6;
  put  mpinit(v4):14:6;
  put w2sol(v1,v2,v3,v4):14:6;
  put b2sol(v1,v2,v3,v4):14:6;
  put  solstatus(v1,v2,v3,v4):14:6;
  put /;
);

************************************************
* Step 3: Approximation step

