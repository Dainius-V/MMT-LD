%Matematinio modeliavimo technologijos. LD2 Dainius Varna EKSfm20 20200033
%1 Dalis
function LD_Antras_pirma()
clear all
close all
a=mod(20200033,6); %varianto parinkimas
fprintf('mod(20200033) = %d', a);
%Error skaiciavimam masyvai:
le=[]
ne=[]
ce=[]
pe=[]
spe=[]
%Max error skaiciavimam masyvai:
lme=[]
nme=[]
cme=[]
pme=[]
spme=[]
N=[2 3 5 7 9] %Aproksimacijos eiles
for i=1:5 %Ciklas prasukamas 5 kartus, nes lagaranzo, niutono, cebysevo, pade, spline
    
    x=linspace(-10,10,N(i)); %Sugeneruojame -10 <= x <= 10 skaicius 
    yN=f1(x) %Apskaiciuojame lygti pagal sugeneruotus skaicius 
    %f1 - lygties funkcija

    %Aproksimacijos:
    l=lagranp(x,yN); %Lagaran?o
    n=newtonp(x,yN); %Niutono
    c=cheby('f1',N(i),min(x),max(x)) %?eby?evo
    xo=-10;10;
    Mp = 3; Np = 2;
    p=padeap('f1',xo,Mp,Np) %Pade
    
    xx=linspace(-10,10,100);
    yl = polyval(l,xx); %Lagaran?o
    y2 = polyval(n,xx); %Niutono
    y3 = polyval(c,xx); %?eby?evo
    y4 = polyval(p,xx); %Pade
    y5 = spline(x,yN,xx);%Splainai
    ymax = f1(xx)
    
    %Max error
    lme(i)= max(ymax-yl);
    nme(i)= max(ymax-y2);
    cme(i)= max(ymax-y3);
    spme(i)= max(ymax-y5);
    
    
    %MSE
    le(i)=immse(xx,yl);
    ne(i)=immse(xx,y2);
    ce(i)=immse(xx,y3);
    pe(i)=immse(xx,y4);
    spe(i)=immse(xx,y5);
    
    
    %Atvaizdavimas 5 grafikai (total 25)
    %geriau atvaizduoti ant vieno be subplot
    subplot(3,2,1)
    plot(xx,y3,'b-',x,yN,'*', xx, f1(xx));grid
    title('Cebysevo daugianaris ');

    subplot(3,2,2)
    plot(xx,y2,'r-',x,yN,'*', xx, f1(xx));grid
    title('Niutono daugianaris');

    subplot(3,2,3)
    plot(xx,yl,'g-',x,yN,'*', xx, f1(xx));grid
    title('Lagaranzo daugianaris');

    subplot(3,2,4)
    plot(xx,y4,'g-',x,yN,'*', xx, f1(xx));grid
    title('Pade daugianaris');

    subplot(3,2,5)
    plot(xx,y5);grid
    title('Spline daugianaris');
    fprintf('N = %d\n',N(i));
    pause % padarom pauze, paspaudus enter pasikeis aproksimacijos eile
    close all
end
fprintf('finished\n');

xxe=linspace(-10,10,100);
%Visu aproksimaciju klaidu duomenims issaugoti masyvai:
errl=[]
errn=[]
errc=[]
errp=[]
errs=[]
errxmax=[]

%paleidziame cikla uzpildyti klaidu kintamuosius (MSE)
for ie=1:5
    fprintf('%d\n\',N(ie));
    errl(ie)=le(ie);errn(ie)=ne(ie);errc(ie)=ce(ie);errp(ie)=pe(ie);
    errs(ie)=spe(ie);

end

%paleidziame cikla uzpildyti klaidu kintamuosius (MAX)
errml=[]
errmn=[]
errmc=[]
errms=[]
for ie=1:5
    fprintf('%d\n\',N(ie));
    errml(ie)=lme(i);errmn(ie)=nme(ie);errmc(ie)=cme(ie);
    errms(ie)=spme(ie);

end

fprintf('arrays filled \n\');
pause % Pauze, kad zinoti kada uzsipildys masyvai

Ne=[2 3 5 7 9] % aproksimacijos eiles klaidu palyginimo grafikui
%Atvaizdavimas:
%MSE
plot(Ne,errl,Ne,errn,Ne,errc,Ne,errs)
grid;
fprintf('errors shown \n\');
pause
close all %uzdarom grafika
%MAX
plot(Ne,errml,Ne,errmn,Ne,errmc,Ne,errms)
grid;
fprintf('max shown \n\');
pause

end

%Gauta funkcija 1 varianto
function y=f1(x)

y=1./(1+exp(x))
end

%Lagaran?o
function [l,L]=lagranp(x,y)

%Input : x = [x0 x1 ... xN], y = [y0 y1 ... yN]
%Output: l = N eiles Lagranzo daugianaris
% L = Lagranzo koeficientu daugianaris
N = length(x)-1; %daugianario eile
l = 0;
for m = 1:N + 1
    P = 1;
    for k = 1:N + 1
        if k ~= m, P = conv(P,[1 -x(k)])/(x(m)-x(k)); end
        %MATLAB dirba su daugianariu koeficientu vektoriumi (mazejancia
        %tvarka), o dvieju daugianariu daugyba atitinka koeficientu
        %vektoriu konvoliucijai...
    end
    L(m,:) = P; %Lagranzo koeficientu daugianaris
    l = l + y(m)*P; %Lagranzo daugianaris (3 lytis)
end
end

%Niutono:
function [n,DD] = newtonp(x,y)
%Output: n = Newton polynomial coefficients of degree N
N = length(x)-1;
DD = zeros(N + 1,N + 1);
DD(1:N + 1,1) = y';
for k = 2:N + 1
    for m = 1: N + 2 - k %Divided Difference Table
        DD(m,k) = (DD(m + 1,k - 1) - DD(m,k - 1))/(x(m + k - 1)- x(m));
    end
end
a = DD(1,:); %Eq.(3.2.6)
n = a(N+1); %Begin with Eq.(3.2.7)
for k = N:-1:1 %Eq.(3.2.7)
    n = [n a(k)] - [0 n*x(k)]; %n(x)*(x - x(k - 1))+a_k - 1
end
end

%?eby?evo:
function [c,x,y] = cheby(f,N,a,b)
%Input : f = function name on [a,b]
%Output: c = Newton polynomial coefficients of degree N
% (x,y) = Chebyshev nodes
if nargin == 2, a = -1; b = 1; end
k = [0: N];
theta = (2*N + 1 - 2*k)*pi/(2*N + 2);
xn = cos(theta); %Eq.(3.3.1a)
x = (b - a)/2*xn +(a + b)/2; %Eq.(3.3.1b)
y = feval(f,x);
d(1) = y*ones(N + 1,1)/(N+1);
for m = 2: N + 1
    cos_mth = cos((m-1)*theta);
    d(m) = y*cos_mth'*2/(N + 1); %Eq.(3.3.6b)
end
xn = [2 -(a + b)]/(b - a); %the inverse of (3.3.1b)
T_0 = 1; T_1 = xn; %Eq.(3.3.3b)
c = d(1)*[0 T_0] +d(2)*T_1; %Eq.(3.3.5)
for m = 3: N + 1
    tmp = T_1;
    T_1 = 2*conv(xn,T_1) -[0 0 T_0]; %Eq.(3.3.3a)
    T_0 = tmp;
    c = [0 c] + d(m)*T_1; %Eq.(3.3.5)
end
end

%Pade
function [num,den,t] = padeap(f,xo,M,N)
%Input : f = function to be approximated around xo in [xmin, xmax]
%Output: num = numerator coeffs of Pade approximation of degree M
% den = denominator coeffs of Pade approximation of degree N
%Tylor series coefficients
a(1) = feval(f,xo);
h = .01; tmp = 1;
for i = 1:M + N
    tmp = tmp*i*h; %i!h^i
    dix = difapx(i,[-i i])*feval(f,xo+[-i:i]*h)'; %derivative(Section 5.3)
    a(i + 1) = dix/tmp; %Taylor series coefficient
end
for m = 1:N
    n = 1:N; A(m,n) = a(M + 1 + m - n);
    b(m) = -a(M + 1 + m);
end
d = A\b'; %Eq.(3.4.4b)
for m = 1: M + 1
    mm = min(m - 1,N);
    q(m) = a(m:-1:m - mm)*[1; d(1:mm)]; %Eq.(3.4.4a)
end
num = q(M + 1:-1:1)/d(N); den = [d(N:-1:1)' 1]/d(N); %descending order

t=a(M + N + 1:-1:1);
end

function [c,err,eoh,A,b] = difapx(N,points)
%difapx.m to get the difference approximation for the Nth derivative
l = max(points);
L = abs(points(1)-points(2))+ 1;
if L < N + 1, error('More points are needed!'); end
for n = 1: L
    A(1,n) = 1;
    for m = 2:L + 2, A(m,n) = A(m - 1,n)*l/(m - 1); end %Eq.(5.3.5)
    l = l-1;
end
b = zeros(L,1); b(N + 1) = 1;
c =(A(1:L,:)\b)'; %coefficients of difference approximation formula
err = A(L + 1,:)*c'; eoh = L-N; %coefficient & order of error term
if abs(err) < eps, err = A(L + 2,:)*c'; eoh = L - N + 1; end
if points(1) < points(2), c = fliplr(c); end
end