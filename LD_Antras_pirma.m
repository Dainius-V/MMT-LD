%Matematinio modeliavimo technologijos. LD2 Dainius Varna EKSfm20 20200033
%1 Dalis
function LD_Antras_pirma()
clear all
close all
a=mod(20200033,6);
fprintf('mod(20200033) = %d', a);
N=21 %galima varijuoti N, didinant paklaida dideja
x=linspace(-10,10,N); %Sugeneruojame -10 <= x <= 10 skaicius 
yN=f1(x) %Apskaiciuojame lygti pagal sugeneruotus skaicius

%Aproksimacijos:
l=lagranp(x,yN); %Lagaranþo
n=newtonp(x,yN); %Niutono
c=cheby('f1',N,min(x),max(x)) %Èebyðevo
xo=-10;10;
Mp = 3; Np = 2;
p=padeap('f1',xo,Mp,Np) %Pade


xx=linspace(-10,10,100);
yl = polyval(l,xx); %Lagaranþo
y2 = polyval(n,xx); %Niutono
y3 = polyval(c,xx); %Èebyðevo
y4 = polyval(p,xx); %Pade
y5 = spline(x,yN,xx);%Splainai


%Atvaizdavimas
subplot(3,2,1)
plot(xx,y3,'b-',x,yN,'*', xx, f1(xx));grid
title('Cebysevo daugianaris');

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

%subplot(3,2,6)
%plot(err1,err2);grid
%title('MSE');

%MSE (vidutinës kvadratinës klaidos) 
ye=f1(xx)
err1 = immse(ye,yl);
err2 = immse(ye,y2);
err3 = immse(ye,y3);
err4 = immse(ye,y4);
err5 = immse(ye,y5);

fprintf('err1 - %f\n\n', err1);
fprintf('err2 - %f\n\n', err2);
fprintf('err3 - %f\n\n', err3);
fprintf('err4 - %f\n\n', err4);
fprintf('err5 - %f\n\n', err5);


end

%Gauta funkcija
function y=f1(x)

y=1./(1+exp(x))
end

%Lagaranþo
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

%Èebyðevo:
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