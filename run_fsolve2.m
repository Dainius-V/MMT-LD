% I = I0 ( exp(Ue/kT)-1 )
% k = 1.38*10^-23; e = 1.602*10^-19;
% U(V) 0.1 I(A) 0.8139*10^-4
% U(V) 0.12 I(A) 1.9808*10^-4

function run_fsolve2()
clc;
x0 = [4.8007*10^-4 1.9808*10^-4];
options=optimset('Display','iter'); 
x = fsolve(@f46,x0,options)
%x = fsolve(@(x) Teth3(x, P, H, L), x0);
disp('Sprendimas rastas')


function y = f46(x)
e=1.602*10^-19;
k=1.38*10^-23;
y(1) = x(1)*(exp((0.14*e)/(k*x(2)))-1);%f-ija kai U=0.14
y(2) = x(1)*(exp((0.12*e)/(k*x(2)))-1);%f-ija kai U=0.12