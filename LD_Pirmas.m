%Matematinio modeliavimo technologijos. LD1 Dainius Varna EKSfm20 20200033
function LD_Pirmas()
clear
Af=[55 58 63 65 66 59 56 87 90 85 81];% temperatûros Farenheitais
bC=[13 14 17 18 19 15 13 30 32 29 27];% temperatûros Celsijais
NA=length(2); % 2, nes reikia surasti k1, k2
x=zeros(NA,1);P=1000000*eye(NA,NA);
xx=[]
for k=1:length(Af) % Paleidziame cikla 11 k. (pagal temperaturu masyvus)
    A(k,:) = [Af(k) Af(11)];
    b(k,:) = A(k,:)*((Af(k)'*Af(k))^(-1)*Af(k)'*bC(k))+x(1);
    [x,K,P]=rlse_online(A(k,:),b(k,:),x,P);
end
k1=x(1); %k1 pasirenkame is gauto rlse algoritmo
k2=bC(1)-Af(1)*k1; %apskaiciuojame k2
fprintf('k1 ivertis taikant RLSE = %f\n',k1);
fprintf('k2 ivertis taikant RLSE = %f\n',k2);
x=A\b; %Palyginame A ir b pagal TLS
fprintf('Palyginimui A\\b = %f\n%f\n',x(1),k2);
for i=1:length(Af) %Ciklas apskaiciuoti Celsius naudojant gautus k1, k2
    b_x(i)=[Af(i)*k1+k2];
end
fprintf('\n\n');
%Nupiesiame linijinigrafika:
plot(b_x')%RLSE
hold on
plot(bC')%Tyrimo
grid;

function[x,K,P]=rlse_online(aT_k1,b_k1,x,P) %RLSE algoritmas
%aT_k1 - Farenheitai, b_k1 - Celsijai, x - iverciai
K=P*aT_k1'/(aT_k1*P*aT_k1'+1);%Stiprinimo matricos sandauga su paklaida
x=x+K*(b_k1-aT_k1*x);%Gaunant naujus duomenis, atnaujiname parametrus
P=P-K*aT_k1*P;%Inversines matricos atnaujinimas, kai gaunami nauji duomenys

