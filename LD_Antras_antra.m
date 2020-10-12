%Matematinio modeliavimo technologijos. LD2 Dainius Varna EKSfm20 20200033
%2 Dalis
%Nr.1 y = a*x*(e^(bx)) ; 'istiesinta funkcija' y'=a'x+b' ; pakeitimas
%y'=ln(y/x), a=e^b', b=a'
%duomenu failas data2.mat
function LD_Antras_antra()
close all;
clc;
%data2.mat duomenys:
yd2 = [0.0323064003719004 0.231348047553774 0.121322571150384 0.257228407554052 0.259299860639396 0.344122540104969 0.121794277221822 0.166236450167924 0.292890761402347 0.157926674919028 0.216361712214747 0.288911308729504 0.282933987298767 0.403248733361967 0.272092901548121 0.409378184530305 0.315305462730464 0.295440257560431 0.575755185740552 0.418293378543285 0.430347528433920 0.403804415892895 0.445010203391100 0.611424322179399 0.455278964702241 0.459386191296984 0.530964948046310 0.435973945122411 0.649907608155333 0.576623536826968 0.892451831766217 0.776308894043234 0.618649005539376 0.949797361638746 0.779914958895666 0.742945250540432 0.761958188217124 0.881623636359153 0.893510835103563 1.02282626308904 0.975707715314429 1.01676631931371 0.915329390636685 0.994584079923118 1.06927026510411 1.31702676831184 1.01634010012880 1.13115467714731 1.23709466041136 1.29643380499495 1.33057002285573 1.27711780848243 1.71144329323348 1.61419794840596 1.36375840650238 1.60232949161230 1.84777783361976 1.77389484907317 1.81494032849274 1.89134371343091 1.88525280210668 1.86231379712229 2.09893495570363 2.32069350152099 2.00987483225021 2.25307181623569 2.35165255302776 2.46644109464859 2.45148431222523 2.58667248318290 2.80704765114765 2.69121893591594 2.62045117334415 2.67311636906179 2.85379295908643 3.08853424060475 3.18091246791208 3.28963889391109 3.26982607485719 3.24655790420656 3.61635245895249 3.64007914581112 3.67856918538396 3.93039621853483 4.01744729980357 4.01453015246790 4.24947235681183 4.34332998694981 4.44312851929957 4.36068431520984 4.60714963012858 4.78722910959019 5.10096386819180 5.14413074429537 5.30459256412172 5.36190784977724 5.46625122425692 5.81482080739339 5.98647261952859 6.02412740028718];
x = [1 1.04040404040404 1.08080808080808 1.12121212121212 1.16161616161616 1.20202020202020 1.24242424242424 1.28282828282828 1.32323232323232 1.36363636363636 1.40404040404040 1.44444444444444 1.48484848484848 1.52525252525253 1.56565656565657 1.60606060606061 1.64646464646465 1.68686868686869 1.72727272727273 1.76767676767677 1.80808080808081 1.84848484848485 1.88888888888889 1.92929292929293 1.96969696969697 2.01010101010101 2.05050505050505 2.09090909090909 2.13131313131313 2.17171717171717 2.21212121212121 2.25252525252525 2.29292929292929 2.33333333333333 2.37373737373737 2.41414141414141 2.45454545454545 2.49494949494950 2.53535353535354 2.57575757575758 2.61616161616162 2.65656565656566 2.69696969696970 2.73737373737374 2.77777777777778 2.81818181818182 2.85858585858586 2.89898989898990 2.93939393939394 2.97979797979798 3.02020202020202 3.06060606060606 3.10101010101010 3.14141414141414 3.18181818181818 3.22222222222222 3.26262626262626 3.30303030303030 3.34343434343434 3.38383838383838 3.42424242424242 3.46464646464646 3.50505050505051 3.54545454545455 3.58585858585859 3.62626262626263 3.66666666666667 3.70707070707071 3.74747474747475 3.78787878787879 3.82828282828283 3.86868686868687 3.90909090909091 3.94949494949495 3.98989898989899 4.03030303030303 4.07070707070707 4.11111111111111 4.15151515151515 4.19191919191919 4.23232323232323 4.27272727272727 4.31313131313131 4.35353535353535 4.39393939393939 4.43434343434343 4.47474747474748 4.51515151515152 4.55555555555556 4.59595959595960 4.63636363636364 4.67676767676768 4.71717171717172 4.75757575757576 4.79797979797980 4.83838383838384 4.87878787878788 4.91919191919192 4.95959595959596 5];

aa=[4;0.5];
y1=myfun1(aa,x); 
a1 = lsqcurvefit(@myfun1,[0;0],x,y1)%randam parametrus tiesiogiai

A=[ones(size(x))' x'];
aa1=A\log(y1)'; %randam parametrus istiesinus lygti...
aa1(1)=exp(aa1(1))
fprintf('%f \n',aa1);

%y2=myfun2(aa,x)+randn(size(x));
y2=myfun2(aa,x);
a2 = lsqcurvefit(@myfun2,[0;0],x,y2)%randam parametrus tiesiogiai
A=[log(x)' ones(size(x))'];%randam parametrus sprendziant TLS...
aa2=A\y2' %randam parametrus istiesinus lygti...


subplot(121);
plot(x,y1,'*',x,myfun1(a1,x),x,myfun1(aa1,x));grid
legend('duomenys','lsqcurvefit','istiesinus');
title('Aproksimacija funkcija {\it y = ax}e^{{\it bx}}');

subplot(122);
plot(x,y2,'*',x,myfun2(a2,x),x,myfun2(aa2,x));grid
legend('duomenys','lsqcurvefit','istiesinus');
title('Aproksimacija funkcija {\it y = a }ln{\it x + b}');

function y = myfun1(a,x)
y = a(1)*exp(a(2)*x);

function y = myfun2(a,x)
y = a(1)'*log(x)+a(2)';
