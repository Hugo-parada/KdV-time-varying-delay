% Variable delay tau(t)=1+0.5 sin(t) versus its extrema
clear;clc;

[t1,E1]=kdvint_timedelay(1.5);
[t2,E2]=kdvint_constant_delay(1.5);
[t3,E3]=kdvint_constant_delay(0.5);

plot(t1,log(E1),t2,log(E2),t3,log(E3),'LineWidth',1.2)
legend('\tau(t)','\tau_{max}=1.5','\tau_{min}=0.5')
xlabel('time');ylabel('ln(E(t))');grid on
