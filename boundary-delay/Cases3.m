% Variable delay versus the corresponding minimum and maximum constants
clear;clc;
alpha=0.2;
beta=0.2;
d=0.2;

[t1,E1]=kdv_time_delay(alpha,beta,d);
[t2,E2]=kdv_constant_delay(alpha,beta,2.5*d);
[t3,E3]=kdv_constant_delay(alpha,beta,0.5*d);

plot(t1,log(E1),t2,log(E2),t3,log(E3),'LineWidth',1.2)
legend('\tau(t)','\tau_{max}=2.5d','\tau_{min}=0.5d')
xlabel('time');ylabel('ln(E(t))');grid on
