% Influence of the size of the time-varying delay
clear;clc;
alpha=0.1;
beta=0.1; 

[t1,E1]=kdv_time_delay(alpha,beta,0.1);
[t2,E2]=kdv_time_delay(alpha,beta,0.7);
[t3,E3]=kdv_time_delay(alpha,beta,1.3);

plot(t1,log(E1),t2,log(E2),t3,log(E3),'LineWidth',1.2)
legend('d=0.1','d=0.7','d=1.3')
xlabel('time');ylabel('ln(E(t))');grid on
