% Influence of beta
clear;clc;
alpha=0.2;
d=0.2;

[t1,E1]=kdv_time_delay(alpha,0.2,d);
[t2,E2]=kdv_time_delay(alpha,0.4,d); 
[t3,E3]=kdv_time_delay(alpha,0.5,d);

plot(t1,log(E1),t2,log(E2),t3,log(E3),'LineWidth',1.2)
legend('\beta=0.2','\beta=0.4','\beta=0.5')
xlabel('time');ylabel('ln(E(t))');grid on
