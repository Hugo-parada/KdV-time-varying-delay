% Influence of the maximum delay M
clear;clc;

[t1,E1]=kdvint_timedelay(1.5);
[t2,E2]=kdvint_timedelay(2);
[t3,E3]=kdvint_timedelay(2.5);

plot(t1,log(E1),t2,log(E2),t3,log(E3),'LineWidth',1.2)
legend('M=1.5','M=2','M=2.5')
xlabel('time');ylabel('ln(E(t))');grid on
