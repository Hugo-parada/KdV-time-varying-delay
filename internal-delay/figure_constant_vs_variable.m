% Internal delay: tau(t)=1+0.5*sin(t) versus constants 1.5 and 0.5.
clear; clc;
[t1,E1] = kdvint_timedelay(1.5);
[t2,E2] = kdvint_ctedelay(1.5);
[t3,E3] = kdvint_ctedelay(0.5);
figure;
plot(t1,log(E1),'LineWidth',1.2); hold on;
plot(t2,log(E2),'LineWidth',1.2);
plot(t3,log(E3),'LineWidth',1.2); hold off;
legend('\tau(t)','\tau_{max}=1.5','\tau_{min}=0.5','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
