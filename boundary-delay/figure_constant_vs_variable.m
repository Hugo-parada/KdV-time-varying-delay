% Boundary delay: variable delay versus its maximum and minimum constants.
clear; clc;
alpha = 0.2; beta = 0.2; d = 0.2;
[t1,E1] = kdv_time_delay(alpha,beta,d);
[t2,E2] = kdv_cte_delay(alpha,beta,2.5*d);
[t3,E3] = kdv_cte_delay(alpha,beta,0.5*d);
figure;
plot(t1,log(E1),'LineWidth',1.2); hold on;
plot(t2,log(E2),'LineWidth',1.2);
plot(t3,log(E3),'LineWidth',1.2); hold off;
legend('\tau(t)','\tau_{max}=2.5d','\tau_{min}=0.5d','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
