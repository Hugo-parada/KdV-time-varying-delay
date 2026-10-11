% Boundary delay: influence of the parameter d.
clear; clc;
alpha = 0.1;
beta = 0.1;
d_values = [0.1,0.7,1.3];
figure; hold on;
for j = 1:numel(d_values)
    [t,E] = kdv_time_delay(alpha,beta,d_values(j));
    plot(t,log(E),'LineWidth',1.2);
end
hold off;
legend('d=0.1','d=0.7','d=1.3','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
% d=1.3 deliberately violates sup_t dot(tau)(t)<1.
