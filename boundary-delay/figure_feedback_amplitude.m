% Additional recovered experiment: influence of beta.
clear; clc;
alpha = 0.2; d = 0.2; beta_values = [0.2,0.4,0.5];
figure; hold on;
for j = 1:numel(beta_values)
    [t,E] = kdv_time_delay(alpha,beta_values(j),d);
    plot(t,log(E),'LineWidth',1.2);
end
hold off;
legend('\beta=0.2','\beta=0.4','\beta=0.5','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
