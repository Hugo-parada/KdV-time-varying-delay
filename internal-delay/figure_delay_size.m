% Internal delay: influence of the maximum delay M.
clear; clc;
M_values = [1.5,2,2.5];
figure; hold on;
for j = 1:numel(M_values)
    [t,E] = kdvint_timedelay(M_values(j));
    plot(t,log(E),'LineWidth',1.2);
end
hold off;
legend('M=1.5','M=2','M=2.5','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
