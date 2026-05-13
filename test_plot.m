clear 
close all

data = readmatrix('2412 xfoil.xlsx');

alpha_1 = data(2:end, 1);
cl_1    = data(2:end, 2);

alpha_3 = data(2:end, 4);
cl_3    = data(2:end, 5);

alpha_6 = data(2:end, 7);
cl_6    = data(2:end, 8);

alpha_9 = data(2:end, 10);
cl_9    = data(2:end, 11);

alpha = -8:0.5:12 ;
cl_list = zeros(1,size(alpha,2));
for i=1:size(alpha,2)
    Cl = main_kth_chw(alpha(i));
    cl_list(i) = Cl;
end
CL_TAT = (pi / 180) * 2 * pi * (alpha - (-2.0772404));
fig =figure();
ax = uiaxes(fig,'Position', [10,10,  1200, 1000]);
hold(ax, "on");
plot(ax, alpha, cl_list, DisplayName='vortex panel method',LineWidth= 2);
plot(ax, alpha, CL_TAT, DisplayName='TAT results',LineWidth= 2);

plot(ax,alpha_1, cl_1, 'DisplayName', 'Xfoil Re = 1e6','LineWidth', 2)
plot(ax,alpha_3, cl_3, 'DisplayName', 'Xfoil Re = 3e6','LineWidth', 2)
plot(ax,alpha_6, cl_6, 'DisplayName', 'Xfoil Re = 6e6','LineWidth', 2)
plot(ax,alpha_9, cl_9, 'DisplayName', 'Xfoil Re = 9e6','LineWidth', 2)

ax.set('Title', title('C_l vs \alpha'), ...
    'Xlabel', xlabel('\alpha'), ...
    'Ylabel', ylabel('C_l'), ...
    'FontSize', 20, ...
    'LineWidth', 2)
grid on
legend(ax, 'Location', 'southeast')