i_T = cell(5,1);
for i = 1 : 5
    V_T = 4 + ((i-1) * 1);
    R = 7; %4 + ((i-1) * 1);
    sim_out = sim('parameter_sweep_LR',10);
    var_out = sim_out.yout{3}.Values;
    var_out = [var_out.Time var_out.Data];
    i_T{i} = var_out;
    plot(var_out(:,1),var_out(:,2))
    hold on
end
hold off