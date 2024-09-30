i_T = cell(5,1);
for i = 1 : 20
    V_T = 350;
    R = 7; %4 + ((i-1) * 1);
    Rel = 0 + ((i-1) * 10);
    %phi = 0; % + ((i-1) * 1);
    sim_out = sim('parameter_sweep_M',3);
    var_out = sim_out.yout{3}.Values;
    var_out = [var_out.Time var_out.Data];
    i_T{i} = var_out;
    plot(var_out(:,1),var_out(:,2))
    hold on
end
hold off