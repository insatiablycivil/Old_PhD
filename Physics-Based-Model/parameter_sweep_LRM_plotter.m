loop_Count = 1;
for i = 1 : loop_Count
    V_T = 33.5;
    R = 4.95; % - ((i-1) * 0.7);
    R_Int = 0.55;%0.33; % + ((i-1) * 0.7);
    A_w = 1; %3.141;
    l_w = 1; % 1990;
    Nw = 4; % + ((i-1) * 50);
    Rel = Nw^2/ (0.004); %(1/0.0486) * 0.030;% + ((i-1) * 10);
    %phi = 0; % + ((i-1) * 1);
    sim_out = sim('parameter_sweep_M',0.18);
    var_out = sim_out.yout{3}.Values;
    var_out = [var_out.Time var_out.Data];
    plot(var_out(:,1),var_out(:,2))
    hold on
    
    sim_out = sim('parameter_sweep_LR',0.18);
    var_out = sim_out.yout{3}.Values;
    var_out = [var_out.Time var_out.Data];
    plot(var_out(:,1),var_out(:,2))
end
hold off