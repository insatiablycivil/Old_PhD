loop_Count = 1;
i_T = cell(loop_Count,1);
v_T = cell(loop_Count,1);
for i = 1 : loop_Count
    V_T = 33.5;
    R = 5; % - ((i-1) * 0.7);
    R_Int = 0.33; % + ((i-1) * 0.7);
    %R_Rate = 0.00001 + ((i-1) * 2);
    A_w = 1;%3.141; %3.141;
    l_w = 1;%1990; % 1990;
    Nw = 300;%300; %300; % + ((i-1) * 50);
    Rel = 0;%(Nw^2)/ (1);%(1*((sqrt(A_w/l_w)*Nw)^2))/ (0.0486); %; %5.5; %  + ((i-1) * 10);
    Rel2 = (Nw^2)/ (0.08); %+ ((i-1) * 0.5);
    
    Mass = 0.04;%88 - ((i-1) * 0.015);
    Mur = 1;%; + ((i-1) * 500);
    K = 0.25 ; % + ((i-1) * 5);
    K_Z = 0; %1 + ((i-1) * 5);
    G = -(Mass * 9.81);
    F_BW = 0.00001; %+ ((i-1) * 50);
    F_C = 0.00001; 
    D_V  = 0.00001;
    Ltch_St= 12;  % - ((i-1) * 0.5);
    Ltch_End = 15; % + ((i-20) * 0.5);
    Ltch_Str = 5; % * (2/(Ltch_End -  Ltch_St));  % 500 - ((i-1) * 0.5)
    %phi = 0; % + ((i-1) * 1);
    sim_out = sim('parameter_sweep_Mechanical',0.18);
    var_out = sim_out.yout{3}.Values;
    var_out = [var_out.Time var_out.Data];
    i_T{i} = var_out;
    plot(var_out(:,1),var_out(:,2))
    hold on
    
%     var_out = sim_out.yout{1}.Values;
%     var_out = [var_out.Time var_out.Data];
%     v_T{i} = var_out;
%     plot(var_out(:,1),var_out(:,2))
end
hold off