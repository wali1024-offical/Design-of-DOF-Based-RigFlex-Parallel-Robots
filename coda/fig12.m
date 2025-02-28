
    clear all; clc;

    %% ================ 1. 参数设定 =================
    tol = 1e-8;        % 收敛容忍度
    max_iter = 300;    % 最大迭代次数  
    E = 2e11;
    r = 0.08;
    h = 0.004;
    w = 0.01;
    W = 0.05*sqrt(3);
    L = pi*r;
    I = (h^3*w)/12;
    d = W/(2*cos(pi/6));

    xa = 0;
    ya = -d;
    xb = d*cos(pi/6);
    yb = d*sin(pi/6);
    xc = d*cos((5*pi)/6);
    yc = d*sin((5*pi)/6);

    n = 3;  % 多项式阶数 - 1

    % 初始猜测 (40 个量, 与 all_vars 顺序对应)
    x_init = [
        -0.5745; 8.7427; -17.3985; 56.8071;          % c1_0...c1_3
         1.5201; 8.7421; -17.3963; 56.7825;          % c2_0...c2_3
         3.6142; 8.7441; -17.3984; 56.7913;          % c3_0...c3_3
       -528.9103; 354.2417; 174.6685; 103.4507; 406.1288; -509.5794;  % Fx1..Fy3 (6个)
       -18.5388; -18.5837; -18.5556;                % Mo1, Mo2, Mo3
       -0.1443;                                     % fai
       -528.9103; 103.4507;                          % FT1x, FT1y
        354.2417; 406.1288;                          % FT2x, FT2y
        174.6685; -509.5794;                         % FT3x, FT3y
       -40.0766; -40.0833; -40.0617;                % M1, M2, M3
        528.9103; 103.4507; 354.2417; 406.1288; 174.6685; -509.5794   % Fm1x..Fm3y (6个)
    ];

    %% =========== 2. 创建符号变量 (含 syms afa) ==========
    syms s real
    syms afa real  % 将 afa 声明为符号量

    % 参考线部分
    sita10 = -pi/2 + ((2*pi)*(1-1))/3 + s/r;
    sita20 = -pi/2 + ((2*pi)*(2-1))/3 + s/r;
    sita30 = -pi/2 + ((2*pi)*(3-1))/3 + s/r;

    xo1 = xa - int(cos(sita10), s, 0, L);
    yo1 = ya - int(sin(sita10), s, 0, L);
    xo2 = xb - int(cos(sita20), s, 0, L);
    yo2 = yb - int(sin(sita20), s, 0, L);
    xo3 = xc - int(cos(sita30), s, 0, L);
    yo3 = yc - int(sin(sita30), s, 0, L);

    Lo = 0.08;
    xm1 = xo1 - Lo*cos(0);
    ym1 = yo1 - Lo*sin(0);
    xm2 = xo2 - Lo*cos(pi/6 + pi/2);
    ym2 = yo2 - Lo*sin(pi/6 + pi/2);
    xm3 = xo3 - Lo*cos((5*pi)/6 + pi/2);
    ym3 = yo3 - Lo*sin((5*pi)/6 + pi/2);

    % 多项式系数 (c1, c2, c3)
    c1 = sym('c1_', [n+1, 1]);
    c2 = sym('c2_', [n+1, 1]);
    c3 = sym('c3_', [n+1, 1]);

    % 其他力、弯矩、角度等未知量
    Fx1 = sym('Fx1'); Fx2 = sym('Fx2'); Fx3 = sym('Fx3');
    Fy1 = sym('Fy1'); Fy2 = sym('Fy2'); Fy3 = sym('Fy3');
    Mo1 = sym('Mo1'); Mo2 = sym('Mo2'); Mo3 = sym('Mo3');
    fai = sym('fai');
    M1 = sym('M1'); M2 = sym('M2'); M3 = sym('M3');
    FT1x = sym('FT1x'); FT1y = sym('FT1y');
    FT2x = sym('FT2x'); FT2y = sym('FT2y');
    FT3x = sym('FT3x'); FT3y = sym('FT3y');
    Fm1x = sym('Fm1x'); Fm1y = sym('Fm1y');
    Fm2x = sym('Fm2x'); Fm2y = sym('Fm2y');
    Fm3x = sym('Fm3x'); Fm3y = sym('Fm3y');

    %% 定义三根杆“根部”角度 (符号表达式, 用 afa)
    sitam1 = 0 - afa;                      
    sitam2 = pi/6 + pi/2 - afa;
    sitam3 = 5*pi/6 + pi/2 - afa;

    %% ============== 3. 搭建方程 F1,F2,F3,...F33 =============
    % 第1根
    sita1    = sum(c1 .* s.^(0:n)');
    dd_sita1 = diff(sita1, s, 2);
    d_sita1  = diff(sita1, s, 1);
    sita1_0  = subs(sita1, s, 0);
    d_sita1_L= subs(d_sita1, s, L);
    fx1 = dd_sita1 + (Fy1/(E*I))*cos(sita1) + (Fx1/(E*I))*sin(sita1);

    F1 = sym(zeros(n+3, 1));
    F1(1)    = sita1_0 - (afa - pi/2 + 2*pi*(1-1)/3);
    F1(n+3)  = d_sita1_L - (Mo1/(E*I) + 1/r);
    for j = 1 : (n+1)
        df_Fc1 = diff(sita1, c1(j));
        F1(j+1) = int(df_Fc1 * fx1, s, 0, L);
    end

    % 第2根
    sita2    = sum(c2 .* s.^(0:n)');
    dd_sita2 = diff(sita2, s, 2);
    d_sita2  = diff(sita2, s, 1);
    sita2_0  = subs(sita2, s, 0);
    d_sita2_L= subs(d_sita2, s, L);
    fx2 = dd_sita2 + (Fy2/(E*I))*cos(sita2) + (Fx2/(E*I))*sin(sita2);

    F2 = sym(zeros(n+3,1));
    F2(1)     = sita2_0 - (afa - pi/2 + 2*pi*(2-1)/3);
    F2(n+3)   = d_sita2_L - (Mo2/(E*I) + 1/r);
    for j = 1 : (n+1)
        df_Fc2 = diff(sita2, c2(j));
        F2(j+1) = int(df_Fc2 * fx2, s, 0, L);
    end

    % 第3根
    sita3    = sum(c3 .* s.^(0:n)');
    dd_sita3 = diff(sita3, s, 2);
    d_sita3  = diff(sita3, s, 1);
    sita3_0  = subs(sita3, s, 0);
    d_sita3_L= subs(d_sita3, s, L);
    fx3 = dd_sita3 + (Fy3/(E*I))*cos(sita3) + (Fx3/(E*I))*sin(sita3);

    F3 = sym(zeros(n+3,1));
    F3(1)     = sita3_0 - (afa - pi/2 + 2*pi*(3-1)/3);
    F3(n+3)   = d_sita3_L - (Mo3/(E*I) + 1/r);
    for j = 1 : (n+1)
        df_Fc3 = diff(sita3, c3(j));
        F3(j+1) = int(df_Fc3 * fx3, s, 0, L);
    end

    % 其余平衡方程
    F4 = Fx1 + Fx2 + Fx3;
    F5 = Fy1 + Fy2 + Fy3;
    F6 = Mo1 + Mo2 + Mo3 + Fx2*W*sin(fai) + Fy2*W*cos(fai) ...
         - Fx1*W*cos(pi/3 - fai) + Fy1*W*sin(pi/3 - fai);

    % 计算端点转角
    sita1_L = subs(sita1, s, L);
    sita2_L = subs(sita2, s, L);
    sita3_L = subs(sita3, s, L);

    F7 = fai - (sita1_L - pi/2);
    F8 = fai - (sita2_L - 7*pi/6);
    F9 = fai - (sita3_L - 11*pi/6);

    % 末端坐标 (含 Lo*cos(...) 等)
    x1L = xm1 + Lo*cos(sitam1) + int(cos(sita1), s, 0, L);
    y1L = ym1 + Lo*sin(sitam1) + int(sin(sita1), s, 0, L);

    x2L = xm2 + Lo*cos(sitam2) + int(cos(sita2), s, 0, L);
    y2L = ym2 + Lo*sin(sitam2) + int(sin(sita2), s, 0, L);

    x3L = xm3 + Lo*cos(sitam3) + int(cos(sita3), s, 0, L);
    y3L = ym3 + Lo*sin(sitam3) + int(sin(sita3), s, 0, L);

    % 计算柔性杆起点坐标
    xT1 = xm1 + Lo*cos(sitam1);
    yT1 = ym1 + Lo*sin(sitam1);
    xT2 = xm2 + Lo*cos(sitam2);
    yT2 = ym2 + Lo*sin(sitam2);
    xT3 = xm3 + Lo*cos(sitam3);
    yT3 = ym3 + Lo*sin(sitam3);

    F10 = x2L - x1L - W*cos(pi/3 + fai);
    F11 = y2L - y1L - W*sin(pi/3 + fai);
    F12 = -(x3L - x2L) - W*cos(fai);
    F13 = -(y3L - y2L) - W*sin(fai);
    F14 = x1L - x3L - W*cos(pi/3 - fai);
    F15 = -(y1L - y3L) - W*sin(pi/3 - fai);

    % 杆1
    F16 = FT1x - Fx1;
    F17 = FT1y - Fy1;
    F18 = M1 - Mo1 - FT1y*(x1L-xT1) - FT1x*(y1L-yT1);
    F19 = FT1x - Fm1x;
    F20 = FT1y - Fm1y;
    F21 = M1 + FT1y*(xT1-xm1) + FT1x*(yT1-ym1);

    % 杆2
    F22 = FT2x - Fx2;
    F23 = FT2y - Fy2;
    F24 = M2 - Mo2 - FT2y*(x2L-xT2) - FT2x*(y2L-yT2);
    F25 = FT2x - Fm2x;
    F26 = FT2y - Fm2y;
    F27 = M2 + FT2y*(xT2-xm2) + FT2x*(yT2-ym2);

    % 杆3
    F28 = FT3x - Fx3;
    F29 = FT3y - Fy3;
    F30 = M3 - Mo3 - FT3y*(x3L-xT3) - FT3x*(y3L-yT3);
    F31 = FT3x - Fm3x;
    F32 = FT3y - Fm3y;
    F33 = M3 + FT3y*(xT3-xm3) + FT3x*(yT3-ym3);

    F_all = [F1; F2; F3; F4; F5; F6; F7; F8; F9; ...
             F10; F11; F12; F13; F14; F15; ...
             F16; F17; F18; F19; F20; F21; ...
             F22; F23; F24; F25; F26; F27; ...
             F28; F29; F30; F31; F32; F33];

    all_vars = [
        c1; c2; c3; 
        Fx1; Fx2; Fx3; Fy1; Fy2; Fy3; 
        Mo1; Mo2; Mo3; fai; 
        FT1x; FT1y; FT2x; FT2y; FT3x; FT3y;
        M1; M2; M3; 
        Fm1x; Fm1y; Fm2x; Fm2y; Fm3x; Fm3y
    ];

    %% ============ 4. 构造数值函数 ============
    J_all = jacobian(F_all, all_vars);
    F_func = matlabFunction(F_all, 'Vars', {all_vars, afa}, 'Outputs', {'F'});
    J_func = matlabFunction(J_all, 'Vars', {all_vars, afa}, 'Outputs', {'J'});

    %% ============ 5. 外层循环: afa 从 0 到 1.0, 叠加绘图 ============
    % 采用 parula 色图为每个 α 取不同颜色，并全部用实线绘制
    % “首端刚性构件”（连接各杆末端）随当前迭代的柔性部分颜色变化，
    % “根部刚性构件”固定采用蓝紫色

    alpha_range = 0 :0.1 : 1.0;
    colors = parula(length(alpha_range));  % 使用 parula 色图
    baseColor = [0.3, 0.3, 0.7];           % 固定蓝紫色，用于根部刚性构件

    figure; 
    hold on;
    set(gcf, 'Color', 'w');  % 白色背景

    % 设置坐标轴标签和范围
    xlabel('X (m)');
    ylabel('Y (m)');
    axis equal; grid on; box on;
    xlim([-0.3, 0.3]); ylim([-0.3, 0.3]);

    % 添加 colorbar 显示 α 值的色图映射
    colormap(parula);
    caxis([min(alpha_range) max(alpha_range)]);
    cb = colorbar;
    cb.Label.String = '';

    % 用于存储每个 α 对应的绘图句柄（只取第一根杆的曲线作为代表）
    h_alpha = gobjects(length(alpha_range),1);

    for i = 1:length(alpha_range)
        alpha = alpha_range(i);
        currentColor = colors(i, :);
        fprintf('======= 当前 afa = %.2f =======\n', alpha);
        
        % 每次都从同一初值开始
        x = x_init;
        
        % 牛顿迭代求解
        for iter = 1:max_iter
            F_val = F_func(x, alpha);
            J_val = J_func(x, alpha);
            
            F_d = double(F_val);
            J_d = double(J_val);
            
            D = -J_d \ F_d;
            x = x + D;
            
            delta_norm = norm(D);
            fprintf('   iter=%d, 增量范数=%.6e\n', iter, delta_norm);
            if delta_norm < tol
                fprintf('   -> 收敛\n\n');
                break;
            end
        end
        if iter == max_iter
            warning('  * 在 %d 次迭代内未收敛.', max_iter);
        end
        
        % 提取多项式系数 C1, C2, C3
        C1 = x(1 : n+1);
        C2 = x(n+2 : 2*(n+1));
        C3 = x(2*(n+1)+1 : 3*(n+1));
        
        % 计算形状并绘图
        s_vals = linspace(0, L, 100);
        syms ss real
        theta_1_expr = sum(C1 .* ss.^(0:n)');
        theta_2_expr = sum(C2 .* ss.^(0:n)');
        theta_3_expr = sum(C3 .* ss.^(0:n)');
        
        theta_1_vals = double(subs(theta_1_expr, ss, s_vals));
        theta_2_vals = double(subs(theta_2_expr, ss, s_vals));
        theta_3_vals = double(subs(theta_3_expr, ss, s_vals));
        
        xm1_val = double(xm1);
        ym1_val = double(ym1);
        xm2_val = double(xm2);
        ym2_val = double(ym2);
        xm3_val = double(xm3);
        ym3_val = double(ym3);
        
        % 将 sitam1,2,3 代入当前的 α
        sitam1_val = double(subs(sitam1, afa, alpha));
        sitam2_val = double(subs(sitam2, afa, alpha));
        sitam3_val = double(subs(sitam3, afa, alpha));
        
        % 三条杆的离散坐标计算
        x_def = zeros(3, length(s_vals));
        y_def = zeros(3, length(s_vals));
        
        % 第1根杆
        x_def(1,1) = xm1_val + Lo*cos(sitam1_val);
        y_def(1,1) = ym1_val + Lo*sin(sitam1_val);
        for k = 2:length(s_vals)
            ds = s_vals(k) - s_vals(k-1);
            x_def(1,k) = x_def(1,k-1) + cos(theta_1_vals(k-1))*ds;
            y_def(1,k) = y_def(1,k-1) + sin(theta_1_vals(k-1))*ds;
        end
        
        % 第2根杆
        x_def(2,1) = xm2_val + Lo*cos(sitam2_val);
        y_def(2,1) = ym2_val + Lo*sin(sitam2_val);
        for k = 2:length(s_vals)
            ds = s_vals(k) - s_vals(k-1);
            x_def(2,k) = x_def(2,k-1) + cos(theta_2_vals(k-1))*ds;
            y_def(2,k) = y_def(2,k-1) + sin(theta_2_vals(k-1))*ds;
        end
        
        % 第3根杆
        x_def(3,1) = xm3_val + Lo*cos(sitam3_val);
        y_def(3,1) = ym3_val + Lo*sin(sitam3_val);
        for k = 2:length(s_vals)
            ds = s_vals(k) - s_vals(k-1);
            x_def(3,k) = x_def(3,k-1) + cos(theta_3_vals(k-1))*ds;
            y_def(3,k) = y_def(3,k-1) + sin(theta_3_vals(k-1))*ds;
        end
        
        % 绘制柔性杆（全部实线），对第1根杆保存句柄用于图例
        for rod = 1:3
            if rod == 1
                h_alpha(i) = plot(x_def(rod,:), y_def(rod,:), '-', 'Color', currentColor, 'LineWidth', 2);
            else
                plot(x_def(rod,:), y_def(rod,:), '-', 'Color', currentColor, 'LineWidth', 2);
            end
        end
        
        % 连接各杆末端刚性构件（首端），采用当前颜色
        X1 = x_def(1,end); Y1 = y_def(1,end);
        X2 = x_def(2,end); Y2 = y_def(2,end);
        X3 = x_def(3,end); Y3 = y_def(3,end);
        plot([X1, X2],[Y1, Y2], '-', 'Color', currentColor, 'LineWidth', 1.5);
        plot([X2, X3],[Y2, Y3], '-', 'Color', currentColor, 'LineWidth', 1.5);
        plot([X3, X1],[Y3, Y1], '-', 'Color', currentColor, 'LineWidth', 1.5);
        
        % 连接“根部”刚性构件——固定蓝紫色
        plot([xm1_val, x_def(1,1)], [ym1_val, y_def(1,1)], 'Color', baseColor, 'LineWidth', 1.5);
        plot([xm2_val, x_def(2,1)], [ym2_val, y_def(2,1)], 'Color', baseColor, 'LineWidth', 1.5);
        plot([xm3_val, x_def(3,1)], [ym3_val, y_def(3,1)], 'Color', baseColor, 'LineWidth', 1.5);
    end

    % 构造图例字符串，使用希腊字母 α 表示
    legend_entries = cell(length(alpha_range), 1);
    for i = 1:length(alpha_range)
        legend_entries{i} = ['\alpha = ', num2str(alpha_range(i), '%.1f')];
    end
    legend(h_alpha, legend_entries, 'Location', 'best', 'Interpreter', 'tex');

    hold off;

