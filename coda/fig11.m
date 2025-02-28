
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

    % 原脚本给出的初始猜测 (末尾附加 ones(9,1))
    x_init = [
       -1.3714;  8.7137; 19.7603; -27.4374;
        0.7230;  8.7046; 19.7838; -27.3872;
        2.8174;  8.7166; 19.7001; -27.2442;
      161.6942; -349.0686; 187.3744; -308.6276;
       15.2253;  293.4022; 10.1011;  10.2320;
       10.2009;   0.0606; ones(9,1)
    ];

    %% =========== 2. 创建符号变量 (含 syms afa) ==========
    syms s real
    syms afa real  % 将 afa 声明为符号量

    % =========== 原脚本中的几何/参考线部分 ===========
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
    sitaT1 = sym('sitaT1'); sitaT2 = sym('sitaT2'); sitaT3 = sym('sitaT3');
    FT1x = sym('FT1x'); FT1y = sym('FT1y');
    FT2x = sym('FT2x'); FT2y = sym('FT2y');
    FT3x = sym('FT3x'); FT3y = sym('FT3y');

    %% ============== 3. 搭建方程 F1,F2,F3,...F15 =============
    % ---------- F1 ----------
    sita1    = sum(c1 .* s.^(0:n)');
    dd_sita1 = diff(sita1, s, 2);
    d_sita1  = diff(sita1, s, 1);
    sita1_0  = subs(sita1, s, 0);
    d_sita1_L= subs(d_sita1, s, L);
    fx1 = dd_sita1 + (Fy1/(E*I))*cos(sita1) + (Fx1/(E*I))*sin(sita1);

    F1 = sym(zeros(n+3, 1));
    F1(1) = sita1_0 - sitaT1;
    F1(n+3) = d_sita1_L - (Mo1/(E*I) + 1/r);
    for j = 1:(n+1)
        df_Fc1 = diff(sita1, c1(j));
        F1(j+1) = int(df_Fc1 * fx1, s, 0, L);
    end

    % ---------- F2 ----------
    sita2    = sum(c2 .* s.^(0:n)');
    dd_sita2 = diff(sita2, s, 2);
    d_sita2  = diff(sita2, s, 1);
    sita2_0  = subs(sita2, s, 0);
    d_sita2_L= subs(d_sita2, s, L);
    fx2 = dd_sita2 + (Fy2/(E*I))*cos(sita2) + (Fx2/(E*I))*sin(sita2);

    F2 = sym(zeros(n+3,1));
    F2(1) = sita2_0 - sitaT2;
    F2(n+3) = d_sita2_L - (Mo2/(E*I) + 1/r);
    for j = 1:(n+1)
        df_Fc2 = diff(sita2, c2(j));
        F2(j+1) = int(df_Fc2 * fx2, s, 0, L);
    end

    % ---------- F3 ----------
    sita3    = sum(c3 .* s.^(0:n)');
    dd_sita3 = diff(sita3, s, 2);
    d_sita3  = diff(sita3, s, 1);
    sita3_0  = subs(sita3, s, 0);
    d_sita3_L= subs(d_sita3, s, L);
    fx3 = dd_sita3 + (Fy3/(E*I))*cos(sita3) + (Fx3/(E*I))*sin(sita3);

    F3 = sym(zeros(n+3,1));
    F3(1) = sita3_0 - sitaT3;
    F3(n+3) = d_sita3_L - (Mo3/(E*I) + 1/r);
    for j = 1:(n+1)
        df_Fc3 = diff(sita3, c3(j));
        F3(j+1) = int(df_Fc3 * fx3, s, 0, L);
    end

    % ---------- 其余平衡/几何方程 ----------
    F4 = Fx1 + Fx2 + Fx3;
    F5 = Fy1 + Fy2 + Fy3;
    F6 = Mo1 + Mo2 + Mo3 + Fx2*W*sin(fai) + Fy2*W*cos(fai) ...
         - Fx1*W*cos(pi/3 - fai) + Fy1*W*sin(pi/3 - fai);

    sita1_L = subs(sita1, s, L);
    sita2_L = subs(sita2, s, L);
    sita3_L = subs(sita3, s, L);

    F7 = fai - (sita1_L - pi/2);
    F8 = fai - (sita2_L - 7*pi/6);
    F9 = fai - (sita3_L - 11*pi/6);

    % 计算末端坐标(含 Lo*cos(afa)等)
    x1L = xm1 + Lo*cos(afa) + int(cos(sita1), s, 0, L);
    y1L = ym1 + Lo*sin(afa) + int(sin(sita1), s, 0, L);

    x2L = xm2 + Lo*cos(afa + pi/6 + pi/2) + int(cos(sita2), s, 0, L);
    y2L = ym2 + Lo*sin(afa + pi/6 + pi/2) + int(sin(sita2), s, 0, L);

    x3L = xm3 + Lo*cos(afa + (5*pi)/6 + pi/2) + int(cos(sita3), s, 0, L);
    y3L = ym3 + Lo*sin(afa + (5*pi)/6 + pi/2) + int(sin(sita3), s, 0, L);

    % 计算柔性杆起点坐标
    xT1 = xm1 + Lo*cos(afa);
    yT1 = ym1 + Lo*sin(afa);
    xT2 = xm2 + Lo*cos(afa+pi/6+pi/2);
    yT2 = ym2 + Lo*sin(afa+pi/6+pi/2);
    xT3 = xm3 + Lo*cos(afa+(5*pi)/6+pi/2);
    yT3 = ym3 + Lo*sin(afa+(5*pi)/6+pi/2);

    F10 = x2L - x1L - W*cos(pi/3 + fai);
    F11 = y2L - y1L - W*sin(pi/3 + fai);
    F12 = -(x3L - x2L) - W*cos(fai);
    F13 = -(y3L - y2L) - W*sin(fai);
    F14 = x1L - x3L - W*cos(pi/3 - fai);
    F15 = -(y1L - y3L) - W*sin(pi/3 - fai);
    F16 = FT1x - Fx1;
    F17 = FT1y - Fy1;
    F18 = Mo1 + FT1y*(x1L-xT1) + FT1x*(y1L-yT1);
    F19 = FT2x - Fx2;
    F20 = FT2y - Fy2;
    F21 = Mo2 + FT2y*(x2L-xT2) + FT2x*(y2L-yT2);
    F22 = FT3x - Fx3;
    F23 = FT3y - Fy3;
    F24 = Mo3 + FT3y*(x3L-xT3) + FT3x*(y3L-yT3);

    F_all = [F1; F2; F3; F4; F5; F6; F7; F8; F9; F10; F11; F12; F13; F14; F15; ...
             F16; F17; F18; F19; F20; F21; F22; F23; F24];

    all_vars = [c1; c2; c3; Fx1; Fx2; Fx3; Fy1; Fy2; Fy3; Mo1; Mo2; Mo3; ...
                fai; sitaT1; sitaT2; sitaT3; FT1x; FT1y; FT2x; FT2y; FT3x; FT3y];

    %% ============ 4. 构造数值函数: 把 afa 也列入 'Vars' ============
    J_all = jacobian(F_all, all_vars);
    F_func = matlabFunction(F_all, 'Vars', {all_vars, afa}, 'Outputs', {'F'});
    J_func = matlabFunction(J_all, 'Vars', {all_vars, afa}, 'Outputs', {'J'});

    %% ============ 5. 外层循环: afa 在多个值上, 叠加绘图 ============
    % 将 alpha_range 取更多点, 例如 0:0.2:1.0
    alpha_range = 0:0.1:1.0;

    % 创建 figure 并配置坐标轴
    figure('Name','goodlookmore1 - Enhanced');
    hold on; grid on; box on;
    xlabel('X (m)');
    ylabel('Y (m)');
    axis square;
    set(gca, 'Position', [0.15 0.15 0.7 0.7]);
    xlim([-0.3 0.3]); 
    ylim([-0.3 0.3]);

    % 使用 parula 色表, 并配合 colorbar
    nAlpha = length(alpha_range);
    cmap = parula(nAlpha);
    colormap(parula);
    caxis([min(alpha_range) max(alpha_range)]);
    cb = colorbar;
    cb.Label.String = '';

    % 用于存储图例句柄（每个 \alpha 只保存一条曲线句柄）
    hCurves = gobjects(nAlpha, 1);

    % 循环遍历各个 \alpha 值
    for idx = 1:nAlpha
        alpha = alpha_range(idx);
        fprintf('======= 当前 α = %.2f =======\n', alpha);

        % 每次采用相同初值
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

        % 提取多项式系数
        C1 = x(1:n+1);
        C2 = x(n+2:2*(n+1));
        C3 = x(2*(n+1)+1:3*(n+1));

        % 计算三根杆的数值形状
        s_vals = linspace(0, L, 100);
        syms ss real
        theta_1_expr = sum(C1 .* ss.^(0:n)');
        theta_2_expr = sum(C2 .* ss.^(0:n)');
        theta_3_expr = sum(C3 .* ss.^(0:n)');

        theta_1_vals = double(subs(theta_1_expr, ss, s_vals));
        theta_2_vals = double(subs(theta_2_expr, ss, s_vals));
        theta_3_vals = double(subs(theta_3_expr, ss, s_vals));

        % 获取基点的数值
        xm1_val = double(xm1);
        ym1_val = double(ym1);
        xm2_val = double(xm2);
        ym2_val = double(ym2);
        xm3_val = double(xm3);
        ym3_val = double(ym3);

        % 计算三条杆的离散坐标（积分方式）
        x_def = zeros(3, length(s_vals));
        y_def = zeros(3, length(s_vals));

        % 第1根
        x_def(1,1) = xm1_val + Lo*cos(alpha);
        y_def(1,1) = ym1_val + Lo*sin(alpha);
        for k = 2:length(s_vals)
            ds = s_vals(k) - s_vals(k-1);
            x_def(1,k) = x_def(1,k-1) + cos(theta_1_vals(k-1)) * ds;
            y_def(1,k) = y_def(1,k-1) + sin(theta_1_vals(k-1)) * ds;
        end

        % 第2根
        x_def(2,1) = xm2_val + Lo*cos(alpha + pi/6 + pi/2);
        y_def(2,1) = ym2_val + Lo*sin(alpha + pi/6 + pi/2);
        for k = 2:length(s_vals)
            ds = s_vals(k) - s_vals(k-1);
            x_def(2,k) = x_def(2,k-1) + cos(theta_2_vals(k-1)) * ds;
            y_def(2,k) = y_def(2,k-1) + sin(theta_2_vals(k-1)) * ds;
        end

        % 第3根
        x_def(3,1) = xm3_val + Lo*cos(alpha + (5*pi)/6 + pi/2);
        y_def(3,1) = ym3_val + Lo*sin(alpha + (5*pi)/6 + pi/2);
        for k = 2:length(s_vals)
            ds = s_vals(k) - s_vals(k-1);
            x_def(3,k) = x_def(3,k-1) + cos(theta_3_vals(k-1)) * ds;
            y_def(3,k) = y_def(3,k-1) + sin(theta_3_vals(k-1)) * ds;
        end

        % 当前迭代的颜色
        lineColor = cmap(idx, :);

        % 以第1根的绘图句柄存进 hCurves，用于图例
        hCurves(idx) = plot(x_def(1,:), y_def(1,:), '-', 'Color', lineColor, 'LineWidth', 2);

        % 另外两根
        plot(x_def(2,:), y_def(2,:), '-', 'Color', lineColor, 'LineWidth', 2);
        plot(x_def(3,:), y_def(3,:), '-', 'Color', lineColor, 'LineWidth', 2);

        % 末端连接
        X1 = x_def(1,end); Y1 = y_def(1,end);
        X2 = x_def(2,end); Y2 = y_def(2,end);
        X3 = x_def(3,end); Y3 = y_def(3,end);

        plot([X1, X2], [Y1, Y2], '-', 'Color', lineColor, 'LineWidth', 2);
        plot([X2, X3], [Y2, Y3], '-', 'Color', lineColor, 'LineWidth', 2);
        plot([X3, X1], [Y3, Y1], '-', 'Color', lineColor, 'LineWidth', 2);

        % 连接“根部”（保持蓝紫色）
        plot([xm1_val, x_def(1,1)], [ym1_val, y_def(1,1)], ...
             '-', 'Color', [0.3 0.3 0.7], 'LineWidth', 1.5);
        plot([xm2_val, x_def(2,1)], [ym2_val, y_def(2,1)], ...
             '-', 'Color', [0.3 0.3 0.7], 'LineWidth', 1.5);
        plot([xm3_val, x_def(3,1)], [ym3_val, y_def(3,1)], ...
             '-', 'Color', [0.3 0.3 0.7], 'LineWidth', 1.5);
    end

    % 添加图例：显示不同 \alpha 值的颜色
    legend(hCurves, arrayfun(@(a)sprintf('\\alpha = %.1f', a), alpha_range, ...
           'UniformOutput', false), 'Location', 'bestoutside');

    hold off;

