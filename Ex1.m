close all
t_max = 200;
%%%%%%%% Part 1 %%%%%%%%%%%
dt = 0.01;
[V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, 200);
figure(1);
plotAll(time, V, i_m, sodium, pot, i_e, [0, t_max])
%%%%%% Part 2 %%%%%%%%%%%%%

figure(2);
plotAll(time, V, i_m, sodium, pot, i_e, [t_spk(1) - (dt * 100), t_spk(1) + (dt * 400)])

%%%%%% Part 3 %%%%%%%%%%%%%
%%%% I_1 = 19 %%%%%%
[V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, 19);
figure(3);
plotAll(time, V, i_m, sodium, pot, i_e, [0, t_max])

%%%%%% Part 4 %%%%%%%%%%%%%
figure(4)
%%%%%%%%%%%%
rate_thresh = 20; % 20 seconds per firing considered "repetitive firing."
%%%%%%%%%%%%
i_0 = 100;
mean_spk_int = 0.0000001; % 
while 1 / mean_spk_int > rate_thresh
    [V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, i_0);
    disp(i_0)
    i_0 = i_0 + 1;
end
plotAll(time, V, i_m, sodium, pot, i_e, [0, t_max])
%%%%%% I_rh = 109 %%%%%%%

%%%%%% Part 5 %%%%%%%%%%%%%
figure(5);
i_0 = 50;
it = 100;
[current, freq] = deal(zeros(1, it));
t_max = 3000;
for i = 1 : it
    [V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, frq] = calculateAll(t_max, dt, i_0);
    disp(i_0)
    if frq < rate_thresh
       freq(i) = 0; 
    else
        freq(i) = frq;
    end
    current(i) = i_0;
    i_0 = i_0 + 1;
end
plot(current, freq);
xlabel("I_0 (pA)")
ylabel("Firings per second", "FontSize", 12, "FontWeight", "bold")
title("Ex 1 Part 5", "FontSize", 12, "FontWeight", "bold")

%%%%%% Ex 2 Part 3 %%%%%%%%%%%%%
figure(6)
dt = 0.05;
t_max = 5000;
[V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, 200);
plot(time, V)
xlabel("time (ms)", "FontSize", 12, "FontWeight", "bold")
ylabel("V (mV)", "FontSize", 12, "FontWeight", "bold")
title("HH Model with dt = 0.05")
dt = 0.01;
tic
[V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, 200);
disp(dt)
toc
dt = 0.002;
tic
[V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, 200);
disp(dt)
toc





function [V, i_m, sodium, pot, i_e, t_spk, mean_spk_int, time, freq] = calculateAll(t_max, dt, i_0)
    t_0 = 40;
    g_na = 400;
    g_k = 200;
    g_l = 2;
    e_na = 99;
    e_k = -85;
    v_l = -65;
    C = 2;
    [V, m, h, n, i_m, i_e, sodium, pot, time] = deal(zeros(1, t_max/dt));
    % Initialize variables.
    V(1) = v_l;
    
    i_e(1) = 0;
    m(1) = inf_V(V(1), "m");
    h(1) = inf_V(V(1), "h");
    n(1) = inf_V(V(1), "n");
    sodium(1) = sodiumCond(g_na, m(1), h(1));
    pot(1) = potCond(g_k, n(1));
    
    i_m(1) = memCurrent(g_l, V(1), v_l, sodium(1), e_na, pot(1), e_k);
    % Done initializing.

    v_spk = -15;
    t_spk = zeros(1, t_max / dt);
    spk_int = zeros(1, t_max / dt);
    t_spk(:) = -1;
    i_spk = 1;
    spike_state = 0;
    for i = 1 : t_max / dt
        if V(i) > v_spk % Detect spike.
            if ~spike_state
               spike_state = 1;
               if i_spk == 1 % Calculate interval from last spike.
                   spk_int(i_spk) = time(i) - t_0;  
               else
                   spk_int(i_spk) = time(i) - t_spk(i_spk - 1);
               end
               t_spk(i) = time(i);
               i_spk = i_spk + 1;
            end
        else
           spike_state = 0; 
        end
        m(i+1) = m(i) + (dgate_dt(V(i), m(i), "m") * dt);
        h(i+1) = h(i) + (dgate_dt(V(i), h(i), "h") * dt);
        n(i+1) = n(i) + (dgate_dt(V(i), n(i), "n") * dt);
        sodium(i+1) = sodiumCond(g_na, m(i+1), h(i+1));
        pot(i+1) = potCond(g_k, n(i+1));
        if dt * i < t_0
           i_e(i+1) = 0;
        else
            i_e(i+1) = i_0;
        end
        i_m(i+1) = memCurrent(g_l, V(i), v_l, sodium(i), e_na, pot(i), e_k);
         
        V(i+1) = V(i) + memPot_dt(i_m(i), i_e(i), C) * dt;
        time(i+1) = time(i) + dt;
    end
    t_spk = t_spk(t_spk >= 0);
    mean_spk_int = mean(spk_int);
    length(t_spk)
    freq = (length(t_spk) / t_max) * 1000;
end

function plotAll(time, V, i_m, sodium, pot, ...
    i_e, limits)
    subplot(5,1,1);
    plot(time, V)
    xlabel("time (ms)", "FontSize", 12, "FontWeight", "bold")
    ylabel("V (mV)", "FontSize", 12, "FontWeight", "bold")
    xlim(limits)
    
    subplot(5,1,2);
    plot(time, i_m)
    xlabel("time (ms)", "FontSize", 12, "FontWeight", "bold")
    ylabel("I_m (pA)", "FontSize", 12, "FontWeight", "bold")
    xlim(limits)
    
    subplot(5,1,3);
    plot(time, sodium)
    xlabel("time (ms)", "FontSize", 12, "FontWeight", "bold")
    ylabel("Na Cond. (nS)", "FontSize", 12, "FontWeight", "bold")
    xlim(limits)
    
    subplot(5,1,4);
    plot(time, pot)
    xlabel("time (ms)", "FontSize", 12, "FontWeight", "bold")
    ylabel("K Cond. (nS)", "FontSize", 12, "FontWeight", "bold")
    xlim(limits)
    
    subplot(5,1,5);
    plot(time, i_e)
    xlabel("time (ms)", "FontSize", 12, "FontWeight", "bold")
    ylabel("I_e (pA)", "FontSize", 12, "FontWeight", "bold")
    xlim(limits)
end

function y = inf_V(V, type)
    if strcmp(type, "m")
        [a_mV, b_mV] = mTrans(V);
        y = a_mV / (a_mV + b_mV);
    elseif strcmp(type, "h")
        [a_hV, b_hV] = hTrans(V);
        y = a_hV / (a_hV + b_hV);
    else % n (not checking for errors).
        [a_nV, b_nV] = nTrans(V);
        y = a_nV / (a_nV + b_nV);
    end
end

function [a_mV, b_mV] = mTrans(V)
    num = 0.1 * (V + 40);
    a_exp = -0.1 * (V + 40);
    den = 1 - exp(a_exp);
    a_mV = num / den;
    b_exp = -0.0556 * (V + 65);
    b_mV = 4 * exp(b_exp);
end

function [a_hV, b_hV] = hTrans(V)
    a_exp = -0.05 * (V + 65);
    a_hV = 0.07 * exp(a_exp);
    num = 1;
    b_exp = -0.1 * (V + 35);
    den = 1 + exp(b_exp);
    b_hV = num / den;
end

function [a_nV, b_nV] = nTrans(V)
    num = 0.01 * (V + 55);
    a_exp = -0.1 * (V + 55);
    den = 1 - exp(a_exp);
    a_nV = num / den;
    b_exp = -0.0125 * (V + 65);
    b_nV = 0.125 * exp(b_exp);
end

function y = dgate_dt(V, prev, type)
    if strcmp(type, "m")
        [a, b] = mTrans(V); 
    elseif strcmp(type, "h")
        [a, b] = hTrans(V); 
    else
        [a, b] = nTrans(V); 
    end
    comp1 = a * (1 - prev);
    comp2 = b * prev;
    y = comp1 - comp2;
end


function y = sodiumCond(g_na, m, h)
    y = g_na * (m .^ 3) * h;
end

function y = potCond(g_k, n)
    y = g_k * (n .^ 4);
end

function y = memCurrent(g_l, V, v_l, sodiumCond, e_na, potCond, e_k)
    comp1 = g_l * (V - v_l);
    comp2 = sodiumCond * (V - e_na);
    comp3 = potCond * (V - e_k);
    y = (comp1 + comp2 + comp3);
end

function y = memPot_dt(memCurrent, i_e, C)
    y = (i_e - memCurrent) / C;
end