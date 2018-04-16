close all

t_max = 200;
t_0 = 40;
dt = 0.01;
g_na = 400;
g_k = 200;
g_l = 2;
e_na = 99;
e_k = -85;
v_l = -65;
C = 2;
i_0 = 200;
[V, m, h, n, i_m, i_e, sodium, pot, time] = deal(zeros(1, T_MAX/dt));
% Initialize variables.
V(1) = v_l;
i_m(1) = memCurrent(g_l, V(1), v_l, sodium(1), e_na, pot(1), e_k);
i_e(1) = 0;
m(1) = inf_V(V(1), "m");
h(1) = inf_V(V(1), "h");
n(1) = inf_V(V(1), "n");
sodium(1) = sodiumCond(g_na, m(1), h(1));
pot(1) = potCond(g_k, n(1));
% Done initializing.

v_spk = -15;
t_spk = -1;
%%%%%%%%%% Part 1 %%%%%%%%%%%%%%

figure(1);
for i = 1 : t_max / dt
    if t_spk < 0 && V(i) > v_spk
       t_spk = time(i); 
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

plotAll(time, V, i_m, sodium, pot, i_e, [0, t_max])

%%%%%% Part 2 %%%%%%%%%%%%

figure(2);
plotAll(time, V, i_m, sodium, pot, i_e, [t_spk - dt * 500, t_spk + dt * 500])


function plotAll(time, V, i_m, sodium, pot, ...
    i_e, limits)
    ax1 = subplot(5,1,1);
    plot(time, V)
    xlabel("time (ms)")
    ylabel("V (mV)")
    xlim(limits)
    ax2 = subplot(5,1,2);
    plot(time, i_m)
    xlabel("time (ms)")
    ylabel("I_m (pA)")
    xlim(limits)
    ax3 = subplot(5,1,3);
    plot(time, sodium)
    xlabel("time (ms)")
    ylabel("Na Conductance (mS)")
    xlim(limits)
    ax4 = subplot(5,1,4);
    plot(time, pot)
    xlabel("time (ms)")
    ylabel("K Conductance (mS)")
    xlim(limits)
    ax5 = subplot(5,1,5);
    plot(time, i_e)
    xlabel("time (ms)")
    ylabel("I_e (pA)")
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
    den = 1 - (exp(1) .^ a_exp);
    a_mV = num / den;
    b_exp = -0.0556 * (V + 65);
    b_mV = 4 * (exp(1) .^ b_exp);
end

function [a_hV, b_hV] = hTrans(V)
    a_exp = -0.05 * (V + 65);
    a_hV = 0.07 * (exp(1) .^ a_exp);
    num = 1;
    b_exp = -0.1 * (V + 35);
    den = 1 + (exp(1) .^ b_exp);
    b_hV = num / den;
end

function [a_nV, b_nV] = nTrans(V)
    num = 0.01 * (V + 55);
    a_exp = -0.1 * (V + 55);
    den = 1 - (exp(1) .^ a_exp);
    a_nV = num / den;
    b_exp = -0.0125 * (V + 65);
    b_nV = 0.125 * (exp(1) .^ b_exp);
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