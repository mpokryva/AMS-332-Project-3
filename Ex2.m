close all

g_l = 50 * 10 .^ -3;
v_l = -65;
v_0 = v_l;
C = 1;
thresh = -45;
v_r = -65;
ref_per = 2;
i_e = 1.1;
dt = 0.1;
t_max = 200;

%%%%% Part 1 %%%%%%%%%
[memPot, extCurr, time] = forward_euler(g_l, v_0, v_l, i_e, C, dt, t_max, thresh, v_r, ref_per);

figure(1)

subplot(1, 2, 1)
plot(time, memPot)
xlabel("time (ms)")
ylabel("V (mV)")

subplot(1, 2, 2)
plot(time, extCurr)
xlabel("time (ms)")
ylabel("I_e (nA)")

%%%%% Part 2 %%%%%%%%%
it = 1;
%i_e = 
for i = 1 : it
    [memPot, extCurr, time] = forward_euler(g_l, v_0, v_l, i_e, C, dt, t_max, thresh, v_r, ref_per);
end


function y = dv_dt(g_l, v, v_l, i_e, C)
    y = (-g_l * (v - v_l) + i_e) / C;
end

function [memPot, extCurr, time] = forward_euler(g_l, v_0, v_l, i_e, C, dt, t_max, thresh, v_r, ref_per)
    [memPot, extCurr, time, ref_counter] = deal(zeros(1, t_max / dt));
    memPot(1) = v_0;
    extCurr(1) = i_e;
    for i = 1 : t_max / dt
        dvdt = dv_dt(g_l, memPot(i), v_l, i_e, C);
        if ref_counter(i) <= 0
           if memPot(i) >= thresh
               memPot(i+1) = v_r;
               ref_counter(i+1) = ref_per;
               
           else
               memPot(i+1) = memPot(i) + (dvdt * dt);
           end
        else
            memPot(i+1) = v_r;
            ref_counter(i+1) = ref_counter(i) - dt;
        end
        extCurr(i+1) = C * dvdt + (g_l * (memPot(i) - v_l));
        time(i+1) = time(i) + dt;
    end
end