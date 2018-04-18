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
[memPot, extCurr, time, t_spk, mean_spk_int, freq] = forward_euler(g_l, v_0, v_l, i_e, C, dt, t_max, thresh, v_r, ref_per);

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
figure(2)
it = 4;
di = 0.1;
[i_e, freq] = deal(zeros(1, it / di));
i_e(1) = 0;
freq(1) = 0;
t_max = 4000;
for i = 1 : it / di 
    [memPot, extCurr, time, t_spk, mean_spk_int, frq] = forward_euler(g_l, v_0, v_l, i_e(i), C, dt, t_max, thresh, v_r, ref_per);
    i_e(i+1) = i_e(i) + di;
    freq(i+1) = frq;
end

plot(i_e, freq)
xlabel("I_e (nA)")
ylabel("Frequency (Hz)")




function y = dv_dt(g_l, v, v_l, i_e, C)
    y = (-g_l * (v - v_l) + i_e) / C;
end

function [memPot, extCurr, time, t_spk, mean_spk_int, freq] = forward_euler(g_l, v_0, v_l, i_e, C, dt, t_max, thresh, v_r, ref_per)
    [memPot, extCurr, time, ref_counter, t_spk, spk_int] = deal(zeros(1, t_max / dt));
    memPot(1) = v_0;
    extCurr(1) = i_e;
    t_spk(:) = -1;
    i_spk = 1;
    for i = 1 : t_max / dt
        dvdt = dv_dt(g_l, memPot(i), v_l, i_e, C);
        if ref_counter(i) <= 0
           if memPot(i) >= thresh % spiked
               if i_spk == 1 % Calculate interval from last spike.
                   spk_int(i_spk) = time(i);  
               else
                   spk_int(i_spk) = time(i) - t_spk(i_spk - 1);
               end
               t_spk(i) = time(i);
               i_spk = i_spk + 1;
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
    t_spk = t_spk(t_spk >= 0);
    mean_spk_int = mean(spk_int);
    freq = (length(t_spk) / t_max) * 1000;
end