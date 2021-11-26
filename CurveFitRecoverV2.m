function [diff1, diff2, time_index_integration] = CurveFitRecoverV2(dechirp_result, chirp_samples, sample_rate, BW,CFO_rate)
    pkt_num = size(dechirp_result,1);
    diff1 = [];
    diff2 = [];
    last_diff = 1 + 0j;
    zero_cluster_ref = 1 + 0j;
    wavelength = 3e8 / 902e6;
    phase_change_threshold = 10 / wavelength * 2 * pi / sample_rate;
    cnt = 1;
    time_index_integration = [];
    for pkt_idx = 1 : pkt_num
        recover_signal = dechirp_result{pkt_idx,1};
        flip_index = dechirp_result{pkt_idx,2};
        time_index = dechirp_result{pkt_idx,3};
        f0_ls = dechirp_result{pkt_idx,4};
        chirp_num = size(recover_signal,1);
        for i = 1: chirp_num
            recover_signal(i,1:100) = zeros([1 100]);
            recover_signal(i,(end - 99):end) = zeros([1 100]);
            st = max(1,flip_index(i) - 100);
            ed = min(flip_index(i) + 100, size(recover_signal,2));
            recover_signal(i,st:ed) = zeros([1 (ed - st + 1)]);
            if flip_index(i) > 2048
                xdata = 1:flip_index(i);
            else
                xdata = (flip_index(i)+1):4096;
            end
            recover_signal_unit = recover_signal(i,xdata);
            zero_idx = recover_signal_unit == 0;
            xdata(zero_idx) = [];
            recover_signal_unit(zero_idx) = [];
            ydata = unwrap(angle(recover_signal_unit));
            fun = @(x,xdata)x(3) + x(2)*xdata + x(1)*(xdata.^2);
            k0 = (ydata(end) - ydata(1)) / (xdata(end) - xdata(1));
            b0 = ydata(1) - k0 * xdata(1);
            x0=[0,k0,b0];
            options = optimoptions('lsqcurvefit','Display', 'off');
            x = lsqcurvefit(fun,x0,xdata,ydata,[],[],options);
            rotate = exp(-1j* fun([x(1),x(2),0],xdata));
            chirp_recover = recover_signal_unit.* rotate;
            zeros_idx_chirp = chirp_recover == 0;
            chirp_recover(zeros_idx_chirp) = [];
            scaler = std(abs(chirp_recover)) / std(unwrap(angle(chirp_recover)));
            [cluster_idx,C_iq,~, d] = kmeans([abs(chirp_recover)', unwrap(angle(chirp_recover))' * scaler],2,'Replicates',5);
            idx1 = cluster_idx == 1;
            idx2 = cluster_idx == 2;
            iq1 = chirp_recover(idx1);
            iq2 = chirp_recover(idx2);
            C = [mean(iq1), mean(iq2)];
            C_normalized = C ./ abs(C);
            if i > 1
                delta_idx = time_index(i) - time_index(i - 1);
                diff1_phase = abs(imag(log(last_diff / ((C(1) - C(2)) / C_normalized(2) * zero_cluster_ref) )));
                diff2_phase = abs(imag(log(last_diff / ((C(2) - C(1)) / C_normalized(1) * zero_cluster_ref) )));
                if diff1_phase < diff2_phase && diff1_phase / (delta_idx) < phase_change_threshold
                    diff1 = [diff1, (C(1) - C(2)) / C_normalized(2) * zero_cluster_ref];
                    diff2 = [diff2, (C(2) - C(1)) / C_normalized(1) * zero_cluster_ref];
                elseif diff2_phase / (delta_idx) < phase_change_threshold
                    diff1 = [diff1, (C(2) - C(1)) / C_normalized(1) * zero_cluster_ref];
                    diff2 = [diff2, (C(1) - C(2)) / C_normalized(2) * zero_cluster_ref];
                else
                    diff1 = [diff1, last_diff];
                    diff2 = [diff2, -last_diff];
                end
            else 
                if C(1) - C(2) ~= 0  && C(2) ~= 0
                    rotate0 = (last_diff / ((C(1) - C(2))));
                rotate0 = rotate0 / abs(rotate0);
                zero_cluster_ref = C(2) * rotate0;
                zero_cluster_ref = zero_cluster_ref / abs(zero_cluster_ref);
                diff1 = [diff1, (C(1) - C(2)) / C_normalized(2) * zero_cluster_ref];
                diff2 = [diff2, (C(2) - C(1)) / C_normalized(1) * zero_cluster_ref];
                end
            end
            time_index_integration = [time_index_integration, time_index(i)];
            last_diff = diff1(end);
        end
    end
end