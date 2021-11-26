function [filter_result2] = BreathFilter(diff, index_ls, time_window, sample_rate, move_rate,filename) 
    move_data = (unwrap(angle(diff)) - angle(diff(1))) / pi /2* move_rate * 100;
    flip_idx = [];
    flip_flag = ones(size(move_data));
    for i = 2:length(flip_flag)
        flip_flag(i) = flip_flag(i-1);
        if ismember(i, flip_idx)
            flip_flag(i) = flip_flag(i) * -1;
        end
    end
    move_diff_wo_outliers = filloutliers(move_data(2:end) - move_data(1:end - 1), 'clip','movmedian',40);

    move_data_wo_outliers = zeros(size(move_data));
    move_data_wo_outliers(1) = move_data(1);
    for i = 2:length(move_data)
        move_data_wo_outliers(i) = move_data_wo_outliers(i - 1) + move_diff_wo_outliers(i - 1) * flip_flag(i);
    end
    interp_step = 1;
    diff_interp = interp1(index_ls,move_data_wo_outliers, 1:interp_step:index_ls(end));
    Wp = 1.5/5e5;
    Ws = 3/5e5;
    Rp = 1;
    Rs = 10;
    [n,wc] = cheb1ord(Wp,Ws,Rp,Rs);
    [b,a] = cheby1(n,Rp,wc);
    filter_result = filter(b,a,diff_interp(index_ls(1):end));
    Wp2 = 0.5/5e5;
    Ws2 = 0.01/5e5;
    Rp2 = 1;
    Rs2 = 10;
    [n2,wc2] = cheb1ord(Wp2,Ws2,Rp2,Rs2);
    [b2,a2] = cheby1(n2,Rp,wc2,'high');
    filter_result2 = filter(b2,a2,filter_result);
end