function [data, position, dechirp_result,CFO_rate] = PacketContentDetectV3(raw_data_complex, BW, SF, sample_rate, pkt_size)
preknowlege_pkt_length = [35,45];
feature_upchirp_num = 2;
corr_threshold = 0.001;
abs_threshold = 0.02;
Wp = 0.3e4/5e5;
Ws = 0.6e4/5e5;
Rp = 1;
Rs = 40;
[n,wc2] = cheb1ord(Wp,Ws,Rp,Rs);
[b,a] = cheby1(n,Rp,wc2);
chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');
data = raw_data_complex;
featureChirp =[repmat(up_chirp,1,feature_upchirp_num),repmat(down_chirp,1,2), down_chirp(1:chirp_samples/4)];
data_phase = unwrap(angle(data));
[corr_result, lag] = xcorr(data, featureChirp);
corr_result = corr_result / length(featureChirp);
[~, sorted_idx] = sort(corr_result, 'descend');
CFO_rate_candidate = zeros([1 100]);
for i = 1 : 100
    CFO_detect_idx = lag(sorted_idx(i)) - 2 * chirp_samples;
    if CFO_detect_idx < 0
        continue;
    end
    CFO_rate_candidate(i) = (data_phase(CFO_detect_idx + 2 * chirp_samples) - data_phase(CFO_detect_idx)) / 2 / chirp_samples;
end
[~, ~, ~, ~, CFO_rate] = filloutliers(CFO_rate_candidate, 'center');
data = raw_data_complex .*  exp(-1j * (1:length(raw_data_complex)) * CFO_rate);
[corr_result, lag] = xcorr(data, featureChirp);
corr_result = corr_result / length(featureChirp);
position = [];
i = length(data) + 1000;
while abs(data(lag(i))) > abs_threshold 
    i = i + 1;
end
pkt_start_compensate = (10 - feature_upchirp_num) * chirp_samples;
guard = 1000;
pkt_length_ls = [];
while (i <= length(corr_result) -  pkt_size * chirp_samples + pkt_start_compensate - 500)
    if abs(data(lag(i))) < abs_threshold || min(abs(data(lag(i) : lag(i) + 1000))) < abs_threshold
        i = i + 1;
        continue;
    end
    [~,I] = max(abs(corr_result(i + pkt_start_compensate: i+ pkt_start_compensate+ 800)));
    i = i - 1 + I;
    ed = i + 12.25 * chirp_samples;
    ed_cnt = 0;
    while ed < length(corr_result)
        if abs(data(lag(ed))) > abs_threshold
            ed_cnt = 0;
        else
            ed_cnt = ed_cnt + 1;
            if ed_cnt > 1000
                break;
            end
        end
        ed = ed + 1;
    end
    if ed >= length(corr_result)
        break;
    end
    pkt_size_empirical = round((ed - i) / chirp_samples - 0.25) + 0.25;
    if pkt_size_empirical < 15
        i = i + pkt_size_empirical * chirp_samples - 1;
        continue;
    end
    pkt_length_ls = [pkt_length_ls, pkt_size_empirical];
    pkt_position = lag(i);
    y = fft(data(pkt_position + chirp_samples: pkt_position + chirp_samples * 2 - 1).* down_chirp);
    [~,idx] = max(abs(y));
    if idx > chirp_samples / 2
        idx = idx - chirp_samples;
    end
    offset = [2 * (1 - idx) - 1, 2 * (1 - idx), 2 * (1 - idx) + 1];
    position_candidate = [0,0,0];
    for j = 1:3
        temp = abs(fft(data(pkt_position + chirp_samples + offset(j): pkt_position + chirp_samples * 2 - 1 + offset(j)).* down_chirp)); 
        position_candidate(j) = temp(1);
    end
    [~, idx] = max(position_candidate);
    position = [position (pkt_position + offset(idx))];
    pkt_abs = abs(data(position(end) : position(end) + chirp_samples * pkt_size_empirical - 1));
    pkt_abs(pkt_abs == 0) = 1;
    abs_mean = mean(pkt_abs);
    zero_mean_abs = pkt_abs - abs_mean;
    filter_result = filter(b,a,zero_mean_abs);
    filterd_abs = pkt_abs;
    filterd_abs(1:end - 200) =pkt_abs(1:end - 200) - filter_result(201:end);
    data(position(end) : position(end) + chirp_samples * pkt_size_empirical - 1) = data(position(end) : position(end) + chirp_samples * pkt_size_empirical - 1) ./ pkt_abs .* filterd_abs;
    i = i + pkt_size_empirical * chirp_samples - 1;
end
dechirp_result = {};
for i = 1 : length(position) - 1
    start_chirp = 2;
    dechirp_signal = [];
    time_index = [];
    flip_index = [];
    f0_ls = [];
    break_flag = false;
    for j = start_chirp : 10
        data_unnormalized = data(position(i) + (j - 1) * chirp_samples: position(i) + j * chirp_samples - 1);
        data_normalized = data_unnormalized ./ abs(data_unnormalized);
        y = fft(data_normalized.* down_chirp);
        [M1,idx] = max(abs(y));
        y(idx) = 0;
        M2 = max(abs(y));
        
        if (M1+M2) / chirp_samples < 0.6
            break_flag = true;
            break;
        end
        if idx > chirp_samples / 2
            idx = idx - chirp_samples / 2;
        end
        f_start = BW / 2 - (idx - 1) * sample_rate / chirp_samples;
        f_k = -BW / chirp_duration;
        t_flip = (- BW / 2 - f_start) / f_k;
        conj_chirp = [chirp(1/sample_rate:1/sample_rate:t_flip, f_start, chirp_duration, f_start - BW, 'linear',0,'complex') chirp(t_flip+ 1/sample_rate:1/sample_rate:chirp_duration, f_start + BW, chirp_duration, f_start, 'linear',0,'complex')];
        dechirp_signal(end + 1, :) = data_unnormalized .* conj_chirp;
        flip_index(end + 1) = (chirp_samples - (idx - 1) * 2);
        f0_ls(end + 1) =  -f_start;
        time_index(end + 1) = position(i) + (j - 1) * chirp_samples;
    end
    if break_flag
        dechirp_result(end + 1,:) = {dechirp_signal, flip_index, time_index, f0_ls};
        continue;
    end
    for j = 1:2
        data_unnormalized = data(position(i) + (9 + j) * chirp_samples: position(i) + (10 + j) * chirp_samples - 1);
        data_normalized = data_unnormalized ./ abs(data_unnormalized);
        y = fft(data_normalized.* up_chirp);
        [M1,idx] = max(abs(y));
        y(idx) = 0;
        M2 = max(abs(y));
        
        if (M1+M2) / chirp_samples < 0.6
            break_flag = true;
            break;
        end 
        dechirp_signal(end + 1, :) = data_unnormalized .* up_chirp;
        flip_index(end + 1) = 0;
        f0_ls(end + 1) =  f_start;
        time_index(end + 1) = position(i) + (9 + j) * chirp_samples;
    end
    if break_flag
        dechirp_result(end + 1,:) = {dechirp_signal, flip_index, time_index, f0_ls};
        continue;
    end
    for j = 1 : pkt_length_ls(i) - 12.25
        data_unnormalized = data(position(i) + (j - 1 + 12.25) * chirp_samples: position(i) + (j + 12.25) * chirp_samples - 1);
        data_normalized = data_unnormalized ./ abs(data_unnormalized);
        y = fft(data_normalized.* down_chirp);
        [M1,idx] = max(abs(y));
        y(idx) = 0;
        M2 = max(abs(y));
        
        if (M1+M2) / chirp_samples < 0.6
            break;
        end
        if idx > chirp_samples / 2
            idx = idx - chirp_samples / 2;
        end
        f_start = BW / 2 - (idx - 1) * sample_rate / chirp_samples;
        f_k = -BW / chirp_duration;
        t_flip = (- BW / 2 - f_start) / f_k;
        conj_chirp = [chirp(0:1/sample_rate:t_flip - 1/sample_rate, f_start, chirp_duration, f_start - BW, 'linear',0,'complex') chirp(t_flip:1/sample_rate:chirp_duration - 1/sample_rate, f_start + BW, chirp_duration, f_start, 'linear',0,'complex')];
        dechirp_signal(end + 1, :) = data_unnormalized .* conj_chirp;
        flip_index(end + 1) = (chirp_samples - (idx - 1) * 2);
        f0_ls(end + 1) =  -f_start;
        time_index(end + 1) = position(i) + (j - 1 + 12.25) * chirp_samples;
    end
    dechirp_result(end + 1,:) = {dechirp_signal, flip_index, time_index, f0_ls};
end

    
