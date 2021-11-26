blue = [0 0.216 0.463];
red = [0.694 0 0.110];
BW = 500000;
SF = 11;
sample_rate = 1e6;
chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
pkt_size = 12.25 + 18;
move_rate = 3e8 / 902e6;
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');

filename = "";
fid = fopen(filename,'rb');
raw_data =fread(fid,[2 inf],'float');
fclose(fid);
raw_data_complex = complex(raw_data(1,:),raw_data(2,:));
[data, position,dechirp_result,CFO_rate] = PacketContentDetectV3(raw_data_complex,BW,SF,sample_rate,pkt_size);
[diff1, diff2,time_index_integration] = CurveFitRecoverV2(dechirp_result, chirp_samples, sample_rate, BW, CFO_rate);
time_window = 1;
diff = diff1;
index_ls = time_index_integration;
BreathFilter(diff, index_ls, time_window, sample_rate, move_rate,filename);