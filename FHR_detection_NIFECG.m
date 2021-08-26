close all
%% Deteksi FHR pada sinyal NIFECG
% dengan metode ekstraksi Template Subtracting
% dengan deteksi QRS Pan-Tompkins
% Algoritma Pan-Tompkins menggunakan program MATLAB yang dibuat oleh
% Sedghamiz, H (2014). Matlab Implementation of Pan Tompkins ECG QRS detector. 
% Dataset: Matonia dkk.(2006). 

%% Parameter Kebutuhan
% Atur ID dan kanal dari Dataset
nama_file = 'Uji4';
id = '01';                          % Subjek ke-
kanal = '1';                        % Kanal ke-
batas_SNR = 0;                      % Batas SNR Seleksi Kanal
fs = 500;                           % frekuensi sampling dataset

% Data visual yang ingin dikeluarkan
cetak = 1;                          % Mencetak hasil deteksi FQRS
cetak_qrs = 0;                      % Mencetak hasil deteksi MQRS
cetak_rekons = 0;                   % Mencetak hasil rekonstruksi TS
cetak_SNR = 0;                      % Mencetak PF dan PN
cetak_proses = 	0;                  % Mencetak proses TS
cetak_raw = 0;                      % Mencetak sinyal AECG
simpan_data = 0;                
cetak_preSNR = 0;                   
cetak_residu = 0;                   % Mencetak sinyal residu
list = 0;

% Path dokumen untuk menyimpan data
PathPLOT = 'E:\TUGAS AKHIR\Processing\Uji Coba Dataset';

%% DATA AECG
%%Import data (dari format .txt)
filename = ['B1_abSignals_', id, '.txt'];
delimiter = '\t';
startRow = 1;
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
comma2point_overwrite(filename);
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

%Tabel berisi 8 kolom :
%4 kolom pertama merupakan data AECG 4 kanal setelah preprocessing (PP)
%4 kolom sinyal residu referensi
inp = table(dataArray{1:end-1}, 'VariableNames', {'PP_1','PP_2','PP_3','PP_4','R_1','R_2','R_3','R_4'});

PP_1 = inp.PP_1;
PP_2 = inp.PP_2;
PP_3 = inp.PP_3;
PP_4 = inp.PP_4;
R_1 = inp.R_1;
R_2 = inp.R_2;
R_3 = inp.R_3;
R_4 = inp.R_4;

% Segmen yang diproses untuk setiap subjek
dstart  = 1;                        % data sampel ke-
dend    = 5000;                     % sampai data sampel ke- (durasi 10 detik, 5000 sampel)
start_loop = 1;                     % Segmen ke-
end_loop = 120;                     % sampai Segmen ke-

%% DATA R MECG
%%Import data (dari format .txt)
filename = ['B1_Maternal_R_', id, '.txt'];
delimiter = '\t';
startRow = 1;
formatSpec = '%f%f%[^\n\r]';
comma2point_overwrite(filename);
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

inp = table(dataArray{1:end-1}, 'VariableNames', {'locs_RM','Mval'});

RM_anot = inp.locs_RM;
RM_val = inp.Mval;

NR_MECG = length(RM_anot);

%% DATA R FECG

%%Import data (dari format .txt)
filename = ['B1_Fetal_R_', id, '.txt'];
delimiter = '\t';
startRow = 1;
formatSpec = '%f%f%[^\n\r]';
comma2point_overwrite(filename);
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

inp = table(dataArray{1:end-1}, 'VariableNames', {'locs_RF','val'});

R_anot = inp.locs_RF;
Rval_anot = inp.val;

NR_FECG = length(R_anot);

%% Inisiasi variabel template subtraction
MHR=[];
MHR_anot=[];
FHR=[];
FHR_anot=[];
N_Rpeak=[];
tMB=[];
tMP=[];
tFD=[];
SNR_1=[];
SNR_2=[];
SNR_3=[];
SNR_4=[];
RR_period_3=[];
RR_period_2=[];
MAE=[];
MSE=[];
cek_FHR=0;
SNRf=[];

%% Iterasi (segmen 10 detik)
for loop=(start_loop:end_loop)
tSNR = [];
tkanal = [];
residu_1 = [];
residu_1 = [];
residu_1 = [];
residu_1 = [];
sum_SNR = 0;

%% RF Anotation
sum_RF = 1;
locs_R_fecg=[];
locs_R_fecg2=[];
R_fecg_val=[];

for i = 1:NR_FECG
    if R_anot(i) >= dstart && R_anot(i) <= dend
        if loop==1
           locs_R_fecg(sum_RF) = R_anot(i)-dstart-50;
        else
           locs_R_fecg(sum_RF) = R_anot(i)-dstart;
        end
        sum_RF = sum_RF+1;
    end
end

%% Loop Setiap Kanal

for kanal=1:4
    if kanal==1
        n_kanal='Kanal 1';
    elseif kanal==2
        n_kanal='Kanal 2';
    elseif kanal==3
        n_kanal='Kanal 3';
    else
        n_kanal='Kanal 4';
    end
    
if kanal==1
    ecg_raw = PP_1;
elseif kanal==2
    ecg_raw = PP_2;
elseif kanal==3
    ecg_raw = PP_3;
else
    ecg_raw = PP_4;
end

ecg = [];
if loop==1
    ecg=ecg_raw((dstart+50):dend);                                                      % tidak dimulai dari satu karena 
elseif loop==120
    ecg=ecg_raw(dstart:length(ecg_raw));
else
    ecg=ecg_raw(dstart:dend);
end

%% Cetak EKG Raw
if cetak_raw
t=(1:1:length(ecg))/fs;
figure
plot(t,ecg);
title('Sinyal EKG Abdomen (EKGa)');
xlabel('Time (s)');
ylabel('Amplitudo');
end

%% Cek Orientasi
neg=0;
pos=0;
process=1;
cek_start=1;
cek_end=1000;
N_cek=5;

PS=[];
if loop<120
for i=1:N_cek
    if N_cek==5&&loop==1
        cek_end=4949;
    end
    if abs(min(ecg(cek_start:cek_end)))>max(ecg(cek_start:cek_end))
        neg = neg+1;
    else pos = pos+1;
    end
    %cari daya setiap subsegment
    sum_PS = 0;
    for j=cek_start:cek_end
        sum_PS = sum_PS+(ecg(j)*ecg(j));
    end
    PS(i) = sum_PS/(cek_end-cek_start);
    
    cek_start=cek_start+1000;
    cek_end=cek_end+1000;
end

if neg>pos
    ecg=-ecg;
end

for i=1:N_cek
if PS(i)>1.5*median(PS)
        process=0;
end
end

else

N_cek = round((length(ecg))/1000);
for i=1:N_cek-1
    if abs(min(ecg(cek_start:cek_end)))>max(ecg(cek_start:cek_end))
        neg = neg+1;
    else pos = pos+1;
    end
    %cari daya setiap subsegment
    sum_PS = 0;
    for j=cek_start:cek_end
        sum_PS = sum_PS+(ecg(j)*ecg(j));
    end
    PS(i) = sum_PS/(cek_end-cek_start);
    
    cek_start=cek_start+1000;
    cek_end=cek_end+1000;
end

if neg>pos
    ecg=-ecg;
end

for i=1:N_cek-1
if PS(i)>1.5*median(PS)
    %disp(['Deteksi Artifak :', n_kanal]);
    %    disp(string(loop));
        process=0;
end
end
end

%% R peak detection
[R_val,locs_R,delay]=pan_tompkin(ecg,fs,0);             %algoritma Pan-Tompkins dengan perubahan

if length(locs_R)>8&&process
%% Q and S detection
%Inisiasi titik Q dan S
locs_Q = [];
locs_S = [];

ecg_inv = [];
ecg_inv = -ecg;                                                        %ecg_h = data ecg hasil BPF
N_pks = length(locs_R);
for i = 1:N_pks
     if i == 1 && (locs_R(1)-round(0.04*fs))<=0
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(1:locs_R(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R(i):locs_R(i)+round(0.04*fs)));  %0.06 = max durasi QR
            locs_Q(i) = locs_Q_temp;
            locs_S(i) = locs_R(i)+locs_S_temp-1;
     elseif i == 1 && (locs_R(1)-round(0.04*fs))>0 
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R(i)-round(0.04*fs):locs_R(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R(i):locs_R(i)+round(0.04*fs)));
            locs_Q(i) = locs_R(i)-round(0.04*fs)+locs_Q_temp-1;
            locs_S(i) = locs_R(i)+locs_S_temp-1;
     elseif i == N_pks && (locs_R(N_pks)+round(0.04*fs))>length(ecg)
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R(i)-round(0.04*fs):locs_R(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R(i):end));
            locs_Q(i) = locs_R(i)-round(0.04*fs)+locs_Q_temp-1;
            locs_S(i) = locs_R(i)+locs_S_temp-1;
     elseif i == N_pks && (locs_R(N_pks)+round(0.04*fs))<=length(ecg)
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R(i)-round(0.04*fs):locs_R(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R(i):(locs_R(N_pks)+round(0.04*fs))));
            locs_Q(i) = locs_R(i)-round(0.04*fs)+locs_Q_temp-1;
            locs_S(i) = locs_R(i)+locs_S_temp-1;
     else                
          [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R(i)-round(0.04*fs):locs_R(i)));   
          [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R(i):locs_R(i)+round(0.04*fs)));
          locs_Q(i) = locs_R(i)-round(0.04*fs)+locs_Q_temp-1;
          locs_S(i) = locs_R(i)+locs_S_temp-1; 
     end       
end

    locs_Q = locs_Q';
    locs_S = locs_S';
    locs_R = locs_R';

%disp([n_kanal, ': QRS MECG done'])
%% R MECG Anotation
k = 1;
locs_R_mecg=[];

for i = 1:NR_MECG
    if RM_anot(i) >= dstart && RM_anot(i) <= dend
        if loop==1
           locs_R_mecg(k) = RM_anot(i)-dstart-50+2;
        else
           locs_R_mecg(k) = RM_anot(i)-dstart+2;
        end
        k = k+1;
    end
end

%% R-R MECG Interval & MHR
%Inisiasi MHR dan MHR anotasi
RRM_period = [];
RRM_period_anot = [];

RRM_period = diff(locs_R);
RRM_mean = mean(RRM_period)/fs;
MHR(loop) = 60/RRM_mean;

RRM_period_anot = diff(locs_R_mecg);
RRM_mean_anot = mean(RRM_period_anot)/fs;
MHR_anot(loop) = 60/RRM_mean_anot;

%% QRS Plotting
if cetak_qrs
figure
hold on
plot(ecg);
plot(locs_R,ecg(locs_R),'rv','MarkerFaceColor','r');
plot(locs_Q,ecg(locs_Q),'rs','MarkerFaceColor','g');
plot(locs_S,ecg(locs_S),'o','MarkerFaceColor','y');
hold off
legend('EKG Abdomen','R','Q','S');
title('Deteksi kompleks MQRS');
xlabel('Unit Sampel');
ylabel('Unit Amplitudo');
end

%% Segment Extraction and Averaging
w = round(mean(RRM_period));                                                % interval R-R untuk resample
ecg_seg = [];
%figure
%hold on
for i = 1:(N_pks-1)
    temp_seg = [];
    temp_seg=ecg(locs_R(i):locs_R(i+1));
    z(i)= length(temp_seg);                                                % interval masing-masing R-R
    ecg_seg(:,i)= resample(temp_seg,w,z(i)); 
%    plot(ecg_seg(:,i));
end
%hold off
%% Segmen and Mean Plotting
%Inisiasi ECG AVR
ecg_avr = [];

%TSA
for j = 1:w
    ecg_avr(j)= mean(ecg_seg(j,1:N_pks-1));
end

%% Reconstruction
mecg_temp = [];
mecg_temp = zeros([length(ecg) 1]);                                     % membuat array kosong sepanjang ecg

k=1;
if locs_R(1)~=1 && locs_R(1)<=(0.85*w)                                            % rekonstruksi mECG sebelum R-peak 1
        for j = w-(locs_R(1)-1):w
            mecg_temp(k) = ecg_avr(j);
            k = k+1;
        end
elseif locs_R(1)~=1 && locs_R(1)>(0.85*w)
        ecg_f2 = [];
        %ecg_f2 = resample(ecg_avr, (locs_R(1)-1),w);
        for j = 1:(locs_R(1)-1)
            mecg_temp(k) = ecg(j);
            k = k+1;
        end
end

for i = 1:N_pks-1
        ecg_f2 = [];
        ecg_f2 = resample(ecg_avr,(z(i)-1),w);                                 % melakukan resample balik untuk menyesuaikan durasi template dengan durasi masing-masing R-R
        for j = 1:(z(i)-1)
            mecg_temp(k) = ecg_f2(j);
            k = k+1;
        end
end
if locs_R(N_pks)<length(ecg) && length(ecg)-locs_R(N_pks)<=(0.85*w)                             % rekonstruksi mECG setelah R-peak terakhir
       for j = 1:(length(ecg)-locs_R(N_pks)-1)
           mecg_temp(k) = ecg_avr(j);
           k = k+1;
       end
elseif locs_R(N_pks)<length(ecg) && length(ecg)-locs_R(N_pks)>(0.85*w)
    ecg_f2 = [];
    %ecg_f2 = resample(ecg_avr, (length(ecg)-locs_R(N_pks)-1),w);
    for j = locs_R(N_pks):(length(ecg)-1)
        mecg_temp(k) = ecg(j);
        k = k+1;
        end
end
% replace qrs with original
for i=1:length(locs_R)
    mecg_temp(locs_Q(i):locs_S(i)) = ecg(locs_Q(i):locs_S(i));
end


        
%% Reconstruction Plotting
if cetak_rekons
figure
hold on
plot(ecg);
plot(mecg_temp,'linewidth',1.5); 
hold off
title('Sinyal MECG tiruan');
legend('Raw AECG','MECG Tiruan');
xlabel('Unit sampel');
ylabel('Unit amplitudo');
end

%% Subtracting
residu = [];
residu = ecg-mecg_temp;

%% Cek Orientasi [Pre-SNR]
cek_start = 1;
cek_end = 1000;
neg=0;
pos=0;
N_cek=5;

%disp([n_kanal, ': Subtraction done'])

if loop==120
    N_cek=3;
end

%Cek orientasi
for i=1:N_cek
    if loop==1&&i==N_cek
        cek_end=4950;
    end
    if abs(min(residu(cek_start:cek_end)))>max(residu(cek_start:cek_end))
        neg = neg+1;
    else pos = pos+1;
    end
    
cek_start=cek_start+1000;
cek_end=cek_end+1000;
end

if neg>pos
    residu=-residu;
end

%% Cetak Proses

t = 1/fs:1/fs:(length(ecg)/fs);
if cetak_proses
figure
subplot(3,1,1);
plot(t,ecg);
title('Sinyal Raw AECG');
xlabel('Waktu (s)');
ylabel('Amplitudo');
subplot(3,1,2);
plot(t,mecg_temp);
title('Sinyal MECG Tiruan');
xlabel('Waktu (s)');
ylabel('Amplitudo');
subplot(3,1,3);
plot(t,residu);
title('Sinyal FECG (Residu)');
xlabel('Waktu (s)');
ylabel('Amplitudo');
end

%% FECG Enhancement and Plotting
%fs = 500;
%f1=10;                                                                      % cuttoff low frequency to get rid of baseline wander
%f2=50;                                                                     % cuttoff frequency to discard high frequency noise
%Wn=[f1 f2]*2/fs;                                                           % cutt off based on fs
%N = 5;                                                                     % order of 3 less processing
%[a,b] = butter(N,Wn);                                                      % bandpass filtering
%ecg_h = filtfilt(a,b,residu);
 
%figure
%subplot(2,1,1);
%plot(residu);
%title('residu signal');
%subplot(2,1,2);
%plot(ecg_h);
%title('residu filtered (BPF[10 50]) signal');

%% Residu Channel Assignment

if kanal==1
    residu_1 = residu;
elseif kanal==2
    residu_2 = residu;
elseif kanal==3
    residu_3 = residu;
else
    residu_4 = residu;
end

%% QRS FECG Detection
[THR, R_val,locs_R_residu,delay]=pan_tompkin_fecg(residu,fs,0);

%% R anotation mathching
if cetak_preSNR
for i = 1:sum_RF-1
    if i==1 && locs_R_fecg(1)-5<=0
        [R_fecg_val(i),locs_R_temp]=max(residu(1:locs_R_fecg(i)+5));
        locs_R_fecg2(i)=locs_R_temp;
    elseif i==(sum_RF-1) && locs_R_fecg(sum_RF-1)+5>length(residu)
        [R_fecg_val(i),locs_R_temp]=max(residu(locs_R_fecg(i)-5:length(residu)));
        locs_R_fecg2(i)=locs_R_temp+locs_R_fecg(i)-5-1;
    else
    [R_fecg_val(i),locs_R_temp]=max(residu(locs_R_fecg(i)-5:locs_R_fecg(i)+5));
    locs_R_fecg2(i)=locs_R_temp+locs_R_fecg(i)-5-1;
    end
end


figure
%hold on
plot(residu);
%plot(locs_R_residu,residu(locs_R_residu),'o','MarkerFaceColor','g');
%plot(locs_R_fecg2,residu(locs_R_fecg2),'rv','MarkerFaceColor','r');
%hold off
title('Sinyal Residu');
xlabel('Unit Sampel');
ylabel('Unit Amplitudo');
end

%% SNR
%Inisiasi titik Q dan S
locs_Qf = [];
locs_Sf = [];
wP_Si=0;
wP_Ni=0;
sum_PSi=0;
sum_PNi=0;

ecg_inv = [];
ecg_inv = -residu;                                                        %ecg_h = data ecg hasil BPF
N_pks = length(locs_R_residu);
for i = 1:N_pks
     if i == 1 && (locs_R_residu(1)-round(0.03*fs))<=0
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(1:locs_R_residu(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R_residu(i):locs_R_residu(i)+round(0.03*fs)));  %0.06 = max durasi QR
            locs_Qf(i) = locs_Q_temp;
            locs_Sf(i) = locs_R_residu(i)+locs_S_temp-1;
     elseif i == 1 && (locs_R_residu(1)-round(0.03*fs))>0 
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R_residu(i)-round(0.03*fs):locs_R_residu(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R_residu(i):locs_R_residu(i)+round(0.03*fs)));
            locs_Qf(i) = locs_R_residu(i)-round(0.03*fs)+locs_Q_temp-1;
            locs_Sf(i) = locs_R_residu(i)+locs_S_temp-1;
     elseif i == N_pks && (locs_R_residu(N_pks)+round(0.03*fs))>length(ecg)
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R_residu(i)-round(0.03*fs):locs_R_residu(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R_residu(i):end));
            locs_Qf(i) = locs_R_residu(i)-round(0.03*fs)+locs_Q_temp-1;
            locs_Sf(i) = locs_R_residu(i)+locs_S_temp-1;
     elseif i == N_pks && (locs_R_residu(N_pks)+round(0.03*fs))<=length(ecg)
            [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R_residu(i)-round(0.03*fs):locs_R_residu(i)));
            [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R_residu(i):(locs_R_residu(N_pks)+round(0.03*fs))));
            locs_Qf(i) = locs_R_residu(i)-round(0.03*fs)+locs_Q_temp-1;
            locs_Sf(i) = locs_R_residu(i)+locs_S_temp-1;
     else                
          [Q_pks_temp,locs_Q_temp] = max(ecg_inv(locs_R_residu(i)-round(0.03*fs):locs_R_residu(i)));   
          [S_pks_temp,locs_S_temp] = max(ecg_inv(locs_R_residu(i):locs_R_residu(i)+round(0.03*fs)));
          locs_Qf(i) = locs_R_residu(i)-round(0.03*fs)+locs_Q_temp-1;
          locs_Sf(i) = locs_R_residu(i)+locs_S_temp-1; 
     end       
end

% Menghitung power EKGf
P_Si = [];
P_Si = residu;
for i=1:N_pks
    
    if i==1&&locs_Qf(1)>1
        P_Si(1:locs_Qf(1))=0;
    end
    
    if i==N_pks && locs_Sf(N_pks)<length(P_Si)
        P_Si(locs_Sf(N_pks):length(P_Si))=0;
    end
    
    if i<N_pks
    P_Si(locs_Sf(i):locs_Qf(i+1))=0;
    end
    
    wP_Si=wP_Si+(locs_Sf(i)-locs_Qf(i));
end

% Squaring
P_Si2 = [];
for j=50:length(P_Si)-50
    P_Si2(j)=P_Si(j)*P_Si(j);
    sum_PSi=sum_PSi+P_Si2(j);
end

% Menghitung power noise
P_Ni = [];
P_Ni = residu;
for i=1:N_pks
    P_Ni(locs_Qf(i):locs_Sf(i))=0;
    
    if i==1&&locs_Qf(1)>1
    wP_Ni=wP_Ni+locs_Qf(1)-1;
    end
    if i==N_pks&&locs_Sf(N_pks)<length(P_Ni)
        wP_Ni=wP_Ni+length(P_Ni)-locs_Sf(N_pks)-1;
    end
    
    if i<N_pks
        wP_Ni=wP_Ni+(locs_Qf(i+1)-locs_Sf(i));
    end
end
% Squaring
P_Ni2 = [];
for j=50:length(P_Ni)-50
    P_Ni2(j)=P_Ni(j)*P_Ni(j);
    sum_PNi=sum_PNi+P_Ni2(j);
end

PF=sum_PSi/wP_Si;
PN=sum_PNi/wP_Ni;
SUB=PF-PN;
if SUB<=0
    SNR=0;
else
SNR= 10*log((PF-PN)/PN);
end

if SNR==Inf || SNR==-Inf || length(locs_R_residu)<12
    SNR=0;
end

%validasi nilai SNR dengan puncak R terdeteksi
if kanal==1
    SNR_1(loop) = SNR;
    tSNR(1,1) = SNR*SNR;
    tkanal(1,:) = residu_1;
elseif kanal==2
    SNR_2(loop) = SNR;
    tSNR(1,2) = SNR*SNR;
    tkanal(2,:) = residu_2;
elseif kanal==3
    SNR_3(loop) = SNR;
    tSNR(1,3) = SNR*SNR;
    tkanal(3,:) = residu_3;
else
    SNR_4(loop) = SNR;
    tSNR(1,4) = SNR*SNR;
    tkanal(4,:) = residu_4;
end

if SNR<batas_SNR
    tSNR(1,kanal)=0;
end
sum_SNR = sum_SNR+tSNR(1,kanal);

%disp('SNR done')

if cetak_SNR
figure
subplot(2,1,1)
plot(P_Si2);
title('Power EKGf');
subplot(2,1,2)
plot(P_Ni2);
title('Power Noise');
end
%loop
%disp([n_kanal, ' done'])
else
   if kanal==1
    SNR_1(loop) = 0;
    tSNR(1,1) = 0;
    tkanal(1,:) = zeros(1,length(ecg));
elseif kanal==2
    SNR_2(loop) = 0;
    tSNR(1,2) = 0;
    tkanal(2,:) = zeros(1,length(ecg));
elseif kanal==3
    SNR_3(loop) = 0;
    tSNR(1,3) = 0;
    tkanal(3,:) = zeros(1,length(ecg));
else
    SNR_4(loop) = 0;
    tSNR(1,4) = 0;
    tkanal(4,:) = zeros(1,length(ecg));
   end 
%loop   
disp([n_kanal, ' mengandung artifak!'])
disp(string(loop))
%tSNR(1,kanal)
end
end


%% Data keempat kanal diambil

%Cek SNR pada looping
cek_SNR = 0;
for i=1:4
    if sqrt(tSNR(1,i))<=batas_SNR
        cek_SNR = cek_SNR+1;
    end
end

SNRf(loop) = (tSNR(1,1)/sum_SNR)*SNR_1(loop)+(tSNR(1,2)/sum_SNR)*SNR_2(loop)+(tSNR(1,3)/sum_SNR)*SNR_3(loop)+(tSNR(1,4)/sum_SNR)*SNR_4(loop);

if cek_SNR<4

kanal_add = [];

for i=1:length(residu)
    kanal_add(i) = (tSNR(1,1)/sum_SNR)*tkanal(1,i)+(tSNR(1,2)/sum_SNR)*tkanal(2,i)+(tSNR(1,3)/sum_SNR)*tkanal(3,i)+(tSNR(1,4)/sum_SNR)*tkanal(4,i);
end

residu = [];
residu = kanal_add;

if cetak_residu
figure
t=(dstart:1:dend)/(fs*60);
plot(residu);
title('Sinyal EKGf (Residu)');
xlabel('Waktu (menit)');
ylabel('Amplitudo');
end

MB=0;
MP=0;
FD=0;

cek_start = 1;
cek_end = 1000;
neg=0;
pos=0;
N_cek=5;

%disp([n_kanal, ': Subtraction done'])

if loop==120
    N_cek=3;
end

%Cek orientasi
for i=1:N_cek
    if loop==1&&i==N_cek
        cek_end=4950;
    end
    if abs(min(residu(cek_start:cek_end)))>max(residu(cek_start:cek_end))
        neg = neg+1;
    else pos = pos+1;
    end
    
cek_start=cek_start+1000;
cek_end=cek_end+1000;
end

if neg>pos
    residu=-residu;
end


%% QRS FECG Detection
[THR, R_val,locs_R_residu,delay]=pan_tompkin_fecg(residu,fs,0);

if length(locs_R_residu)>12
%% R peak validation [ perlu diperbaiki lagi. terutama untuk yang misplaced. kalo yang intervalnya terlalu panjang (misdetected sejauh ini aman) ]
RR_period_1 = [];
RR_period_1 = diff(locs_R_residu);
RR_period_n = [];
N_R = length(RR_period_1);
N_R_i = length(RR_period_1);
interval =0;
cek=0;
add=0;
a=1;

RR_med = median(RR_period_1);

for i=1:N_R
% Searching for missed beat
    if RR_period_1(i)>=1.3*RR_med
        add=round(RR_period_1(i)/RR_med);
        if mod(RR_period_1(i),RR_med)>=(0.5*RR_med)&&mod(RR_period_1(i),RR_med)<(0.6*RR_med)
            add=add-1;
            temp=RR_med;
        else
        temp=RR_period_1(i)/add;
        end
            for k=1:add
                RR_period_n(a)=temp;
                a=a+1;
            end
            MB=MB+(add-1);
% Searching for false detection 
    elseif i==N_R||i==1
        if RR_period_1(i)<0.9*RR_med
            FD=FD+1;
        elseif RR_period_1(i)>1.1*RR_med
            RR_period_n(a)=RR_med;
            a=a+1;
            FD=FD+1;
        else
            RR_period_n(a)=RR_period_1(i);
            a=a+1;
        end    
    elseif RR_period_1(i)>=(0.25*fs)
        RR_period_n(a)=RR_period_1(i);
        a=a+1;
    else
        FD = FD+1;
    end
end

N_R=a-1;

for i=1:N_R-1
 %Searching for misplaced beat
    if RR_period_n(i)<=(0.9*RR_med)&&RR_period_n(i+1)>=(1.1*RR_med)||RR_period_n(i)>=(1.1*RR_med)&&RR_period_n(i+1)<=(0.9*RR_med)
        RR_tot = RR_period_n(i)+RR_period_n(i+1);
        temp = RR_tot/2;
        RR_period_n(i) = temp;
        RR_period_n(i+1) = temp;
        MP = MP+1;
    end
end

for i=1:N_R
    % Searching for another missed beat after misplaced correction
    if RR_period_n(i)>=1.3*RR_med
        add=round(RR_period_n(i)/RR_med);
        if add==1
            RR_period_n(i)=RR_med;
        end
    MP=MP+1;
    end
end

%% FHR Calculation

RR_period_3(loop)= mean(RR_period_n)/fs;
FHR(loop) = (1/RR_period_3(loop))*60;
else
FHR(loop) = 0;
end
else
FHR(loop) = 0;
end

%if loop>1&&cek_FHR
%    if FHR(loop)>(25+FHR_var)||FHR(loop)<(FHR_var-25)
%        disp('Out of FHR Variability')
%        disp(string(loop))
%        FHR(loop)=0;
%    end
%end

%if FHR(loop)~=0
%    FHR_var=FHR(loop);
%    cek_FHR=1;
%end
        

%if length(locs_R_residu)<10
%    FHR(loop) = NaN;
%end

%% R anotation mathching
for i = 1:sum_RF-1
    if i==1 && locs_R_fecg(1)-5<=0
        [R_fecg_val(i),locs_R_temp]=max(residu(1:locs_R_fecg(i)+5));
        locs_R_fecg2(i)=locs_R_temp;
    elseif i==(sum_RF-1) && locs_R_fecg(sum_RF-1)+5>length(residu)
        [R_fecg_val(i),locs_R_temp]=max(residu(locs_R_fecg(i)-5:length(residu)));
        locs_R_fecg2(i)=locs_R_temp+locs_R_fecg(i)-5-1;
    else
    [R_fecg_val(i),locs_R_temp]=max(residu(locs_R_fecg(i)-5:locs_R_fecg(i)+5));
    locs_R_fecg2(i)=locs_R_temp+locs_R_fecg(i)-5-1;
    end
end

%% R FECG Plot 
if cetak
figure
hold on
plot(residu);
plot(locs_R_residu,residu(locs_R_residu),'o','MarkerFaceColor','g');
plot(locs_R_fecg2,residu(locs_R_fecg2),'rv','MarkerFaceColor','r');
hold off
title('Deteksi R FECG');
xlabel('Sampel');
ylabel('Amplitudo');
legend('FECG','Deteksi','Referensi');

chr = int2str(loop);
PlotName = ['\', id, '_', chr, '.png'];
saveas(gcf, [PathPLOT PlotName], 'png')
close
end

RR_period_2(loop) = mean(diff(locs_R_fecg))/fs;
FHR_anot(loop)= (1/RR_period_2(loop))*60;
tMB(loop)=MB;
tFD(loop)=FD;
tMP(loop)=MP;
N_Rpeak(loop)= length(locs_R_fecg);

dstart=dstart+5000;
dend=dend+5000;
end
