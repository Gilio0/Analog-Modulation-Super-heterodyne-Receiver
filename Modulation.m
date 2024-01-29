clear all;
clc;
%carrier frequencies
Fc_BBCArabic2 = 100000;
Fc_FM9090 = 155000;
Fc_QuranPalestine = 210000;
Fc_RussianVoice = 265000;
Fc_SkyNewsArabia = 320000;

% Specify the path to .wav files (input messages)
file_path_BBCArabic2 = 'E:\3rd Year Electronics & Electrical Communications\1st Term\labs\comm project\Short_BBCArabic2.wav';
file_path_FM9090 = 'E:\3rd Year Electronics & Electrical Communications\1st Term\labs\comm project\Short_FM9090.wav';
file_path_QuranPalestine = 'E:\3rd Year Electronics & Electrical Communications\1st Term\labs\comm project\Short_QuranPalestine.wav';
file_path_RussianVoice = 'E:\3rd Year Electronics & Electrical Communications\1st Term\labs\comm project\Short_RussianVoice.wav';
file_path_SkyNewsArabia = 'E:\3rd Year Electronics & Electrical Communications\1st Term\labs\comm project\Short_SkyNewsArabia.wav';

% Use audioread to read the .wav files (input messages)...(same fs for all)
[BBCArabic2, Fm1] = audioread(file_path_BBCArabic2);
[FM9090, fs] = audioread(file_path_FM9090);
[QuranPalestin, fs] = audioread(file_path_QuranPalestine);
[RussianVoice, fs] = audioread(file_path_RussianVoice);
[SkyNewsArabia, fs] = audioread(file_path_SkyNewsArabia);

% adding 2 channels to make single channel
BBCArabic2 = BBCArabic2(:,2) + BBCArabic2(:,1);
FM9090 = FM9090(:,2) + FM9090(:,1);
QuranPalestin = QuranPalestin(:,2) + QuranPalestin(:,1);
RussianVoice = RussianVoice(:,2) + RussianVoice(:,1);
SkyNewsArabia = SkyNewsArabia(:,2) + SkyNewsArabia(:,1);

%increase the sampling freq by factor & using interp function
 message_BBCArabic2 = interp(BBCArabic2,15);
 message_FM9090 = interp(FM9090,15);
 message_QuranPalestin = interp(QuranPalestin,15);
 message_RussianVoice = interp(RussianVoice,15);
 message_SkyNewsArabia = interp(SkyNewsArabia,15);
 fs = fs * 15;
 
% get the max length of all messages
Messages_lengthes = [length(message_BBCArabic2),length(message_FM9090),length(message_QuranPalestin),length(message_RussianVoice),length(message_SkyNewsArabia)];
maxMessageLength = max(Messages_lengthes);

% make all messages have the same length
numZerosToAddIn_FM9090 = (maxMessageLength - length(message_FM9090));
numZerosToAddIn_QuranPalestin = (maxMessageLength - length(message_QuranPalestin));
numZerosToAddIn_RussianVoice = (maxMessageLength - length(message_RussianVoice));
numZerosToAddIn_SkyNewsArabia = (maxMessageLength - length(message_SkyNewsArabia));
message_FM9090 = [message_FM9090; zeros(numZerosToAddIn_FM9090,1)];
message_QuranPalestin = [message_QuranPalestin; zeros(numZerosToAddIn_QuranPalestin,1)];
message_RussianVoice = [message_RussianVoice; zeros(numZerosToAddIn_RussianVoice,1)];
message_SkyNewsArabia = [message_SkyNewsArabia; zeros(numZerosToAddIn_SkyNewsArabia,1)];

% make fast fourier transform
BBCArabic2_spectrum = fft(message_BBCArabic2,maxMessageLength);
FM9090_spectrum = fft(message_FM9090,maxMessageLength);
QuranPalestin_spectrum = fft(message_QuranPalestin,maxMessageLength);
RussianVoice_spectrum = fft(message_RussianVoice,maxMessageLength);
SkyNewsArabia_spectrum = fft(message_SkyNewsArabia,maxMessageLength);

%plot the messages spectrum
f = fs * (0:(maxMessageLength/2))/maxMessageLength; % Frequency axis
magnitude_spectrum_BBCArabic2 = 2/maxMessageLength * abs(BBCArabic2_spectrum(1:maxMessageLength/2+1));
magnitude_spectrum_FM9090 = 2/maxMessageLength * abs(FM9090_spectrum(1:maxMessageLength/2+1));
magnitude_spectrum_QuranPalestin = 2/maxMessageLength * abs(QuranPalestin_spectrum(1:maxMessageLength/2+1));
magnitude_spectrum_RussianVoice = 2/maxMessageLength * abs(RussianVoice_spectrum(1:maxMessageLength/2+1));
magnitude_spectrum_SkyNewsArabia = 2/maxMessageLength * abs(SkyNewsArabia_spectrum(1:maxMessageLength/2+1));

figure;
xlabel('Frequency (Hertz)');
ylabel('Magnitude');
subplot(2,3,1);
plot(f, magnitude_spectrum_BBCArabic2);
title('BBCArabic2 Spectrum');
subplot(2,3,2);
plot(f, magnitude_spectrum_FM9090);
title('FM9090 Spectrum');
subplot(2,3,3);
plot(f, magnitude_spectrum_QuranPalestin);
title('QuranPalestin Spectrum');
subplot(2,3,4);
plot(f, magnitude_spectrum_RussianVoice);
title('RussianVoice Spectrum');
subplot(2,3,5); 
plot(f, magnitude_spectrum_SkyNewsArabia);
title('SkyNewsArabia Spectrum');

%AM_DSB_SC modulation using ammod function
dsb_sc_BBCArabic2  = (ammod(message_BBCArabic2,Fc_BBCArabic2,fs));
dsb_sc_FM9090  = (ammod(message_FM9090,Fc_FM9090,fs));
dsb_sc_QuranPalestin  = (ammod(message_QuranPalestin,Fc_QuranPalestine,fs));
dsb_sc_RussianVoice  = (ammod(message_RussianVoice,Fc_RussianVoice,fs));
dsb_sc_SkyNewsArabia  = (ammod(message_SkyNewsArabia,Fc_SkyNewsArabia,fs));
TDM = dsb_sc_BBCArabic2 + dsb_sc_FM9090+ dsb_sc_QuranPalestin + dsb_sc_RussianVoice + dsb_sc_SkyNewsArabia;

%plot the modulated signals spectrum
fft_of_modulated_BBCArabic2 = fft(dsb_sc_BBCArabic2);
fft_of_modulated_FM9090 = fft(dsb_sc_FM9090);
fft_of_modulated_QuranPalestin = fft(dsb_sc_QuranPalestin);
fft_of_modulated_RussianVoice = fft(dsb_sc_RussianVoice);
fft_of_modulated_SkyNewsArabia = fft(dsb_sc_SkyNewsArabia);

%Plotting FDM
f = linspace(-fs/2, fs/2, length(fft_of_modulated_BBCArabic2)); % Frequency vector(X-axis)
figure;
plot(f, fftshift(abs(fft_of_modulated_BBCArabic2)))
hold on;
plot(f, fftshift(abs(fft_of_modulated_FM9090)))
hold on;
plot(f, fftshift(abs(fft_of_modulated_QuranPalestin)))
hold on;
plot(f, fftshift(abs(fft_of_modulated_RussianVoice)))
hold on;
plot(f, fftshift(abs(fft_of_modulated_SkyNewsArabia)))
legend('BBCArabic2','FM9090','QuranPalestin','RussianVoice','SkyNewsArabia');
title('FDM');

%applying filters to get each message from the FDM
BBCArabic2_BandPassObject = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 82000, 85000, 115000, 118000, 60, 1, 60, fs);
BBCArabic2_BandPassFilter = design(BBCArabic2_BandPassObject);
filtered_BBCArabic2 = filter(BBCArabic2_BandPassFilter, TDM);
figure;
xlabel('Frequency (Hertz)');
ylabel('Magnitude');
plot(f,fftshift(abs(fft((filtered_BBCArabic2)))));
title('BBCArabic2 Spectrum after filtering from FDM');
grid on;

FM9090_BandPassObject = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 125000, 130000, 180000, 181500, 60, 1, 60, fs);
FM9090_BandPassFilter = design(FM9090_BandPassObject);
filtered_FM9090 = filter(FM9090_BandPassFilter, TDM);
figure;
xlabel('Frequency (Hertz)');
ylabel('Magnitude');
plot(f,fftshift(abs(fft((filtered_FM9090)))));
title('FM9090 Spectrum after filtering from FDM');
grid on;

QuranPalestin_BandPassObject = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 176800, 179000, 237800, 240900, 60, 1, 60, fs);
QuranPalestin_BandPassFilter = design(QuranPalestin_BandPassObject);
filtered_QuranPalestin = filter(QuranPalestin_BandPassFilter, TDM);
figure;
xlabel('Frequency (Hertz)');
ylabel('Magnitude');
plot(f,fftshift(abs(fft((filtered_QuranPalestin)))));
title('QuranPalestin Spectrum after filtering from FDM');
grid on;

RussianVoice_BandPassObject = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 231500, 232600, 300000, 305000, 60, 1, 60, fs);
RussianVoice_BandPassFilter = design(RussianVoice_BandPassObject);
filtered_RussianVoice = filter(RussianVoice_BandPassFilter, TDM);
figure;
xlabel('Frequency (Hertz)');
ylabel('Magnitude');
plot(f,fftshift(abs(fft((filtered_RussianVoice)))));
title('RussianVoice Spectrum after filtering from FDM');
grid on;

SkyNewsArabia_BandPassObject = fdesign.highpass('Fst,Fp,Ast,Ap', 300000, 305000, 60, 1, fs);
SkyNewsArabia_BandPassFilter = design(SkyNewsArabia_BandPassObject);
filtered_SkyNewsArabia = filter(SkyNewsArabia_BandPassFilter, TDM);
figure;
xlabel('Frequency (Hertz)');
ylabel('Magnitude');
plot(f,fftshift(abs(fft((filtered_SkyNewsArabia)))));
title('SkyNewsArabia Spectrum after filtering from FDM');
grid on;
%IF stage

IF = 27500; %Hz
IF_carrier_BBCArabic2 = Fc_BBCArabic2 + IF; % carrier frequency1
IF_carrier_FM9090 = Fc_FM9090 + IF; % carrier frequency2
IF_carrier_QuranPalestin = Fc_QuranPalestine + IF; % carrier frequency3
IF_carrier_RussianVoice = Fc_RussianVoice + IF; % carrier frequency4
IF_carrier_SkyNewsArabia = Fc_SkyNewsArabia + IF; % carrier frequency5

t_BBCArabic2 = (1:1:length(filtered_BBCArabic2))';
t_FM9090 = (1:1:length(filtered_FM9090))';
t_QuranPalestin = (1:1:length(filtered_QuranPalestin))';
t_RussianVoice = (1:1:length(filtered_RussianVoice))';
t_SkyNewsArabia = (1:1:length(filtered_SkyNewsArabia))';

carrier_BBCArabic2_IF = cos(2 * pi * IF_carrier_BBCArabic2 * t_BBCArabic2 * (1 / (fs)));
carrier_FM9090_IF2 = cos(2 * pi * IF_carrier_FM9090 * t_FM9090 * (1 / (fs)));
carrier_QuranPalestin_IF3 = cos(2 * pi * IF_carrier_QuranPalestin * t_QuranPalestin * (1 / (fs)));
carrier_RussianVoice_IF4 = cos(2 * pi * IF_carrier_RussianVoice * t_RussianVoice * (1 / (fs)));
carrier_SkyNewsArabia_IF5 = cos(2 * pi * IF_carrier_SkyNewsArabia * t_SkyNewsArabia * (1 / (fs)));

carrier_BBCArabic2_IF(end + length(filtered_BBCArabic2) - length(carrier_BBCArabic2_IF), 1) = 0;
carrier_FM9090_IF2(end + length(filtered_FM9090) - length(carrier_FM9090_IF2), 1) = 0;
carrier_QuranPalestin_IF3(end + length(filtered_QuranPalestin) - length(carrier_QuranPalestin_IF3), 1) = 0;
carrier_RussianVoice_IF4(end + length(filtered_RussianVoice) - length(carrier_RussianVoice_IF4), 1) = 0;
carrier_SkyNewsArabia_IF5(end + length(filtered_SkyNewsArabia) - length(carrier_SkyNewsArabia_IF5), 1) = 0;

IF_received_BBCArabic2 = filtered_BBCArabic2 .* carrier_BBCArabic2_IF;
IF_received_FM9090 = filtered_FM9090 .* carrier_FM9090_IF2;
IF_received_QuranPalestin = filtered_QuranPalestin .* carrier_QuranPalestin_IF3;
IF_received_RussianVoice = filtered_RussianVoice .* carrier_RussianVoice_IF4;
IF_received_SkyNewsArabia = filtered_SkyNewsArabia .* carrier_SkyNewsArabia_IF5;

%prepearing for spectrum plotting
spectrum_IF_RECEIVED_BBCArabic2 = fftshift(fft(IF_received_BBCArabic2));
spectrum_IF_RECEIVED_FM9090 = fftshift(fft(IF_received_FM9090));
spectrum_IF_RECEIVED_QuranPalestin = fftshift(fft(IF_received_QuranPalestin));
spectrum_IF_RECEIVED_RussianVoice = fftshift(fft(IF_received_RussianVoice));
spectrum_IF_RECEIVED_SkyNewsArabia = fftshift(fft(IF_received_SkyNewsArabia));

%adjusting X-axis for each then plotting each message after IF oscillator
f_RECEIVED_BBCArabic2 = (-length(spectrum_IF_RECEIVED_BBCArabic2) / 2:1:length(spectrum_IF_RECEIVED_BBCArabic2) / 2 - 1)';
f_RECEIVED_FM9090 = (-length(spectrum_IF_RECEIVED_FM9090) / 2:1:length(spectrum_IF_RECEIVED_FM9090) / 2 - 1)';
f_RECEIVED_QuranPalestin = (-length(spectrum_IF_RECEIVED_QuranPalestin) / 2:1:length(spectrum_IF_RECEIVED_QuranPalestin) / 2 - 1)';
f_RECEIVED_RussianVoice = (-length(spectrum_IF_RECEIVED_RussianVoice) / 2:1:length(spectrum_IF_RECEIVED_RussianVoice) / 2 - 1)';
f_RECEIVED_SkyNewsArabia = (-length(spectrum_IF_RECEIVED_SkyNewsArabia) / 2:1:length(spectrum_IF_RECEIVED_SkyNewsArabia) / 2 - 1)';

figure;
plot(f_RECEIVED_BBCArabic2 * fs / length(spectrum_IF_RECEIVED_BBCArabic2), abs(spectrum_IF_RECEIVED_BBCArabic2));
title("BBCArabic2 after multiplying by IF carrier");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

figure;
plot(f_RECEIVED_FM9090 * fs / length(spectrum_IF_RECEIVED_FM9090), abs(spectrum_IF_RECEIVED_FM9090));
title("FM9090 after multiplying by IF carrier");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

figure;
plot(f_RECEIVED_QuranPalestin * fs / length(spectrum_IF_RECEIVED_QuranPalestin), abs(spectrum_IF_RECEIVED_QuranPalestin));
title("QuranPalestin after multiplying by IF carrier");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

figure;
plot(f_RECEIVED_RussianVoice * fs / length(spectrum_IF_RECEIVED_RussianVoice), abs(spectrum_IF_RECEIVED_RussianVoice));
title("RussianVoice after multiplying by IF carrier");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

figure;
plot(f_RECEIVED_SkyNewsArabia * fs / length(spectrum_IF_RECEIVED_SkyNewsArabia), abs(spectrum_IF_RECEIVED_SkyNewsArabia));
title("SkyNewsArabia after multiplying by IF carrier");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

%% IF stage BPF

BBCArabic2_BandPassObject_IFStage = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 4000, 5000, 50000, 55000, 60, 1, 60, fs);
BBCArabic2_BandPassFilter_IfStage = design(BBCArabic2_BandPassObject_IFStage,'equiripple');
filtered_BBCArabic2_IFStage = filter(BBCArabic2_BandPassFilter_IfStage, IF_received_BBCArabic2);
spectrum_filtered_BBCArabic2_IFStage = fftshift(abs(fft(filtered_BBCArabic2_IFStage)));
f_RECEIVED_BBCArabic2_BPF = (-length(spectrum_filtered_BBCArabic2_IFStage) / 2:1:length(spectrum_filtered_BBCArabic2_IFStage) / 2 - 1)';
figure;
plot(f_RECEIVED_BBCArabic2_BPF * fs / length(spectrum_filtered_BBCArabic2_IFStage),abs(spectrum_filtered_BBCArabic2_IFStage),'r');
title("BBCArabic2 after IF stage BPF");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

FM9090_BandPassObject_IFStage = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 4000, 5000, 50000, 55000, 60, 1, 60, fs);
FM9090_BandPassFilter_IfStage = design(FM9090_BandPassObject_IFStage,'equiripple');
filtered_FM9090_IFStage = filter(FM9090_BandPassFilter_IfStage, IF_received_FM9090);
spectrum_filtered_FM9090_IFStage = fftshift(abs(fft(filtered_FM9090_IFStage)));
f_RECEIVED_FM9090_BPF = (-length(spectrum_filtered_FM9090_IFStage) / 2:1:length(spectrum_filtered_FM9090_IFStage) / 2 - 1)';
figure;
plot(f_RECEIVED_FM9090_BPF * fs / length(spectrum_filtered_FM9090_IFStage),abs(spectrum_filtered_FM9090_IFStage),'r');
title("FM9090 after IF stage BPF");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

QuranPalestin_BandPassObject_IFStage = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 4000, 5000, 50000, 55000, 60, 1, 60, fs);
QuranPalestin_BandPassFilter_IfStage = design(QuranPalestin_BandPassObject_IFStage,'equiripple');
filtered_QuranPalestin_IFStage = filter(QuranPalestin_BandPassFilter_IfStage, IF_received_QuranPalestin);
spectrum_filtered_QuranPalestin_IFStage = fftshift(abs(fft(filtered_QuranPalestin_IFStage)));
f_RECEIVED_QuranPalestin_BPF = (-length(spectrum_filtered_QuranPalestin_IFStage) / 2:1:length(spectrum_filtered_QuranPalestin_IFStage) / 2 - 1)';
figure;
plot(f_RECEIVED_QuranPalestin_BPF * fs / length(spectrum_filtered_QuranPalestin_IFStage),abs(spectrum_filtered_QuranPalestin_IFStage),'r');
title("QuranPalestin after IF stage BPF");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

RussianVoice_BandPassObject_IFStage = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 4000, 5000, 50000, 55000, 60, 1, 60, fs);
RussianVoice_BandPassFilter_IfStage = design(RussianVoice_BandPassObject_IFStage,'equiripple');
filtered_RussianVoice_IFStage = filter(RussianVoice_BandPassFilter_IfStage, IF_received_RussianVoice);
spectrum_filtered_RussianVoice_IFStage = fftshift(abs(fft(filtered_RussianVoice_IFStage)));
f_RECEIVED_RussianVoice_BPF = (-length(spectrum_filtered_RussianVoice_IFStage) / 2:1:length(spectrum_filtered_RussianVoice_IFStage) / 2 - 1)';
figure;
plot(f_RECEIVED_RussianVoice_BPF * fs / length(spectrum_filtered_RussianVoice_IFStage),abs(spectrum_filtered_RussianVoice_IFStage),'r');
title("RussianVoice after IF stage BPF");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

SkyNewsArabia_BandPassObject_IFStage = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 16370, 16720, 50000, 55000, 60, 1, 60, fs);
SkyNewsArabia_BandPassFilter_IfStage = design(SkyNewsArabia_BandPassObject_IFStage,'equiripple');
filtered_SkyNewsArabia_IFStage = filter(SkyNewsArabia_BandPassFilter_IfStage, IF_received_SkyNewsArabia);
spectrum_filtered_SkyNewsArabia_IFStage = fftshift(abs(fft(filtered_SkyNewsArabia_IFStage)));
f_RECEIVED_SkyNewsArabia_BPF = (-length(spectrum_filtered_SkyNewsArabia_IFStage) / 2:1:length(spectrum_filtered_SkyNewsArabia_IFStage) / 2 - 1)';
figure;
plot(f_RECEIVED_SkyNewsArabia_BPF * fs / length(spectrum_filtered_SkyNewsArabia_IFStage),abs(spectrum_filtered_SkyNewsArabia_IFStage),'r');
title("SkyNewsArabia after IF stage BPF");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
%% back to base band
T_BBCArabic2 = (1:1:length(filtered_BBCArabic2_IFStage))';
T_FM9090 = (1:1:length(filtered_FM9090_IFStage))';
T_QuranPalestin = (1:1:length(filtered_QuranPalestin_IFStage))';
T_RussianVoice = (1:1:length(filtered_RussianVoice_IFStage))';
T_SkyNewsArabia = (1:1:length(filtered_SkyNewsArabia_IFStage))';

carrier_BBCArabic2_Base_Band = cos(2 * pi * IF * T_BBCArabic2 * (1 / (fs)));
carrier_BBCArabic2_Base_Band(end + length(filtered_BBCArabic2) - length(carrier_BBCArabic2_Base_Band), 1) = 0;
Base_Band_received_BBCArabic2 = filtered_BBCArabic2_IFStage .* carrier_BBCArabic2_Base_Band;
spectrum_Base_Band_received_BBCArabic2 = fftshift(fft(Base_Band_received_BBCArabic2));
f_BBCArabic2_BASE_BAND = (-length(spectrum_Base_Band_received_BBCArabic2) / 2:1:length(spectrum_Base_Band_received_BBCArabic2) / 2 - 1)';

figure;
plot(f_BBCArabic2_BASE_BAND * fs / length(spectrum_Base_Band_received_BBCArabic2), abs(spectrum_Base_Band_received_BBCArabic2));
title("Base Band stage of BBCArabic2");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

BBCArabic2_lowPassObjectBaseBandFilter = fdesign.lowpass('Fp,Fst,Ap,Ast',20000,30000,1,80,fs);
BBCArabic2_lowPassFilterBaseBandStage = design(BBCArabic2_lowPassObjectBaseBandFilter,'equiripple');
filtered_BBCArabic2_BaseBandStage = filter(BBCArabic2_lowPassFilterBaseBandStage, Base_Band_received_BBCArabic2);
spectrum_filtered_BBCArabic2_BaseBandStage = fftshift(abs(fft(filtered_BBCArabic2_BaseBandStage)));

figure;
plot(f_BBCArabic2_BASE_BAND * fs / length(spectrum_Base_Band_received_BBCArabic2), abs(spectrum_filtered_BBCArabic2_BaseBandStage));
title("filtered Base Band stage of BBCArabic2");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
% Play the audio
Demodulated_BBCArabic2 = downsample(filtered_BBCArabic2_BaseBandStage,15);
audiowrite('Demodulated_BBCArabic2.wav',Demodulated_BBCArabic2,fs/15);


carrier_FM9090_Base_Band = cos(2 * pi * IF * T_FM9090 * (1 / (fs)));
carrier_FM9090_Base_Band(end + length(filtered_FM9090) - length(carrier_FM9090_Base_Band), 1) = 0;
Base_Band_received_FM9090 = filtered_FM9090_IFStage .* carrier_FM9090_Base_Band;
spectrum_Base_Band_received_FM9090 = fftshift(fft(Base_Band_received_FM9090));
f_FM9090_BASE_BAND = (-length(spectrum_Base_Band_received_FM9090) / 2:1:length(spectrum_Base_Band_received_FM9090) / 2 - 1)';

figure;
plot(f_FM9090_BASE_BAND * fs / length(spectrum_Base_Band_received_FM9090), abs(spectrum_Base_Band_received_FM9090));
title("Base Band stage of FM9090");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

FM9090_lowPassObjectBaseBandFilter = fdesign.lowpass('Fp,Fst,Ap,Ast',20000,30000,1,80,fs);
FM9090_lowPassFilterBaseBandStage = design(FM9090_lowPassObjectBaseBandFilter,'equiripple');
filtered_FM9090_BaseBandStage = filter(FM9090_lowPassFilterBaseBandStage, Base_Band_received_FM9090);
spectrum_filtered_FM9090_BaseBandStage = fftshift(abs(fft(filtered_FM9090_BaseBandStage)));

figure;
plot(f_FM9090_BASE_BAND * fs / length(spectrum_Base_Band_received_FM9090), abs(spectrum_filtered_FM9090_BaseBandStage));
title("filtered Base Band stage of FM9090");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
% Play the audio
Demodulated_FM9090 = downsample(filtered_FM9090_BaseBandStage,15);
audiowrite('Demodulated_FM9090.wav',Demodulated_FM9090,fs/15);


carrier_QuranPalestin_Base_Band = cos(2 * pi * IF * T_QuranPalestin * (1 / (fs)));
carrier_QuranPalestin_Base_Band(end + length(filtered_QuranPalestin) - length(carrier_QuranPalestin_Base_Band), 1) = 0;
Base_Band_received_QuranPalestin = filtered_QuranPalestin_IFStage .* carrier_QuranPalestin_Base_Band;
spectrum_Base_Band_received_QuranPalestin = fftshift(fft(Base_Band_received_QuranPalestin));
f_QuranPalestin_BASE_BAND = (-length(spectrum_Base_Band_received_QuranPalestin) / 2:1:length(spectrum_Base_Band_received_QuranPalestin) / 2 - 1)';

figure;
plot(f_QuranPalestin_BASE_BAND * fs / length(spectrum_Base_Band_received_QuranPalestin), abs(spectrum_Base_Band_received_QuranPalestin));
title("Base Band stage of QuranPalestin");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

QuranPalestin_lowPassObjectBaseBandFilter = fdesign.lowpass('Fp,Fst,Ap,Ast',20000,30000,1,80,fs);
QuranPalestin_lowPassFilterBaseBandStage = design(QuranPalestin_lowPassObjectBaseBandFilter,'equiripple');
filtered_QuranPalestin_BaseBandStage = filter(QuranPalestin_lowPassFilterBaseBandStage, Base_Band_received_QuranPalestin);
spectrum_filtered_QuranPalestin_BaseBandStage = fftshift(abs(fft(filtered_QuranPalestin_BaseBandStage)));

figure;
plot(f_QuranPalestin_BASE_BAND * fs / length(spectrum_Base_Band_received_QuranPalestin), abs(spectrum_filtered_QuranPalestin_BaseBandStage));
title("filtered Base Band stage of QuranPalestin");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
% Play the audio
Demodulated_QuranPalestin = downsample(filtered_QuranPalestin_BaseBandStage,15);
audiowrite('Demodulated_QuranPalestin.wav',Demodulated_QuranPalestin,fs/15);

carrier_RussianVoice_Base_Band = cos(2 * pi * IF * T_RussianVoice * (1 / (fs)));
carrier_RussianVoice_Base_Band(end + length(filtered_RussianVoice) - length(carrier_RussianVoice_Base_Band), 1) = 0;
Base_Band_received_RussianVoice = filtered_RussianVoice_IFStage .* carrier_RussianVoice_Base_Band;
spectrum_Base_Band_received_RussianVoice = fftshift(fft(Base_Band_received_RussianVoice));
f_RussianVoice_BASE_BAND = (-length(spectrum_Base_Band_received_RussianVoice) / 2:1:length(spectrum_Base_Band_received_RussianVoice) / 2 - 1)';

figure;
plot(f_RussianVoice_BASE_BAND * fs / length(spectrum_Base_Band_received_RussianVoice), abs(spectrum_Base_Band_received_RussianVoice));
title("Base Band stage of RussianVoice");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

RussianVoice_lowPassObjectBaseBandFilter = fdesign.lowpass('Fp,Fst,Ap,Ast',20000,30000,1,80,fs);
RussianVoice_lowPassFilterBaseBandStage = design(RussianVoice_lowPassObjectBaseBandFilter,'equiripple');
filtered_RussianVoice_BaseBandStage = filter(RussianVoice_lowPassFilterBaseBandStage, Base_Band_received_RussianVoice);
spectrum_filtered_RussianVoice_BaseBandStage = fftshift(abs(fft(filtered_RussianVoice_BaseBandStage)));

figure;
plot(f_RussianVoice_BASE_BAND * fs / length(spectrum_Base_Band_received_RussianVoice), abs(spectrum_filtered_RussianVoice_BaseBandStage));
title("filtered Base Band stage of RussianVoice");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
% Play the audio
Demodulated_RussianVoice = downsample(filtered_RussianVoice_BaseBandStage,15);
audiowrite('Demodulated_RussianVoice.wav',Demodulated_RussianVoice,fs/15);

carrier_SkyNewsArabia_Base_Band = cos(2 * pi * IF * T_SkyNewsArabia * (1 / (fs)));
carrier_SkyNewsArabia_Base_Band(end + length(filtered_SkyNewsArabia) - length(carrier_SkyNewsArabia_Base_Band), 1) = 0;
Base_Band_received_SkyNewsArabia = filtered_SkyNewsArabia_IFStage .* carrier_SkyNewsArabia_Base_Band;
spectrum_Base_Band_received_SkyNewsArabia = fftshift(fft(Base_Band_received_SkyNewsArabia));
f_SkyNewsArabia_BASE_BAND = (-length(spectrum_Base_Band_received_SkyNewsArabia) / 2:1:length(spectrum_Base_Band_received_SkyNewsArabia) / 2 - 1)';

figure;
plot(f_SkyNewsArabia_BASE_BAND * fs / length(spectrum_Base_Band_received_SkyNewsArabia), abs(spectrum_Base_Band_received_SkyNewsArabia));
title("Base Band stage of SkyNewsArabia");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

SkyNewsArabia_lowPassObjectBaseBandFilter = fdesign.lowpass('Fp,Fst,Ap,Ast',20000,30000,1,80,fs);
SkyNewsArabia_lowPassFilterBaseBandStage = design(SkyNewsArabia_lowPassObjectBaseBandFilter,'equiripple');
filtered_SkyNewsArabia_BaseBandStage = filter(SkyNewsArabia_lowPassFilterBaseBandStage, Base_Band_received_SkyNewsArabia);
spectrum_filtered_SkyNewsArabia_BaseBandStage = fftshift(abs(fft(filtered_SkyNewsArabia_BaseBandStage)));

figure;
plot(f_SkyNewsArabia_BASE_BAND * fs / length(spectrum_Base_Band_received_SkyNewsArabia), abs(spectrum_filtered_SkyNewsArabia_BaseBandStage));
title("filtered Base Band stage of SkyNewsArabia");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
% Play the audio
Demodulated_SkyNewsArabia = downsample(filtered_SkyNewsArabia_BaseBandStage,15);
audiowrite('Demodulated_SkyNewsArabia.wav',Demodulated_SkyNewsArabia,fs/15);