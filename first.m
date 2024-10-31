clc;
clearvars;
close all;

%---------------------------------------------------------------%

                         % Parameters

filename = 'pulse.txt'; 
delimiterIn = ' ';       % Delimiter (space)
headerlinesIn = 1;       % Number of header lines
Fs = 100;                % Sampling frequency (100 Hz)
max_pulse_rate_hz = 3;   % Maximum pulse rate in Hz (180 bpm = 3 Hz)
min_breath_rate_hz = 0.5; % Minimum frequency to preserve (remove below 0.5 Hz)

%---------------------------------------------------------------%
                      % Load the data

Data_struct = importdata(filename, delimiterIn, headerlinesIn);
Led_R = Data_struct.data(:, 1);  % Red LED data (first column)
Led_IR = Data_struct.data(:, 2); % Infrared LED data (second column)

%---------------------------------------------------------------%
                  % Skip the first 10 seconds

skip_samples = 10 * Fs;  % Number of samples to skip
Led_R = Led_R(skip_samples + 1:end);  % Red signal after 10 seconds
Led_IR = Led_IR(skip_samples + 1:end); % Infrared signal after 10 seconds

                % Select a 30-second observation window
window_samples = 30 * Fs;  % Number of samples in 30 seconds
Led_R = Led_R(1:window_samples);  %  first 30 seconds of Red signal
Led_IR = Led_IR(1:window_samples); %  first 30 seconds of Infrared signal                 
time = (0:window_samples - 1) / Fs;

%---------------------------------------------------------------%
        %  Low-pass filtering (in the frequency domain)

                     %  FFT of both signals
Led_R_fft = fft(Led_R);
Led_IR_fft = fft(Led_IR);

                      % Frequency vector
frequencies = (0:window_samples - 1) * (Fs / window_samples);

         %  low-pass filter:  frequencies above 3 Hz to zero

cutoff_idx = find(frequencies > max_pulse_rate_hz, 1); %  cutoff index
Led_R_fft(cutoff_idx:end - cutoff_idx + 1) = 0;  
Led_IR_fft(cutoff_idx:end - cutoff_idx + 1) = 0;  

       % Inverse FFT to get the filtered signals back in time domain
Led_R_filtered = real(ifft(Led_R_fft));
Led_IR_filtered = real(ifft(Led_IR_fft));

%---------------------------------------------------------------%
              %  Plots of original and filtered signals

figure;

                  %  Original Red Signal
subplot(2, 2, 1);
plot(time, Led_R, 'r');
title('Original Red LED Signal (30 seconds)');
xlabel('Time (seconds)');
ylabel('Amplitude');

              %  Filtered Red Signal (Low-pass)
subplot(2, 2, 2);
plot(time, Led_R_filtered, 'r');
title('Filtered Red LED Signal (Low-pass at 3 Hz)');
xlabel('Time (seconds)');
ylabel('Amplitude');

               %  Original Infrared Signal
subplot(2, 2, 3);
plot(time, Led_IR, 'b');
title('Original Infrared LED Signal (30 seconds)');
xlabel('Time (seconds)');
ylabel('Amplitude');

            %  Filtered Infrared Signal (Low-pass)
subplot(2, 2, 4);
plot(time, Led_IR_filtered, 'b');
title('Filtered Infrared LED Signal (Low-pass at 3 Hz)');
xlabel('Time (seconds)');
ylabel('Amplitude');
%---------------------------------------------------------------%

  %  High-pass filtering to remove frequencies below 0.5 Hz for Red signal

%  high-pass filter: frequencies below 0.5 Hz to zero for Red signal
cutoff_idx_highpass = find(frequencies < min_breath_rate_hz, 1, 'last'); 
Led_R_fft(1:cutoff_idx_highpass) = 0;  % Zero out low frequencies for Red signal

% Inverse FFT to get the high-pass filtered signal back in time domain
Led_R_highpass = real(ifft(Led_R_fft));

            %  plot for High-pass Filtered Red Signal
figure;
plot(time, Led_R_highpass, 'r');
title('High-pass Filtered Red Signal (Above 0.5 Hz)');
xlabel('Time (seconds)');
ylabel('Amplitude');
%---------------------------------------------------------------%
         %   pulse rate from the filtered Red signal


        % detect the peaks of the high-pass filtered signal

[pks, locs] = findpeaks(Led_R_highpass, time);

                %  the time differences between peaks

time_diffs = diff(locs);  

                  %  the pulse rate (BPM)

mean_peak_interval = mean(time_diffs);  % Average time between peaks
pulse_rate_bpm = 60 / mean_peak_interval;  % Convert to beats per minute
%---------------------------------------------------------------%
                          %  Interpolation

       % peaks (max) and valleys (min) for Red and Infrared signals
[pks_R, locs_R] = findpeaks(Led_R_filtered, time);  % Maxima of Red signal
[pks_R_min, locs_R_min] = findpeaks(-Led_R_filtered, time);  % Minima of Red signal (inverted peaks)
pks_R_min = -pks_R_min;  % Restoring the negative peaks to their original value

[pks_IR, locs_IR] = findpeaks(Led_IR_filtered, time);  % Maxima of Infrared signal
[pks_IR_min, locs_IR_min] = findpeaks(-Led_IR_filtered, time);  % Minima of Infrared signal (inverted peaks)
pks_IR_min = -pks_IR_min;  
               % Interpolate using 'spline' method
interp_R_max = interp1(locs_R, pks_R, time, 'spline');
interp_R_min = interp1(locs_R_min, pks_R_min, time, 'spline');
interp_IR_max = interp1(locs_IR, pks_IR, time, 'spline');
interp_IR_min = interp1(locs_IR_min, pks_IR_min, time, 'spline');

          %  R and SaO2 using interpolated curves
I_AC_R = (interp_R_max - interp_R_min) / 2;
I_AC_IR = (interp_IR_max - interp_IR_min) / 2;

I_DC_R = (interp_R_max + interp_R_min) / 2;
I_DC_IR = (interp_IR_max + interp_IR_min) / 2;

             %  R for each sample
R_values = (I_AC_R ./ I_DC_R) ./ (I_AC_IR ./ I_DC_IR);

            %   mean R value
R_mean = mean(R_values);

             % SaO2 using the given formula
SaO2 = 110 - 25 * R_mean;

              %  plot for Interpolation and SaO2
figure;
subplot(2, 1, 1);
plot(time, Led_R_filtered, 'b'); hold on;
plot(time, interp_R_max, 'r', 'LineWidth', 1.5);  % Max curve in red
plot(time, interp_R_min, 'r', 'LineWidth', 1.5);  % Min curve in red
title(sprintf('Red LED Signal (Filtered) with Interpolation, SaO_2 = %.2f%%', SaO2));
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('Filtered Signal', 'Interpolated Max', 'Interpolated Min');

             %  Infrared Signal with Interpolation
subplot(2, 1, 2);
plot(time, Led_IR_filtered, 'b'); hold on;
plot(time, interp_IR_max, 'r', 'LineWidth', 1.5);  % Max curve in red
plot(time, interp_IR_min, 'r', 'LineWidth', 1.5);  % Min curve in red
title(sprintf('Infrared LED Signal (Filtered) with Interpolation, SaO_2 = %.2f%%', SaO2));
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('Filtered Signal', 'Interpolated Max', 'Interpolated Min');

                    % Display SaO2
fprintf('Estimated SaO2: %.2f%%\n', SaO2);

            %  plot for Pulse Rate based on peaks
figure;
plot(time, Led_R_highpass, 'b'); hold on;
plot(locs, pks, 'rv', 'MarkerFaceColor', 'r');  
title(sprintf('Filtered Red Signal - BPM = %.2f', pulse_rate_bpm));
xlabel('Time (s)');
ylabel('Amplitude');
legend('Filtered Red Signal', 'Detected Peaks');

%---------------------------------------------------------------%
          % Alternative Pulse Rate Computation using FFT

       %  the FFT of the high-pass filtered Red signal
Led_R_highpass_fft = fft(Led_R_highpass);

               %  the magnitude of the FFT
fft_magnitude = abs(Led_R_highpass_fft);

                % Frequency vector in BPM (60 * f)
frequencies_bpm = frequencies * 60;

    %  frequency interval between 60 and 90 BPM (1 Hz and 1.5 Hz)
focus_idx = (frequencies_bpm >= 60 & frequencies_bpm <= 90);

        %  relevant part of the FFT and frequency vectors
fft_magnitude_focus = fft_magnitude(focus_idx);
frequencies_bpm_focus = frequencies_bpm(focus_idx);

          % peak in that range (corresponding to pulse rate)
[peak_value, peak_idx] = max(fft_magnitude_focus);
pulse_rate_bpm_fft = frequencies_bpm_focus(peak_idx);

          % Plot of the FFT magnitude within the 60-90 BPM range
figure;
plot(frequencies_bpm_focus, fft_magnitude_focus);
title(sprintf('Filtered Red Signal - Frequency axis - BPM = %.2f', pulse_rate_bpm_fft));
xlabel('BPM');
ylabel('Magnitude');

                % Display the pulse rate from FFT
fprintf('Estimated pulse rate from FFT: %.2f BPM\n', pulse_rate_bpm_fft);

%---------------------------------------------------------------%
