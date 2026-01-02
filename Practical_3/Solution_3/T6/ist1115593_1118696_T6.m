%N = 1000; %Number of Bits
rb = 1e3; 
fs = 100 * rb;
% for b. T_sim = 10e-3; % 10ms
T_sim = 1;
N = rb * T_sim;
oversample_factor = fs/rb;

binary_data = randi([0,1],1, N);

V = 1; 
polar_nrz = (2 * binary_data - 1) * V; 
zoh_signal = repelem(polar_nrz, oversample_factor); 

% Add AWGN (Noise) with power sigma^2 = 0.2 W
sigma2 = 0.2; 
sigma = sqrt(sigma2); 
noisy_signal = zoh_signal + sigma2 * randn(size(zoh_signal)); 

% Sampling at Receiver for 5.1 A
%sampled_signal = noisy_signal(oversample_factor:oversample_factor:end); % Sample once per bit

%Match Filter For 5.1 B
matched_filter = flip(polar_nrz);
matched_filter = repelem(matched_filter, oversample_factor);
filtered_signal = filter(matched_filter, 1, noisy_signal);
sampled_signal = filtered_signal(oversample_factor:oversample_factor:end);


detected_data = sampled_signal > 0; % Decision threshold at 0


errors = sum(binary_data ~= detected_data); 
ber = errors / N; % Bit Error Rate


time = (0:length(zoh_signal)-1) / fs; 
time_detected = (0:N-1) / rb; 
time_filter = (0:length(filtered_signal)-1) / fs;
% Plot Waveforms
figure;

% Transmitted Signal (Polar NRZ)
subplot(4, 1, 1);
plot(time, zoh_signal, 'LineWidth', 1);
title('Transmitted Signal (Polar NRZ)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Noisy Signal (After AWGN)
subplot(4, 1, 2);
plot(time, noisy_signal, 'LineWidth', 1);
title('Received Signal (After AWGN, \sigma^2 = 0.2 W)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Filtered Signal (After Matched Filter)
subplot(4, 1, 3);
plot(time_filter, filtered_signal, 'LineWidth', 1);
title('Filtered Signal (After Matched Filter)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Detected Signal (After Sampling)
subplot(4, 1, 4);
stem(time_detected, detected_data, 'LineWidth', 1);
title('Detected Signal');
xlabel('Time (s)');
ylabel('Binary Data');
grid on;





figure;
bit_start = 1; 
num_bits_to_plot = 10; 
samples_to_plot = num_bits_to_plot * oversample_factor;

time_zoom = (0:samples_to_plot-1) / fs;

plot(time_zoom, zoh_signal(1:samples_to_plot), 'b', 'LineWidth', 1.2); hold on;


sample_times = (0:num_bits_to_plot-1) / rb; 
sample_values = sampled_signal(1:num_bits_to_plot); 
stem(sample_times, sample_values, 'r', 'LineWidth', 1.5);

title('Detailed Observation of ZOH and Sampling');
xlabel('Time (s)');
ylabel('Amplitude');
legend('ZOH Signal', 'Sampled Signal');
grid on;


errors = sum(binary_data ~= detected_data); 
ber_value = errors / N;
display(ber_value);

figure;
eyediagram(noisy_signal, oversample_factor/10); 
title('Eye Diagram for \sigma^2 = 0.2');