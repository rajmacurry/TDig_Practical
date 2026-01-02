% Parameters
T_b = 1; 
V = 1;
data = [0, 1, 0, 0, 1, 1];
polar_nrz = (2 * data - 1) * V; 
N = length(data);
fs = 100; 
dt = 1 / fs;
t = 0:dt:(N * T_b - dt); 


zoh_signal = repelem(polar_nrz, fs); 


sigma2 = 0.0000000002; 
sigma = sqrt(sigma2); 
noisy_signal = zoh_signal + sigma * randn(size(zoh_signal)); % Add noise to signal


integrated_response = zeros(size(noisy_signal)); 
samples_per_bit = fs; 

for i = 1:N
    
    start_sample = (i-1) * samples_per_bit + 1; 
    end_sample = i * samples_per_bit; 
    integrated_response(start_sample:end_sample) = ...
        cumsum(noisy_signal(start_sample:end_sample)) * dt; % Integration with time step
    
    if i < N
        integrated_response(end_sample+1:end) = integrated_response(end_sample+1:end) - integrated_response(end_sample);
    end
end

% Plotting
figure;

% Plot Noisy Signal
subplot(3, 1, 1);
plot(t, noisy_signal, 'LineWidth', 1.5);
title('Noisy Received Signal (\sigma^2 = neligiable)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot Integrated Response
subplot(3, 1, 2);
plot(t, integrated_response, 'LineWidth', 1.5);
title('Response of Integrate & Dump Filter');
xlabel('Time (s)');
ylabel('Integrated Value');
grid on;

% Overlay Detected Bits
subplot(3, 1, 3);
time_detected = (0:N-1) * T_b + T_b/2; % Time vector for bit decisions
stem(time_detected, polar_nrz, 'r', 'LineWidth', 1.5); % Original bits
hold on;
stem(time_detected, integrated_response(1:fs:end), 'b', 'LineWidth', 1.5); % Integrated values
title('Detected Bits and Integrated Response');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Bits', 'Integrated Response');
grid on;


% Given parameters
fs = 100 * 1000; % Sampling frequency (100x bit rate of 1000 bps)
B = fs / 2; % Bandwidth (Hz)
V = 1; % Pulse amplitude
sigma2_values = [20, 40, 80, 160]; % Noise power values

% Initialize results
N0_values = sigma2_values / B; % Noise spectral density
snr_values = V^2 ./ sigma2_values; % SNR for each noise power
Pe_values = zeros(size(sigma2_values)); % Bit error probabilities

% Calculate theoretical error probabilities
for i = 1:length(sigma2_values)
    Pe_values(i) = qfunc(sqrt(snr_values(i))); % Q-function for error probability
end

% Display results
disp('Theoretical Error Probabilities:');
disp(table(sigma2_values', N0_values', snr_values', Pe_values', ...
    'VariableNames', {'NoisePower', 'N0', 'SNR', 'Pe'}));


Rb = 1e3; 
fs = 100 * Rb; 
T_sim = 10; 
N = Rb * T_sim; 
oversample_factor = fs / Rb; 
V = 1; 
sigma2_values = [20, 40, 80, 160]; 


binary_data = randi([0, 1], 1, N);


polar_nrz = (2 * binary_data - 1) * V;


zoh_signal = repelem(polar_nrz, oversample_factor);

simulated_ber = zeros(size(sigma2_values));
theoretical_ber = zeros(size(sigma2_values));

% Loop through each noise power value
for idx = 1:length(sigma2_values)
    sigma2 = sigma2_values(idx); % Noise power
    sigma = sqrt(sigma2); % Noise standard deviation
    
    % Add AWGN to the signal
    noisy_signal = zoh_signal + sigma * randn(size(zoh_signal));
    
    % Integrate & Dump (I&D) filter
    int_dump_output = zeros(1, N); % Initialize I&D output
    for bit = 1:N
        start_sample = (bit - 1) * oversample_factor + 1;
        end_sample = bit * oversample_factor;
        int_dump_output(bit) = sum(noisy_signal(start_sample:end_sample));
    end
    
    % Decision-making (Thresholding)
    detected_data = int_dump_output > 0; % Decision threshold at 0
    
    % Calculate BER
    errors = sum(binary_data ~= detected_data); % Count errors
    simulated_ber(idx) = errors / N; % Simulated BER
    
    % Theoretical BER
    snr = V^2 / sigma2; % Signal-to-noise ratio
    theoretical_ber(idx) = qfunc(sqrt(snr)); % Theoretical BER for I&D filter
end

% Display results
disp('Performance of Integrate & Dump Filter:');
disp(table(sigma2_values', simulated_ber', theoretical_ber', ...
    'VariableNames', {'NoisePower', 'SimulatedBER', 'TheoreticalBER'}));

% Plot results
figure;
semilogy(sigma2_values, simulated_ber, 'r-o', 'LineWidth', 1.5); hold on;
semilogy(sigma2_values, theoretical_ber, 'b--s', 'LineWidth', 1.5);
grid on;
xlabel('Noise Power (\sigma^2)');
ylabel('Bit Error Rate (BER)');
title('BER Comparison: Integrate & Dump Filter');
legend('Simulated BER', 'Theoretical BER', 'Location', 'Best');


