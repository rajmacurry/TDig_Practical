V = 1;
T_b = 1;

oversample_factor = 100;
time_vector = linspace(0, T_b, oversample_factor);

data_sequence = [0, 1, 0, 0, 1, 1];
polar_nrz = V * (2 * data_sequence - 1);


transmitted_signal = [];
for i = 1:length(data_sequence)
    transmitted_signal = [transmitted_signal, polar_nrz(i) * ones(1, oversample_factor)];
end


matched_filter_response = flip(polar_nrz);
time_vector_mf = linspace(0, length(matched_filter_response) * T_b, length(matched_filter_response));
filtered_signal = conv(transmitted_signal, matched_filter_response, 'same');
t_filtered = linspace(0, length(transmitted_signal) * T_b, length(filtered_signal));


figure;

% Plot the transmitted signal
subplot(2,1,1);
plot(linspace(0, length(transmitted_signal) * T_b, length(transmitted_signal)), transmitted_signal, 'LineWidth', 1.5);
title('Transmitted Signal (Polar NRZ)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot the impulse response of the matched filter
subplot(2,1,2);
plot(time_vector_mf, matched_filter_response, 'LineWidth', 1.5);
title('Matched Filter Impulse Response');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;


% Given values of noise power sigma^2
sigma_squared_values = [20, 40, 80, 160]; % in W

% Amplitude of the transmitted signal
V = 1;

% Preallocate an array to store theoretical error probabilities
theoretical_error_probabilities = zeros(1, length(sigma_squared_values));

% Calculate the theoretical BER for each noise power
for i = 1:length(sigma_squared_values)
    sigma2 = sigma_squared_values(i); % Noise power
    sigma = sqrt(sigma2); % Standard deviation of the noise
    
    % Calculate the theoretical error probability using the Q-function
    theoretical_error_probabilities(i) = qfunc(V / sigma); % Q-function
end

% Display the theoretical error probabilities
disp('Theoretical Error Probabilities:');
disp(theoretical_error_probabilities);