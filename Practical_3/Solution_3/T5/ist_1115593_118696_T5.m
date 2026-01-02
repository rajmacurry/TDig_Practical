R = 1000;      
fs = 100 * R;   
Spb = fs / R;    
N = 1000;
B = 4000;


binary_sequence = randi([0, 1], 1, N);


unipolar_nrz = repelem(binary_sequence, Spb); % repeat element


% Manchester encoding
manchester = zeros(1, N * Spb); % Initialize Manchester signal
for i = 1:N
    % Define the time range for each bit
    start_idx = (i - 1) * Spb + 1;
    mid_idx = start_idx + Spb/2 - 1;
    end_idx = start_idx + Spb - 1;
    
    if binary_sequence(i) == 1
        manchester(start_idx:mid_idx) = 1;   % First half: High
        manchester(mid_idx+1:end_idx) = 0;  % Second half: Low
    else
        manchester(start_idx:mid_idx) = 0;   % First half: Low
        manchester(mid_idx+1:end_idx) = 1;  % Second half: High
    end
end

% Design the channel (Low-Pass Filter)
fc = B;                   % Cutoff frequency (Hz)
[b, a] = butter(6, fc/(fs/2)); % 6th order Butterworth filter

% Pass signals through the channel
unipolar_nrz_out = filter(b, a, unipolar_nrz);  % Output of Unipolar NRZ
manchester_out = filter(b, a, manchester);      % Output of Manchester

% Noise powers to test
noise_powers = [0, 0.01, 0.1, 1];

% Time vector for plotting
t = (0:length(unipolar_nrz)-1) / fs;

% Loop through different noise powers
for i = 1:length(noise_powers)
    noise_power = noise_powers(i);

    % Generate Gaussian noise
    noise = sqrt(noise_power) * randn(size(unipolar_nrz));

    % Add noise to the channel outputs
    unipolar_nrz_noisy = unipolar_nrz_out + noise;
    manchester_noisy = manchester_out + noise;

    % Plot results for Unipolar NRZ
    figure;
    subplot(2, 1, 1);
    plot(t, unipolar_nrz_noisy, 'LineWidth', 1.5);
    title(['Unipolar NRZ - Channel Output with Noise Power N = ', num2str(noise_power)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-1, 2]);
    grid on;

    % Plot results for Manchester
    subplot(2, 1, 2);
    plot(t, manchester_noisy, 'LineWidth', 1.5);
    title(['Manchester - Channel Output with Noise Power N = ', num2str(noise_power)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-1, 2]);
    grid on;
end