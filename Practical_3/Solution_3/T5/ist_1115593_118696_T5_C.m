% Parameters
R = 1000;           % Bit rate (bits per second)
fs = 100 * R;       % Sampling frequency (Hz)
Spb = fs / R;       % Samples per bit
N_bits = 1000;      % Number of bits in the binary sequence
B = 4000;           % Channel bandwidth (Hz)
G = 1;              % Channel gain
noise_powers = [0.01, 1]; % Noise powers to test

% Generate random binary sequence
binary_sequence = randi([0, 1], 1, N_bits);

% Unipolar NRZ encoding
x = repelem(binary_sequence, Spb); % Repeat each bit Spb times

% Design the channel (Low-Pass Filter)
fc = B;                   % Cutoff frequency (Hz)
[b, a] = butter(6, fc / (fs / 2)); % 6th order Butterworth filter

% Compute PSD of the input signal x
[psd_x, f] = periodogram(x, [], [], fs, 'power');

% Initialize figure for PSD comparison
figure;
hold on;

% Loop over noise powers
for i = 1:length(noise_powers)
    N = noise_powers(i);

    % Pass signal through the channel
    y = filter(b, a, x);                % Filtered signal
    noise = sqrt(N) * randn(size(x));  % Gaussian noise
    y_noisy = G * y + noise;           % Channel output with noise

    % Compute PSD of the channel output
    [psd_y, ~] = periodogram(y_noisy, [], [], fs, 'power');

    % Plot PSD of the channel output
    plot(f, 10 * log10(psd_y), 'DisplayName', ['PSD (N = ', num2str(N), ')']);
end

% Plot PSD of the input signal
plot(f, 10 * log10(psd_x), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Input Signal (x)');

% Graph Labels
title('Power Spectral Density of Input and Output Signals');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('show');
grid on;
hold off;
