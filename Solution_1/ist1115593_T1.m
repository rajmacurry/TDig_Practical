% read audio file
[fluxo_simbolos, fs] = audioread("ST_G7_T2.wav");

%unique identifies distinct symbols, ~ ignore result, index
[simb, ~, idx] = unique(fluxo_simbolos);

%compute number of bits
n_bits = ceil(log2(max(fluxo_simbolos) + 1));

%count occurance of each symbol
counts = accumarray(idx, 1);

prob_simb = counts / length(fluxo_simbolos);

% finding symbol with max probability of occurence
[prob_simb_max, idx_max] = max(prob_simb);
simb_prob_max = simb(idx_max);
max(prob_simb)


% finding symbol with min probablity of ocurence
[prob_simb_min, idx_min] = min(prob_simb);
simb_prob_min = simb(idx_min);

min(prob_simb)

% plot for the histogram of symbols
figure;
histogram(fluxo_simbolos, 'Normalization', 'probability', 'BinMethod', 'auto');
title('Histogram of Symbols');
xlabel('Symbol Value');
ylabel('Probability');
grid on;


%plot for the pdf of the symbols
figure;
bar(simb, prob_simb, 'FaceColor', 'b', 'EdgeColor', 'k');
title('Probability Density Function of Symbols');
xlabel('Symbol Value');
ylabel('Probability');
grid on;



%removing symbols with probabilities less than min probablity sym
non_zero_idx = prob_simb > 3.6550e-05;

%symbols and probabilities above threshold
simb_no_zero = simb(non_zero_idx);
prob_simb_no_zero = prob_simb(non_zero_idx);

% sorting probabilities and storing indexes
[prob_simb_no_zero, sort_idx] = sort(prob_simb_no_zero);

%sorting symbols on the basis of indexes above
simb_no_zero = simb_no_zero(sort_idx);

%no of relevant symbols
num_sym = length(prob_simb_no_zero);

%Decision content: theoretical max info content
Conteudo_Decisao = log2(num_sym);

%Entropy
Entropia = -sum(prob_simb_no_zero .* log2(prob_simb_no_zero));

%redundancy = diff betn theoritical max and actual entropy
Redundancia = Conteudo_Decisao - Entropia;
Redundancia

%rate of information
R = fs * Entropia;

%theoritical max rate
Rd = fs * Conteudo_Decisao;

%difference (redundant bits per sec)
Db_Redund = Rd - R;

R
Rd
Db_Redund

% Huffman coding

%create probability table for the huffman coding process
prob_table = table(simb_no_zero', prob_simb_no_zero', 'VariableNames', {'Symbol', 'Probability'});

%normalize
prob_table.Probability = prob_table.Probability / sum(prob_table.Probability);


%generate huffman codes
[huffman_table, compression] = huffman(prob_table.Probability);


%initialization for tabela_final table
n_symbols = length(simb_no_zero);
tabela_final = cell(n_symbols, 4);


%populate the table
for i = 1:n_symbols
    tabela_final{i, 1} = prob_table.Symbol(i); 
    tabela_final{i, 2} = prob_table.Probability(i); 
    
   
    binary_code = huffman_table(i, :); 
    binary_code = strtrim(binary_code);
    
    binary_code = char(binary_code);  % Convert cell to a character array
    non_empty_count = sum(binary_code ~= ' ');  % Now compare as a string

    tabela_final{i, 3} = non_empty_count; 
    tabela_final{i, 4} = binary_code; 
end



%convert tablela_final to a MATLAB table
tabela_final = cell2table(tabela_final, 'VariableNames', {'Symbol', 'Probability', 'CodeLength', 'BinaryCode'});

probability = prob_table.Probability;
code_length = tabela_final.CodeLength;


probability = probability(:);
code_length = code_length(:);

% average code length in bits
L_med = sum(probability .* code_length);
L_med

% entropy
H = -sum(prob_table.Probability .* log2(prob_table.Probability));



% efficiency
Eficiencia = H/L_med;
Eficiencia



% now we encode and decode symbols using the huffman codes

%initialize an empty string to store the bin representation of encoded
%symbol
fluxo_bin = ''; 

% converting fluxo_simbolos to cell array of strings, because
% containers.Map requires it
fluxo_simbolos_str = cellstr(num2str(fluxo_simbolos)); 


%creating dictionary for symbol and huffman codes to use later
symbol_to_code = containers.Map(cellstr(num2str(tabela_final.Symbol)), tabela_final.BinaryCode);

%encoding
for i = 1:length(fluxo_simbolos_str)
    symbol = fluxo_simbolos_str{i}; 
    if isKey(symbol_to_code, symbol) %checks if exists in dict
        fluxo_bin = [fluxo_bin, symbol_to_code(symbol)]; %append
    end
end

%total length of encoded bin sequence
n_fluxo_bin = length(fluxo_bin);
n_fluxo_bin

% no of orig symbols in sequence
N = length(fluxo_simbolos); 

%average code rate: avg bits per symbol in the encoded sequence
R_cod_med = n_fluxo_bin / N;
R_cod_med


%create a dict for decoding
code_to_symbol = containers.Map(tabela_final.BinaryCode, cellstr(num2str(tabela_final.Symbol)));


fluxo_simbolo_descod = []; 
current_code = ''; 


for i = 1:length(fluxo_bin)
    current_code = [current_code, fluxo_bin(i)];  %get code 1 bit by bit
    
    
    if isKey(code_to_symbol, current_code) % check for evey new bit added to current_code
        
        decoded_symbol = code_to_symbol(current_code);
        fluxo_simbolo_descod = [fluxo_simbolo_descod; str2double(decoded_symbol)]; 
        current_code = '';  %empty it after found in dict
    end
end



%now plot the original and decoded symbols up to 1 sec
duration = 1;
t_original = (0:length(fluxo_simbolos)-1) / fs; 
t_decoded = (0:length(fluxo_simbolo_descod)-1) / fs; 

t_original = t_original(t_original <= duration);
t_decoded = t_decoded(t_decoded <= duration);

original_signal_plot = fluxo_simbolos(1:length(t_original));
decoded_signal_plot = fluxo_simbolo_descod(1:length(t_decoded));

figure;
hold on; 
plot(t_original, original_signal_plot, 'b', 'DisplayName', 'Original Signal');
plot(t_decoded, decoded_signal_plot, 'r', 'DisplayName', 'Decoded Signal');
hold off;

xlabel('Time (seconds)');
ylabel('Amplitude');
title('Original Signal and Decoded Signal');
legend show;
grid on;