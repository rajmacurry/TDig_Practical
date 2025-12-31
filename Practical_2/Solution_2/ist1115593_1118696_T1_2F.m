HT = H';
[p, n] = size(HT);  
k = n - p;         
P = HT(:, 1:k);      

% Compute the generator matrix G
G = [eye(k), P'];

data = load('ST_Dados_G07_T1.mat');
d = data.ST_Dados_G07_T3;
subset_size = size(G, 1);

num_subsets = floor(length(d) / subset_size);
c = zeros(num_subsets, size(G, 2));

for i = 1:num_subsets
    start_index = (i - 1) * subset_size + 1;
    end_index = start_index + subset_size - 1;
    current_subset = d(start_index:end_index);
    c(i, :) = current_subset * G;
end
 
num_rows = size(c, 2);

e = zeros(num_subsets, num_rows);

for i = 1:num_subsets
    error_position = randi(num_rows); 
    e(i, error_position) = 1;
end


r = zeros(size(c)); 

for i = 1:size(c, 1)
    
    r(i, :) = mod(c(i, :) + e(i, :), 2);
end



num_syndrome_bits = size(H, 2); 
s = zeros(num_subsets, num_syndrome_bits);

for i = 1:num_subsets
    
    s(i, :) = mod(r(i, :) * H, 2); 
end



corrected_r = r; 

for i = 1:num_subsets
    if all(s(i, :) == 0)
        continue; 
    else
        for j = 1:size(H, 2)
            if isequal(s(i, :), H(:, j)')
                corrected_r(i, j) = mod(corrected_r(i, j) + 1, 2);
                break;
            end
        end
    end
end

transmitted_messages = corrected_r(:, 1:p);
syndrome = mod(c * H, 2);

