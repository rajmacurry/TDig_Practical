p = 8;

H = hammgen(p);
H = H';

min_col_sum_zero = inf; 
index_min_col_sum_zero  = [];
[num_rows, ~] = size(H);

for r = 1:num_rows
    combinations = nchoosek(1:num_rows, r);
    
    for i = 1:size(combinations, 1)
        row_sum = sum(H(combinations(i, :), :), 1);
        
        if all(mod(row_sum, 2) == 0) 
            if r < min_col_sum_zero
                min_col_sum_zero = r; 
                index_min_col_sum_zero  = combinations(i, :); 
            end
            break; 
        end
    end
    
    if ~isempty(index_min_col_sum_zero )
        break;  
    end
end


% Q8-9
syndrome = load('Sindroma_G07_T1_F1.mat');

col_sum2_sind = [];

for r = 1:2
    combinations = nchoosek(1:num_rows, r);
    for i = 1:size(combinations, 1)
        row_sum = mod(sum(H(combinations(i, :), :), 1), 2);
        if isequal(row_sum, syndrome.sindroma')
            if r==1
                col_sum2_sind = [col_sum2_sind; combinations(i, 1), 0]; 
            else
                col_sum2_sind = [col_sum2_sind; combinations(i, :)];            end
        end
    end
end
