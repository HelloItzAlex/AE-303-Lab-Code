function [deltaPs] = dP2(data)

    array_data = table2array(data);

    % Get the number of test runs
    num_tests = size(array_data, 1) / 4;
 
    % Preallocate the output matrix
    deltaPs = zeros(num_tests, 1);
    
    % Loop through each test and compute the difference
    for i = 1:num_tests
        row1 = (i - 1) * 4 + 3
        row2 = row1 + 1
        
        % Compute the average of each row before computing the difference
        avg_row1 = mean(array_data(row1, :));
        avg_row2 = mean(array_data(row2, :));
        
        deltaPs(i, 1) = avg_row2 - avg_row1;
    end
end