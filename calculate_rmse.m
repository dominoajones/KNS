function rmse_value = calculate_rmse(matrix1, matrix2) %matrix1 should be modelled, matrix2 should be observed
    % Check if the input matrices have the same dimensions
    if ~isequal(size(matrix1), size(matrix2))
        error('Input matrices must have the same dimensions');
    end

    % Calculate the differences between the matrices
    differences = matrix1 - matrix2;

    % Square the differences
    squared_differences = differences .^ 2;

    % Compute the mean of the squared differences
    mean_squared_error = mean(squared_differences(:));

    % Take the square root of the mean squared error to get RMSE
    rmse_value = sqrt(mean_squared_error);
end