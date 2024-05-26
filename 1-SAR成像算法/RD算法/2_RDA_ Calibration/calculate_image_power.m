function epsilon_p = calculate_image_power(DN, delta_a, delta_b)
    % Check the size of the input DN matrix
    [rows, cols] = size(DN);
    
    % Ensure the matrix is square
    if rows ~= cols
        error('Input DN matrix must be square.');
    end
    
    % Define the original size
    original_size = 15;
    
    % Calculate the scaling factor
    scale_factor = rows / original_size;
    
    if mod(scale_factor, 1) ~= 0
        error('Input DN matrix size must be an integer multiple of %d.', original_size);
    end
    
    % Calculate the dimensions of N_A and N_B based on the scaling factor
    NA_size = 3 * scale_factor; % N_B is a 3x3 region at the center
    NB_size = 4 * scale_factor;   % N_A is 4x4 regions in each corner

    % Identify the indices for N_B regions (4x4 in each corner)
    NB_Top_Left = DN(1:NB_size, 1:NB_size); % Top-left
    NB_Top_Right = DN(1:NB_size, end-NB_size+1:end); % Top-right
    NB_Bottom_Left = DN(end-NB_size+1:end, 1:NB_size); % Bottom-left
    NB_Bottom_right = DN(end-NB_size+1:end, end-NB_size+1:end); % Bottom-right

    % Calculate the sum of squares for N_B regions
    NB_pixel = numel(NB_Top_Left) + numel(NB_Top_Right) + numel(NB_Bottom_Left) + numel(NB_Bottom_right); % NB total pixel
    NB_sum = sum(sum(NB_Top_Left.^2)) + sum(sum(NB_Top_Right.^2)) + sum(sum(NB_Bottom_Left.^2)) + sum(sum(NB_Bottom_right.^2)); % NB total DN

    % Identify the center of the matrix
    center = (rows + 1) / 2;
    
    % Calculate the indices for N_B region (3 rows and 3 columns at the center)
    NA_rows = DN(center-1: center+1, 1: end);
    NA_cols = DN(1:end, center-1:center+1);
    NA_center = DN(center-1:center+1, center-1:center+1);

    % Calculate the sum of squares for N_A regions
    NA_pixel = numel(NA_rows) + numel(NA_cols) - numel(NA_center); % NA total pixel
    NA_sum = sum(sum(NA_rows.^2)) + sum(sum(NA_cols.^2)) - sum(sum(NA_center.^2));
    
    % Calculate epsilon_p using the given formula
    epsilon_p = (NA_sum - (NA_pixel / NB_pixel) * NB_sum) * delta_a * delta_b;
    
end
