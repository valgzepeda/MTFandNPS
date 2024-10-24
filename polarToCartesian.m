function output_image = polarToCartesian(input_matrix, P)
    % input_matrix: Matrix of size NxM
    % P: Desired size of output square image (PxP)

    % Get the size of the input matrix
    [N, M] = size(input_matrix);
    
    % Define the angular positions (from 1 to 180 degrees with N steps)
    theta = linspace(0, 2*pi, N);  % From 0 to pi radians (1 to 180 degrees)
    
    % Define the radial distances (assuming they are pixel units)
    r = linspace(1, M, M);

    % Create a grid for the desired output image
    x = linspace(-M, M, P);
    y = linspace(-M, M, P);
    [X, Y] = meshgrid(x, y);
    
    % Convert Cartesian coordinates to polar coordinates
    R = sqrt(X.^2 + Y.^2);  % Radial distances
    Theta = atan2(Y, X);    % Angles in radians
    Theta(Theta < 0) = Theta(Theta < 0) + 2*pi;  % Convert range from -pi to pi to 0 to 2*pi
    Theta(Theta > pi) = Theta(Theta > pi) - pi;  % Map back to 0 to pi
    
    % Interpolate the input matrix to the output image
    % Convert polar coordinates to indices
    r_idx = interp1(r, 1:M, R, 'linear', NaN);
    theta_idx = interp1(theta, 1:N, Theta, 'linear', NaN);
    
    % Ensure indices are within valid range
    r_idx = max(min(r_idx, M), 1);
    theta_idx = max(min(theta_idx, N), 1);
    
    % Initialize the output image
    output_image = NaN(size(X));
    
    % Fill the output image using the interpolated values from the input matrix
    for i = 1:P
        for j = 1:P
            if ~isnan(r_idx(i,j)) && ~isnan(theta_idx(i,j))
                output_image(i,j) = interp2(input_matrix, r_idx(i,j), theta_idx(i,j), 'linear');
            end
        end
    end
    
    % Replace NaN values with zeros or a default background value
    output_image(isnan(output_image)) = 0;
end
