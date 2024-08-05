function [rootsArray] = roots_cols(A)

% Define the matrix (each column represents polynomial coefficients)

% Get the number of columns
numCols = size(A, 2);
numRows = size(A,1);
% Initialize a cell array to store the roots of each column
rootsArray = zeros(numRows-1, numCols);

% Loop through each column and find the roots
for i = 1:numCols
    % Extract the column vector
    polyCoeffs = A(:, i);
    
    % Find the roots of the polynomial
    rootsArray(:,i) = roots(polyCoeffs);
end

