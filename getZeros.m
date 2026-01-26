function zeros = getZeros(A, B, C, D)
    syms s; % Define symbolic variable
    n = size(A,1);
    M = [s*eye(n)-A, -B; C, D];
    zero_eq = det(M); % Compute determinant symbolically
    zeros = vpasolve(zero_eq, s); % Use vpasolve for numerical solutions
    zeros = double(zeros); % Convert symbolic results to double
end