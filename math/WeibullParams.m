function [k, lambda] = WeibullParams(mu, sigma)

% Define the objective function to minimize
objective = @(k) (mu - gamma(1 + 1/k) * (sigma / sqrt(gamma(1 + 2/k) - (gamma(1 + 1/k))^2)))^2 + ...
                 (sigma^2 - (gamma(1 + 2/k) - gamma(1 + 1/k)^2) * (sigma / sqrt(gamma(1 + 2/k) - gamma(1 + 1/k)^2))^2)^2;

% Initial guess for k
k_initial_guess = 1.5;

% Solve for k using fminsearch
k = fminsearch(objective, k_initial_guess);

% Compute lambda using the solved k
lambda = mu / gamma(1 + 1/k);


end

