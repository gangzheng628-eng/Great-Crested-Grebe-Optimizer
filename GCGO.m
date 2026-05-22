function [global_best, best_fitness, best_fitness_history] = GCGO(fobj, dim, lb, ub, pop_size, max_iter)
%GCGO Great Crested Grebe Optimizer.
%
%   This function implements the Great Crested Grebe Optimizer (GCGO), a
%   seasonally switched nature-inspired metaheuristic for continuous
%   optimization. The algorithm models two plumage phases of the great
%   crested grebe and uses four behavior-inspired mechanisms to coordinate
%   global exploration and local exploitation.
%
%   Main mechanisms:
%       1. Molting-based phase switching
%          The stochastic molting factor controls the transition between the
%          winter plumage phase and the summer plumage phase.
%
%       2. Population density cycle
%          The active population size is dynamically adjusted during the
%          iterative search process.
%
%       3. Winter plumage phase
%          Diving foraging and defense escape are used to enhance global
%          exploration and population diversity.
%
%       4. Summer plumage phase
%          Courtship and reproduction and brood carrying are used to enhance
%          coordinated exploitation and local refinement.

%% Parameter settings
if isscalar(lb)
    lb = repmat(lb, 1, dim);
else
    lb = lb(:)';
end

if isscalar(ub)
    ub = repmat(ub, 1, dim);
else
    ub = ub(:)';
end

N0 = pop_size;
eta_B = 0.03;
alpha_N = 0.05;
molting_threshold = 1.9290;

B = eta_B * N0;

N_func = @(t) round(N0 - B * sin(2 * pi * t / max_iter) * log(1 + alpha_N * t));

N_schedule = zeros(max_iter, 1);
for t = 1:max_iter
    N_schedule(t) = N_func(t);
end

N_schedule = max(N_schedule, 2);
pop_max = max(N_schedule);

%% Population initialization
pop = repmat(struct('position', [], 'fitness', []), pop_max, 1);

for i = 1:pop_max
    pop(i).position = lb + rand(1, dim) .* (ub - lb);
    pop(i).fitness = fobj(pop(i).position);
end

[best_fitness, best_idx] = min([pop.fitness]);
global_best = pop(best_idx).position(:)';

best_fitness_history = zeros(max_iter, 1);

%% Main loop
for t = 1:max_iter

    current_N = N_schedule(t);

    [best_fitness, best_idx] = min([pop(1:current_N).fitness]);
    global_best = pop(best_idx).position(:)';

    r_molting = max(rand, realmin);
    M = 2 * cos((1 - t / max_iter) * log(1 / r_molting));

    active_fitness = [pop(1:current_N).fitness];
    total_fit = sum(abs(active_fitness - min(active_fitness))) + eps;

    for i = 1:current_N

        xi = pop(i).position;
        fi = pop(i).fitness;
        gi = global_best;

        if M < molting_threshold

            %% Mechanism 1: Diving foraging
            if rand < rand
                xr = pop(randi(current_N)).position;
                qi = 0.5 * (xi + xr);

                water_disturbance = (0.3 + 0.4 * rand) .* randn(1, dim);
                x_new = xi + water_disturbance .* abs(2 * rand(1, dim) .* xi - qi);

            %% Mechanism 2: Defense escape
            else
                x1 = pop(randi(current_N)).position;
                x2 = pop(randi(current_N)).position;
                predator = pop(randi(current_N)).position;

                U = randi([0, 1], 1, dim);

                diff_r = xi - predator;
                dist_sq = sum(diff_r .^ 2) + eps;

                F_coulomb = randn(1, dim) .* diff_r ./ dist_sq;

                sigma = 0.5 * norm(ub - lb) + eps;
                escape_coeff = exp(-dist_sq / (2 * sigma ^ 2));

                explore_step = rand(1, dim) .* (x2 - x1);
                x_candidate = xi + explore_step + escape_coeff .* F_coulomb;

                x_new = (1 - U) .* xi + U .* x_candidate;
            end

        else

            %% Mechanism 3: Courtship and reproduction
            if rand > rand
                p1 = pop(randi(current_N)).position;
                p2 = pop(randi(current_N)).position;

                x3 = pop(randi(current_N)).position;
                x4 = pop(randi(current_N)).position;
                x5 = pop(randi(current_N)).position;

                St = 1 - exp(-abs(fi - min(active_fitness)) / total_fit);
                theta_c = 0.8 - 0.3 * cos(pi * t / max_iter);
                theta_a = theta_c * (0.2 + 0.6 * St);

                d = x3 - x4;
                v = x5 - x4;

                v = v - (dot(v, d) / (dot(d, d) + eps)) * d;
                v = v / (norm(v) + eps);

                x_new = xi ...
                    + theta_a .* (p1 + p2 - 2 * xi) ...
                    + theta_c .* (d + rand(1, dim) .* v);

            %% Mechanism 4: Brood carrying
            else
                parent = pop(randi(current_N)).position;
                chick = pop(randi(current_N)).position;

                x6 = pop(randi(current_N)).position;
                x7 = pop(randi(current_N)).position;

                fitness_gap = abs(fi - best_fitness);
                perception = fitness_gap / total_fit;
                rho_i = max((1 + perception * dim) ^ (-0.3), 0.2);

                levy_step = levy_flight(dim);
                beta_t = 0.5 + 0.5 * (1 - t / max_iter);

                x_new = xi ...
                    + rho_i .* rand(1, dim) .* (parent - chick) ...
                    + rho_i .* rand(1, dim) .* (x6 - x7) ...
                    + beta_t .* levy_step .* (gi - xi);
            end
        end

        %% Boundary handling
        x_new = min(max(x_new, lb), ub);

        %% Fitness evaluation
        f_new = fobj(x_new);

        %% Greedy replacement
        if f_new < fi
            pop(i).position = x_new;
            pop(i).fitness = f_new;
        end
    end

    %% Current best update
    [best_fitness, best_idx] = min([pop(1:current_N).fitness]);
    global_best = pop(best_idx).position(:)';

    best_fitness_history(t) = best_fitness;
end

end

%% Levy flight
function L = levy_flight(dim)

beta = 1.5;

sigma_u = (gamma(1 + beta) * sin(pi * beta / 2) / ...
    (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);

u = randn(1, dim) * sigma_u;
v = randn(1, dim);

L = u ./ (abs(v) .^ (1 / beta) + eps);

end