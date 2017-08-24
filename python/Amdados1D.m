function Amdados1D()
%------------------------------------------------------------------------------
% Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
% Copyright : IBM Research Ireland, 2017
%------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
% 1D Advection-diffusion PDE solver with data assimilation based on Kalman filter.
%
% Use the following command-line options to run this script from the project root directory:
%   matlab -nodisplay -nodesktop -nosplash -r "addpath('./python'); Amdados1D(); exit(0);"
%
% On Mac OS the full path to Matlab is, for example, "/Applications/MATLAB_R2016a.app/bin/matlab",
% so the command looks as follows:
%   /Applications/MATLAB_R2016a.app/bin/matlab -nodisplay -nodesktop -nosplash \
%           -r "addpath('./python'); Amdados1D(); exit(0);"
%--------------------------------------------------------------------------------------------------

close all;
figure(1);

N = 100;
domain_size = 1000;
T = 1500;
diffusion_coef = 7.0;
flow_velocity = 1;

noise_sigma_Q = 1;
noise_sigma_R = 1;
cov_sigma = 3;

% Space and time steps.
dx = domain_size / (N - 1);
dt = min( (dx^2)/(2*diffusion_coef), dx/flow_velocity );    % respects stability criteria

A = ModelMatrix(N, diffusion_coef, flow_velocity, dx, dt);
[field, P, H, Q, R] = InitKalmanFilter(N, noise_sigma_Q, noise_sigma_R, cov_sigma);

% Compute and accumulate observations. Also compute and record the "true" field for comparison.
[observations, true_field, max_field] = GenerateObservations(InitialField(N), N, A, H, T, dt);

% Time integration loop.
count = 0;
for t=0:dt:T
    count = count + 1;

    % Make one step forward in time with the model and apply boundary conditions.
    field = A * field;
    field(1) = 0;
    field(N) = 0;

    %S = svd(P); fprintf('lambda_min = %f,  lambda_max = %f\n', min(S(:)), max(S(:)));

    % Get observations obtained from the "true" field (simulation).
    z = observations(:,count);

    % Apply Kalman filter at few points of observations.
    P = A * P * transpose(A) + Q;
    y = z - H * field;
    K = (P * transpose(H)) * pinv(H * P * transpose(H) + R);
    P = (eye(N) - K * H) * P;
    field = field + K * y;

    % Visualizations.
    max_field = max(max_field, max(field(:)));
    if true
        plot(true_field(:,count), '-b');
        hold on;
        plot(field, '-r');
        hold off;
        ylim([0 max(max_field,1)]);
    else
        plot(sqrt(abs(field)));
        ylim([0 sqrt(abs(max_field))]);
    end
    drawnow;
    fprintf('.');
    if mod(count, 80) == 0, fprintf('\n'); end
    %pause;
end
fprintf('\n');

end

%--------------------------------------------------------------------------------------------------
% Compute and accumulate observations. Also compute and record the "true" field for comparison.
%--------------------------------------------------------------------------------------------------
function [observations, true_field_history, max_field] = GenerateObservations( ...
                                                            initial_true_field, N, A, H, T, dt)
    % Count the number of steps.
    count = 0;
    for t=0:dt:T
        count = count + 1;
    end

    % Initialize the placeholders and temporary variables.
    field = initial_true_field;
    true_field_history = zeros(N, count);
    observations = zeros(size(H,1), count);

    % Time-integration of the "true" field.
    count = 0;
    for t=0:dt:T
        count = count + 1;
        field = A * field;
        true_field_history(:,count) = field;
        observations(:,count) = H * field;
    end

    % Maximum field value for visualization.
    max_field = max(true_field_history(:));

end

%--------------------------------------------------------------------------------------------------
% Compute model matrix, which is the inverse of differentiating one: x_{t+1} = A * x_{t}.
%--------------------------------------------------------------------------------------------------
function A = ModelMatrix(N, diffusion_coef, flow_velocity, dx, dt)

    rho = (diffusion_coef * dt) / dx^2;
    vel_0 = 2*dx/dt;

    diag_elem  = 1 + 2*rho;
    sub_diag   = -flow_velocity/vel_0 - rho;
    super_diag = +flow_velocity/vel_0 - rho;

    A = zeros(N, N);
    for k=2:(N-1)
        A(k, k-1) = sub_diag;
        A(k, k)   = diag_elem;
        A(k, k+1) = super_diag;
    end
    A(1,1) = 1;
    A(N,N) = 1;
    A = pinv(A);

end

%--------------------------------------------------------------------------------------------------
% Create initial ("true") field with a spike at some point and zeros elsewhere. Note, the spike
% is not very sharp to make the field differentiable.
%--------------------------------------------------------------------------------------------------
function field = InitialField(N)

    field = zeros(N,1);
    sigma = 3;                              % in logical units (point indices)
    for x=(-4*sigma):(4*sigma)
        field(x + 15) = 10000 * exp(-x^2/(2*sigma^2)) / (sigma * sqrt(2*pi));
    end
    field(1) = 0;
    field(N) = 0;

end

%--------------------------------------------------------------------------------------------------
% Function initializes the Kalman filter, its matrices and correspoding state variables.
%--------------------------------------------------------------------------------------------------
function [field, P, H, Q, R] = InitKalmanFilter(N, noise_sigma_Q, noise_sigma_R, cov_sigma)

    % Initially all zeros because we have no idea what is the "true" field.
    field = zeros(N,1);

    % Covariance matrix. Initially it is diagonal one with quite narrow distribution
    % around the main diagonal.
    P = eye(N);
    for i=1:N
        for j=(i+1):N
            v = exp(-(i-j)^2/(2*cov_sigma^2));
            %v = 0;
            P(i,j) = v;
            P(j,i) = v;
        end
    end

    % Observation matrix H has units on a sub-set of observations and zeros elsewhere.
    % The step of observation grid corresponds to the initial correlation length "cov_sigma".
    assert(N > 10);
    num_obs = min(max(round(N/cov_sigma), 1), N);   % number of observations
    fprintf('num. observations = %d\n', num_obs);
    idx = round(linspace(1, N, num_obs));
    assert(length(idx) == num_obs);
    H = zeros(num_obs, N);
    for k=1:num_obs
        H(k,idx(k)) = 1;
    end

    % Noise covariance matrices of the model and observation respectively.
    Q = eye(N)         * noise_sigma_Q;
    R = eye(size(H,1)) * noise_sigma_R;

end


