function tmpAmdados2D()
%------------------------------------------------------------------------------
% Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
% Copyright : IBM Research Ireland, 2017
%------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
% Advection-diffusion PDE solver with data assimilation based on Kalman filter.
%
% Use the following command-line options to run this script from the project root directory:
%   matlab -nodisplay -nodesktop -nosplash -r "addpath('./python'); Amdados2D(); exit(0);"
%
% On Mac OS the full path to Matlab is, for example, "/Applications/MATLAB_R2016a.app/bin/matlab",
% so the command looks as follows:
%   /Applications/MATLAB_R2016a.app/bin/matlab -nodisplay -nodesktop -nosplash \
%           -r "addpath('./python'); Amdados2D(); exit(0);"
%--------------------------------------------------------------------------------------------------

close all;

% Primary parameters set by user.
conf.output_dir = 'output';           % the output folder for the main application
conf.diffusion_coef = 7.0;            % coefficient of diffusion

conf.nx = 60;                         % number of nodal points in x dimension
conf.ny = 60;                         % number of nodal points in y dimension

conf.domain_size_x = 1000;            % global domain size in x-dimension [meters]
conf.domain_size_y = 1000;            % global domain size in y-dimension [meters]

conf.spot_x = 423;                     % substance spot abscissa in the global domain [meters]
conf.spot_y = 587;                    % substance spot ordinate in the global domain [meters]

conf.spot_density = 10000;            % substance spot concentration at initial time [units???]
conf.observation_nx = 30;             % #observation points in x dimension
conf.observation_ny = 30;             % #observation points in y dimension

conf.integration_period = 3500;       % integration period 0...T [seconds]
conf.integration_nsteps = 200;        % min. number of time steps to be done

conf.flow_model_max_vx = 1.0;         % module of max. flow velocity in x-dimension [meters/seconds]
conf.flow_model_max_vy = 1.0;         % module of max. flow velocity in y-dimension [meters/seconds]

conf.model_ini_var          = 1.0;    % initial variance (diag. elements of model cov. matrix P)
conf.model_ini_covar_radius = 1.0;    % initial radius of correlation among domain points [meters]
conf.model_noise_Q          = 1.0;    % amplitude (here deviation) of model noise
conf.model_noise_R          = 1.0;    % amplitude (here deviation) of measurement noise

global ZERO_BOUNDARY_CONDITIONS;
ZERO_BOUNDARY_CONDITIONS = false;

% Initialize the rest of parameters that can be deduced from the primary ones.
conf = InitDependentParams(conf);
glo_idx = GlobalIndices(conf);
disp(conf);

% Run forward simulation and record the "true" solutions.
true_fields = ComputeTrueFields(conf, glo_idx);

% ***** Simulation with Kalman filter. *****

% Initially the state field is all zeros because we have no idea what is the "true" one.
field = zeros(conf.nx, conf.ny);
P = InitialCovar(conf, glo_idx);
H = ComputeH(conf, glo_idx);

% Open AVI file for writing.
[status, cmdout] = system(['mkdir -p ' conf.output_dir]);
assert(status == 0, cmdout);
outputVideo = VideoWriter([conf.output_dir filesep 'both_fields.avi']);
outputVideo.FrameRate = 15;
open(outputVideo);
figure(1);

try
    % Time integration with data assimilation.
    for k=1:conf.Nt
        fprintf('.');
        if (mod(k,80) == 0), fprintf('\n'); end
        big_image = horzcat(transpose(true_fields(:,:,k)), ...
                                    zeros(size(field,2),10), transpose(field));
        vmin = min(big_image(:));
        vmax = max(big_image(:));
        big_image = abs((big_image - vmin) / max(vmax - vmin, eps));
        imagesc(big_image);
        title('True density (left) vs. Kalman filtering (right)');
        drawnow();
        writeVideo(outputVideo, big_image);

        Q = ComputeQ(conf);
        R = ComputeR(conf);

        % Compute inverse A (i.e. the inverse model matrix).
        t = (k - 1) * conf.dt;
        invA = InverseModelMatrix(conf, glo_idx, t);

        % x_prior <- A * x_previous.
        field1D = reshape(field, conf.problem_size, 1);
        field = invA \ field1D;
        field = reshape(field, conf.nx, conf.ny);
        field = ApplyBoundaryCondition(conf, field);

        % P_prior <- A * P_previous * A^t + Q.
        P = invA \ transpose(invA \ P) + Q;
        P = (P + transpose(P)) * 0.5;           % symmetrize

        % Residual between the "true" observations and the currently simulated ones.
        z = H * reshape(true_fields(:,:,k), conf.problem_size, 1);
        y = z - H * field1D;

        % S = H * P_prior * H^t + R;
        S = H * P * transpose(H) + R;
        S = (S + transpose(S)) * 0.5;           % symmetrize

        % x = x_prior + K * y = x_prior * P * H^t * S^{-1} * y.
        field1D = reshape(field, conf.problem_size, 1);
        field1D = field1D + P * transpose(H) * (S \ y);
        field = reshape(field1D, conf.nx, conf.ny);

        % P_new = (I - K*H) * P_prior = P_prior - P_prior * H^t * S^{-1} * H * P_prior.
        P = P - P * transpose(H) * (S \ (H * P));
    end
    fprintf('\n\n\n');
catch ME
    close(outputVideo);
    disp(ME);
    for q=1:length(ME.stack), disp(ME.stack(q)); end
    fprintf('\n\n\n');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------
% Function initializes dependent parameters given the primary ones specified by user.
%--------------------------------------------------------------------------------------------------
function conf = InitDependentParams(conf)

    D = conf.diffusion_coef;
    assert(D > 0);

    conf.problem_size = conf.nx * conf.ny;
    dx = conf.domain_size_x / (conf.nx-1);
    dy = conf.domain_size_y / (conf.ny-1);
    assert((dx > 0) && (dy > 0));
    conf.dx = dx;
    conf.dy = dy;

    % Deduce the optimal time step from the stability criteria.
    tiny = realmin / eps^3;
    dt_base = conf.integration_period / conf.integration_nsteps;
    max_vx = conf.flow_model_max_vx;
    max_vy = conf.flow_model_max_vy;
    dt = min(dt_base, min( min(dx^2, dy^2)/(2.0*D + tiny), ...
                           1.0/(abs(max_vx)/dx + abs(max_vy)/dy + tiny) ));
    assert(dt > 0);
    conf.dt = dt;
    conf.Nt = ceil(conf.integration_period / dt);

    conf.rho_x = D * dt / dx^2;
    conf.rho_y = D * dt / dy^2;

    conf.v0x = 2.0 * dx / dt;
    conf.v0y = 2.0 * dy / dt;

end

%--------------------------------------------------------------------------------------------------
% Each nodal point gets a unique global index on the grid. Function initializes a 2D array
% of indices:  index of (x,y) = glo_idx(x,y).
%--------------------------------------------------------------------------------------------------
function glo_idx = GlobalIndices(conf)

    glo_idx = reshape((1:conf.problem_size), conf.nx, conf.ny);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advection-diffusion PDE stuff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------
% Function applies zero Dirichlet boundary condition to the state field.
%--------------------------------------------------------------------------------------------------
function field = ApplyBoundaryCondition(conf, field)

    global ZERO_BOUNDARY_CONDITIONS;
	assert((size(field,1) == conf.nx) && (size(field,2) == conf.ny));
    if ZERO_BOUNDARY_CONDITIONS
        field(1,:) = 0;  field(end,:) = 0;
        field(:,1) = 0;  field(:,end) = 0;
    end

end

%--------------------------------------------------------------------------------------------------
% Function initializes inverse matrix of implicit Euler time-integrator:
% B * x_{t+1} = x_{t}, where B = A^{-1} is the matrix returned by this function.
% The matrix must be inverted while iterating forward in time: x_{t+1} = A * x_{t}.
%--------------------------------------------------------------------------------------------------
function invA = InverseModelMatrix(conf, glo_idx, t)

    global ZERO_BOUNDARY_CONDITIONS;
    nx = conf.nx;
    ny = conf.ny;

    if ZERO_BOUNDARY_CONDITIONS
        num_nz = 4 + 2*(nx-2) + 2*(ny-2) + 5*(nx-2)*(ny-2);
    else
        num_nz = 5*nx*ny;
    end
    rows = zeros(num_nz, 1);
    cols = zeros(num_nz, 1);
    vals = zeros(num_nz, 1);

    [vx, vy] = Flow(conf, t);
    rho_x = conf.rho_x;
    rho_y = conf.rho_y;
    vx = vx / conf.v0x;
    vy = vy / conf.v0y;

    N = 0;
    for y=1:ny
    for x=1:nx
        i = glo_idx(x,y);
        if (x == 1) || (x == nx) || (y == 1) || (y == ny)
            if ZERO_BOUNDARY_CONDITIONS
                N = N + 1;  rows(N) = i;  cols(N) = i;  vals(N) = 1;
            else
                N = N + 1;  rows(N) = i;  cols(N) = i;  vals(N) = 1 + 2*(rho_x + rho_y);

                if x == 1
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,  y);  vals(N) = - 2*vx - rho_x;
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x+1,y);  vals(N) = + 2*vx - rho_x;
                elseif x == nx
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x-1,y);  vals(N) = - 2*vx - rho_x;
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,  y);  vals(N) = + 2*vx - rho_x;
                else
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x-1,y);  vals(N) = - vx - rho_x;
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x+1,y);  vals(N) = + vx - rho_x;
                end

                if y == 1
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y  );  vals(N) = - 2*vy - rho_y;
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y+1);  vals(N) = + 2*vy - rho_y;
                elseif y == ny
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y-1);  vals(N) = - 2*vy - rho_y;
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y  );  vals(N) = + 2*vy - rho_y;
                else
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y-1);  vals(N) = - vy - rho_y;
                    N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y+1);  vals(N) = + vy - rho_y;
                end
            end
        else
            N = N + 1;  rows(N) = i;  cols(N) = i;               vals(N) = 1 + 2*(rho_x + rho_y);
            N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x-1,y);  vals(N) = - vx - rho_x;
            N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x+1,y);  vals(N) = + vx - rho_x;
            N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y-1);  vals(N) = - vy - rho_y;
            N = N + 1;  rows(N) = i;  cols(N) = glo_idx(x,y+1);  vals(N) = + vy - rho_y;
        end
    end
    end
    assert(N == num_nz);
    invA = sparse(rows, cols, vals, conf.problem_size, conf.problem_size, num_nz);

end

%--------------------------------------------------------------------------------------------------
% Function computes flow components given a time.
%--------------------------------------------------------------------------------------------------
function [vx, vy] = Flow(conf, t)

    vx = -conf.flow_model_max_vx * sin(0.1 * t / conf.integration_period - pi);
    vy = -conf.flow_model_max_vy * sin(0.2 * t / conf.integration_period - pi);

end

%--------------------------------------------------------------------------------------------------
% Create initial ("true") field with a spike at some point and zeros elsewhere.
% Note, the spike is not very sharp to make the field differentiable.
%--------------------------------------------------------------------------------------------------
function field = InitialField(conf)

    nx = round(conf.nx);
    ny = round(conf.ny);
    field = zeros(nx, ny);
    sigma = round(1);                       % in logical units (point indices)
    cx = round(conf.spot_x / conf.dx);
    cy = round(conf.spot_y / conf.dy);
    a = conf.spot_density / (sigma^2 * 2 * pi);
    b = 1.0 / (2*sigma^2);

    for x=(-4*sigma):(4*sigma)
        ix = min(max(x + cx, 0), nx-1) + 1;
        for y=(-4*sigma):(4*sigma)
            iy = min(max(y + cy, 0), ny-1) + 1;
            field(ix,iy) = field(ix,iy) + a * exp(-b * (x*x + y*y));
        end
    end
    field = ApplyBoundaryCondition(conf, field);

end

%--------------------------------------------------------------------------------------------------
% Using model matrix A, the function integrates advection-diffusion equation forward in time
% (x_{t+1} = A * x_{t}) and records all the solutions - state fields. These fields are considered
% as the "true" state of the nature and the source of the "true" observations.
%--------------------------------------------------------------------------------------------------
function true_fields = ComputeTrueFields(conf, glo_idx)

    figure(1);
    true_fields = zeros(conf.nx, conf.ny, conf.Nt);
    field = InitialField(conf);
    for k=1:conf.Nt
        true_fields(:,:,k) = field;

        fprintf('.');
        if (mod(k,80) == 0), fprintf('\n'); end
        vmin = min(field(:));
        vmax = max(field(:));
        imagesc( transpose( abs( (field - vmin) / max(vmax - vmin, eps) ) ) );
        drawnow();

curr = load(sprintf('../output/curr%05d.txt', k-1));
disp(norm(curr - field)/norm(field));
        
        t = (k - 1) * conf.dt;
        invA = InverseModelMatrix(conf, glo_idx, t);
        
B = load(sprintf('../output/B%05d.txt', k-1));        
disp(norm(B  - invA, 'fro')/norm(invA, 'fro'));
disp(norm(B' - invA, 'fro')/norm(invA, 'fro'));

        field1D = reshape(field, conf.problem_size, 1);
        field = invA \ field1D;
        field = reshape(field, conf.nx, conf.ny);
        field = ApplyBoundaryCondition(conf, field);
        
next = load(sprintf('../output/next%05d.txt', k-1));
disp(norm(next - field)/norm(field));

    end
    fprintf('\n\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filter stuff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------
% Function computes the model noise covariance matrix.
%--------------------------------------------------------------------------------------------------
function Q = ComputeQ(conf)

    N = conf.problem_size;
    Q = spdiags(conf.model_noise_Q * rand(N,1) + 1, 0, N, N);

end

%--------------------------------------------------------------------------------------------------
% Function computes the measurement noise covariance matrix.
%--------------------------------------------------------------------------------------------------
function R = ComputeR(conf)

    N = conf.observation_nx * conf.observation_ny;
    R = spdiags(conf.model_noise_R * rand(N,1) + 1, 0, N, N);

end

%--------------------------------------------------------------------------------------------------
% Function initializes observation matrix H such that observations
% can be available at a (small) subset of domain points.
%--------------------------------------------------------------------------------------------------
function H = ComputeH(conf, glo_idx)

    % Domain and problem sizes, spacial discretization steps.
    nx = conf.nx;
    ny = conf.ny;

    % At least two observations in either dimension, but not more than number of nodal points.
    obs_nx = min(max(round(conf.observation_nx), 2), nx);
    obs_ny = min(max(round(conf.observation_ny), 2), ny);

    num_obs = obs_nx * obs_ny;
    if num_obs < conf.problem_size
        fprintf('WARNING: the number of observations is less than the number of nodal points\n');
        rows = zeros(num_obs,1);
        cols = zeros(num_obs,1);
        vals = ones (num_obs,1);

        % Seed the observation points more or less uniformly over the rectangular domain.
        % TODO: avoid observations at the boundary.
        count = 0;
        for y=1:obs_ny
        for x=1:obs_nx
            count = count + 1;
            rows(count) = count;
            cols(count) = glo_idx(round((x * nx) / double(obs_nx)), ...
                                  round((y * ny) / double(obs_ny)));
        end
        end
        assert(count == num_obs);
        H = sparse(rows, cols, vals, num_obs, conf.problem_size, num_obs);
    else
        H = speye(conf.problem_size);
    end

end

%--------------------------------------------------------------------------------------------------
% Function computes the initial covariance matrix as function of exponential distance.
%--------------------------------------------------------------------------------------------------
function P = InitialCovar(conf, glo_idx)

    nx = conf.nx;
    ny = conf.ny;
    problem_size = conf.problem_size;

    variance = conf.model_ini_var;
    P = zeros(problem_size, problem_size);

    % Here we express correlation distances in logical coordinates of nodal points.
    sx = max(nx * double(conf.model_ini_covar_radius) / conf.domain_size_x, 1);
    sy = max(ny * double(conf.model_ini_covar_radius) / conf.domain_size_y, 1);
    Rx = round(ceil(3*sx));
    Ry = round(ceil(3*sx));

    for v=1:ny
        for u=1:nx
            i = glo_idx(u,v);
            for y=(v-Ry):(v+Ry)
                if (1 <= y) && (y <= ny)
                    dy = (v-y)/sy;
                    for x=(u-Rx):(u+Rx)
                        if (1 <= x) && (x <= nx)
                            dx = (u-x)/sx;
                            j = glo_idx(x,y);
                            covariance = variance * exp(-0.5*(dx*dx + dy*dy));
                            P(i,j) = covariance;
                            P(j,i) = covariance;
                        end
                    end
                end
            end
        end
    end
    P = (P + transpose(P)) * 0.5;           % symmetrize

end

