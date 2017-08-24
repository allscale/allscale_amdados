//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
//             Fearghal O'Donncha, feardonn@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include "allscale/api/user/operator/pfor.h"
#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/geometry.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/sparse_matrix.h"
#include "amdados/app/utils/configuration.h"
#include "amdados/app/utils/cholesky.h"
#include "amdados/app/utils/lu.h"

#undef AMDADOS_ENABLE_GNUPLOT
#include "gnuplot.h"

#include <cstdlib>
#include <unistd.h>
#include <mutex>
#include <cassert>

using ::allscale::api::user::pfor;
using ::amdados::app::utils::Configuration;

// Components to be included
// 1) Grid structures specific to each subdomain
//    Each grid structure contains information
//    a) structures for three different layers of resolution: (100m (1) ; 20m(2); 4m(3))
//    b) Solution on each layer
//
// 2) Mechanism to switch from layers 1, 2 & 3
// 3) Advection diffusion solver translated from api-prototype
//  ) Boundary synchronizations translated from api
//    Iterative check on error convergence across each subdomain
//    until error norm of all 4 boundaries < threshold
// 4) Data assimilation structures translated from api-prototype
// 5) Matrix operation structures for DA
// 6) Data assimilation solution
// 7) File reads for initial conditions (simple assci format)
// 8) File reads for flowfields and observation data (larger files in structured format)

// Model Initialization & file read
// Create structures here for read from file following variables:
// Ndom, nelems;

// Data structures initialization for grid, advection diffusion and DA

// Advection diffusion solution
// Data assimilation check on observation
// Switch to appropriate grid resolution
// Data assimilation solver
// File output at periodic intervals

namespace amdados {
namespace app {

const double TINY = numeric_limits<double>::min() /
           std::pow(numeric_limits<double>::epsilon(),3);

using namespace ::amdados::app::utils;
    //using ::allscale::api::user::pfor;
using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;

// Structure keeps 4-side boundary of a subdomain.
struct Boundary { std::vector<double> side[4]; };

using cube_t         = Grid<double,3>;
using DA_subfield_t  = Matrix<NELEMS_X,NELEMS_Y>;
using DA_vector_t    = Vector<SUB_PROBLEM_SIZE>;
using DA_matrix_t    = Matrix<SUB_PROBLEM_SIZE,SUB_PROBLEM_SIZE>;
using DA_sp_matrix_t = SpMatrix<SUB_PROBLEM_SIZE,SUB_PROBLEM_SIZE>;

using field_grid_t    = Grid<DA_subfield_t,2>;
using vec_grid_t      = Grid<DA_vector_t,2>;
using mat_grid_t      = Grid<DA_matrix_t,2>;
using cholesky_grid_t = Grid<Cholesky<SUB_PROBLEM_SIZE>,2>;
using lu_grid_t       = Grid<LUdecomposition<SUB_PROBLEM_SIZE>,2>;
using boundary_grid_t = Grid<Boundary,2>;

//#################################################################################################
// Utilities.
//#################################################################################################

//-------------------------------------------------------------------------------------------------
// Function rounds the value to the nearest integer.
//-------------------------------------------------------------------------------------------------
inline int Round(double val)
{
    assert_true(std::fabs(val) < numeric_limits<int>::max());
    return static_cast<int>(std::floor(val + 0.5));
}

//-------------------------------------------------------------------------------------------------
// Function creates an output directory inside the current one, which is supposed to be the
// project root folder. If the directory is already exist all its content will be deleted.
// TODO: this will not work on Windows, use STL "experimental" instead.
//-------------------------------------------------------------------------------------------------
void CreateAndCleanOutputDir(const Configuration & conf)
{
    const string & dir = conf.asString("output_dir");
    assert_true(!dir.empty());
    string cmd("mkdir -p ");
    cmd += dir;
    int retval = std::system(cmd.c_str());
    retval = std::system("sync");
    retval = std::system((string("/bin/rm -fr ") + dir + "/*").c_str());
    retval = std::system("sync");
    (void)retval;
}

//-------------------------------------------------------------------------------------------------
// Function maps the subdomain local coordinates (sx,sy) to global ones (gx,gy)
// inside the whole domain given subdomain's position in the grid.
//-------------------------------------------------------------------------------------------------
inline void Sub2Global(int & gx, int & gy, const point2d_t & subdomain, const int sx, const int sy)
{
    assert_true((0 <= sx) && (sx < NELEMS_X));
    assert_true((0 <= sy) && (sy < NELEMS_Y));
    gx = sx + subdomain[0] * NELEMS_X;
    gy = sy + subdomain[1] * NELEMS_Y;
}

//-------------------------------------------------------------------------------------------------
// Functor converts 2D subdomain index (x,y) into a global plain one given point's coordinates
// in a subdomain and subdomain's position in the grid.
//-------------------------------------------------------------------------------------------------
inline int FlatGlobalIndex(int ix, int iy, const point2d_t & subdomain, bool flipY = false)
{
    int gx = 0, gy = 0;
    Sub2Global(gx, gy, subdomain, ix, iy);
    if (flipY) gy = GLOBAL_NELEMS_Y - 1 - gy;
    return Matrix<GLOBAL_NELEMS_X,GLOBAL_NELEMS_Y>::sub2ind(gx, gy);
}

//-------------------------------------------------------------------------------------------------
// Functions for global reduction across all the subdomains.
// TODO: make them MPI friendly.
//-------------------------------------------------------------------------------------------------
double ReduceMean(const Grid<double,2> & grid)
{
    double sum = 0.0;
    for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
    for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { sum += grid[{x,y}]; }}
    return (sum / static_cast<double>(SubDomGridSize[_X_] * SubDomGridSize[_Y_]));
}
double ReduceAbsMin(const Grid<double,2> & grid)
{
    double v = std::fabs(grid[{0,0}]);
    for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
    for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { v = std::min(v, std::fabs(grid[{x,y}])); }}
    return v;
}
double ReduceAbsMax(const Grid<double,2> & grid)
{
    double v = std::fabs(grid[{0,0}]);
    for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
    for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { v = std::max(v, std::fabs(grid[{x,y}])); }}
    return v;
}

//-------------------------------------------------------------------------------------------------
// Function writes the whole property field into a file in binary, grayscaled PGM format,
// where all the values are adjusted to [0..255] interval.
//-------------------------------------------------------------------------------------------------
void WriteImage(const Configuration & conf, const field_grid_t & field,
                const char * title, int time_index,
                std::vector<unsigned char> & image_buffer,
                bool print_image = true, gnuplot::Gnuplot * gp = nullptr)
{
    if (!print_image && (gp == nullptr)) return;
    const bool flipY = true;

    const int ImageSize = GLOBAL_NELEMS_X * GLOBAL_NELEMS_Y;
    if (image_buffer.capacity() < ImageSize) { image_buffer.reserve(ImageSize); }
    image_buffer.resize(ImageSize);

    // Compute the minimum and maximum values of the property field.
    double lower = 0.0, upper = 0.0;
    {
        Grid<double,2> minvalues(SubDomGridSize), maxvalues(SubDomGridSize);
        assert_true((Origin[_X_] == 0) && Origin[_Y_] == 0);
        pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
            minvalues[idx] = *(std::min_element(field[idx].begin(), field[idx].end()));
            maxvalues[idx] = *(std::max_element(field[idx].begin(), field[idx].end()));
        });
        lower = ReduceAbsMin(minvalues);
        upper = ReduceAbsMax(maxvalues);
        assert_true(upper - lower > TINY);
    }

    // Convert the property field into one-byte-per-pixel representation.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        const DA_subfield_t & subfield = field[idx];
        for (int x = 0; x < NELEMS_X; ++x) {
        for (int y = 0; y < NELEMS_Y; ++y) {
            int i = FlatGlobalIndex(x, y, idx, flipY);   // Y-axis points upward
            //bool ok = (i == Matrix<NELEMS_X,NELEMS_Y>::sub2ind(x,y)); assert_true(ok); // XXX
            assert_true((0 <= i) && (i < ImageSize));
            image_buffer[i] = static_cast<unsigned char>(
                    Round(255.0 * (subfield(x,y) - lower) / (upper - lower)) );
        }}
    });

    // Write a file in binary, grayscaled PGM format.
    if (print_image) {
        std::stringstream filename;
        filename << conf.asString("output_dir") << "/" << title
                 << std::setfill('0') << std::setw(5) << time_index << ".pgm";
        std::ofstream f(filename.str(),
                std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        assert_true(f.good()) << "failed to open file for writing: " << filename.str() << endl;
        f << "P5\n" << GLOBAL_NELEMS_X << " " << GLOBAL_NELEMS_Y << "\n" << int(255) << "\n";
        f.write(reinterpret_cast<const char*>(image_buffer.data()), image_buffer.size());
        f << std::flush;
    }

    // Plot image.
    if (gp != nullptr) {
        char title[64];
        snprintf(title, sizeof(title)/sizeof(title[0]), "frame%05d", time_index);
        gp->PlotGrayImage(image_buffer.data(), GLOBAL_NELEMS_X, GLOBAL_NELEMS_Y, title, flipY);
    }
}

//def MakeVideo(conf, filetitle):
    //""" Function creates a single video file from a sequence of field states
        //written into image files.
    //"""
    //print("\n\n\n")
    //framerate = 24
    //if os.system("ffmpeg -y -f image2 -framerate " + str(framerate) + " -pattern_type glob -i '" +
                    //conf.output_dir + "/" + filetitle + "*.png' " +
                    //conf.output_dir + "/" + filetitle + ".avi"):
        //print("WARNING: unable to write video: ffmpeg failed")
    //print("\n\n\n")



//#################################################################################################
// Initialization.
//#################################################################################################

//-------------------------------------------------------------------------------------------------
// Function initializes dependent parameters given the primary ones specified by user.
//-------------------------------------------------------------------------------------------------
void InitDependentParams(Configuration & conf)
{
    // Ensure integer values for certain parameters.
    assert_true(conf.asDouble("num_domains_x")      == conf.asInt("num_domains_x"));
    assert_true(conf.asDouble("num_domains_y")      == conf.asInt("num_domains_y"));
    assert_true(conf.asDouble("num_elems_x")        == conf.asInt("num_elems_x"));
    assert_true(conf.asDouble("num_elems_y")        == conf.asInt("num_elems_y"));
    assert_true(conf.asDouble("observation_nx")     == conf.asInt("observation_nx"));
    assert_true(conf.asDouble("observation_ny")     == conf.asInt("observation_ny"));
    assert_true(conf.asDouble("integration_nsteps") == conf.asInt("integration_nsteps"));

    // Check the geometry. In C++ version the whole domain is divided into subdomains.
    assert_true(conf.asInt("num_domains_x") == NUM_DOMAINS_X) << "num_domains_x mismatch" << endl;
    assert_true(conf.asInt("num_domains_y") == NUM_DOMAINS_Y) << "num_domains_y mismatch" << endl;
    assert_true(conf.asInt("num_elems_x") == NELEMS_X) << "num_elems_x mismatch" << endl;
    assert_true(conf.asInt("num_elems_y") == NELEMS_Y) << "num_elems_y mismatch" << endl;

    const int nx = GLOBAL_NELEMS_X;
    const int ny = GLOBAL_NELEMS_Y;

    const double D = conf.asDouble("diffusion_coef");
    assert_true(D > 0.0);

    conf.SetInt("problem_size", nx * ny);
    const double dx = conf.asDouble("domain_size_x") / (nx-1);
    const double dy = conf.asDouble("domain_size_y") / (ny-1);
    assert_true((dx > 0) && (dy > 0));
    conf.SetDouble("dx", dx);
    conf.SetDouble("dy", dy);

    // Deduce the optimal time step from the stability criteria.
    const double TINY = numeric_limits<double>::min() /
               std::pow(numeric_limits<double>::epsilon(),3);
    const double dt_base = conf.asDouble("integration_period") /
                            conf.asDouble("integration_nsteps");
    const double max_vx = conf.asDouble("flow_model_max_vx");
    const double max_vy = conf.asDouble("flow_model_max_vy");
    const double dt = std::min(dt_base,
                        std::min( std::min(dx*dx, dy*dy)/(2.0*D + TINY),
                                  1.0/(std::fabs(max_vx)/dx + std::fabs(max_vy)/dy + TINY) ));
    assert_true(dt > 0);
    conf.SetDouble("dt", dt);
    conf.SetInt("Nt", static_cast<int>(std::ceil(conf.asDouble("integration_period") / dt)));
}

//-------------------------------------------------------------------------------------------------
// Create initial ("true") field with a spike at some point and zeros elsewhere.
// Note, the spike is not very sharp to make the field differentiable.
// Here 'field' is used as temporary, easy-to-use container which will eventually
// initialize the 'state'.
// \param  state  multi-layered structure that keeps density fields of all the subdomains.
// \param  field  temporary, easy-to-use container that takes fields of all the subdomains.
// \param  conf   configuration parameters.
//-------------------------------------------------------------------------------------------------
void InitialField(domain_t & state, field_grid_t & field, const Configuration & conf)
{
    // Global coordinates of the density spot centre.
    const int cx = Round(conf.asDouble("spot_x") / conf.asDouble("dx"));
    const int cy = Round(conf.asDouble("spot_y") / conf.asDouble("dy"));
    assert_true((0 <= cx) && (cx < GLOBAL_NELEMS_X) &&
                (0 <= cy) && (cy < GLOBAL_NELEMS_Y))
                << "high-concentration spot is not inside the domain" << endl;

    // Parameters of the global 2D Gaussian model of the spike.
    const double sigma = 1.0;                           // in logical units (point indices)
    const double a = conf.asDouble("spot_density") / (std::pow(sigma,2) * 2.0 * M_PI);
    const double b = 1.0 / (2.0 * std::pow(sigma,2));

    // For all the subdomains initialize the spike distribution by parts.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        DA_subfield_t & subfield = field[idx];
        FillMatrix(subfield, 0.0);
        int gx = 0, gy = 0;
        for (int x = 0; x < NELEMS_X; ++x) {
        for (int y = 0; y < NELEMS_Y; ++y) {
            Sub2Global(gx, gy, idx, x, y);      // get global point indices
            int dx = gx - cx;
            int dy = gy - cy;
            if ((std::abs(dx) <= 4*sigma) && (std::abs(dy) <= 4*sigma)) {
                subfield(x,y) += a * std::exp(-b * (dx*dx + dy*dy));
            }
        }}

        // Initialize all cells on the 100m resolution.
        state[idx].setActiveLayer(L_100m);
        AllscaleFromMatrix(state[idx].getLayer<L_100m>(), field[idx]);
    });
}




//def GlobalIndices(conf):
    //""" Each nodal point gets a unique global index on the grid.
        //Function initializes a 2D array of indices:  index of (x,y) = glo_idx(x,y).
    //"""
    //glo_idx = np.arange(int(conf.problem_size)).reshape((conf.nx, conf.ny))
    //return glo_idx


//#################################################################################################
// Advection-diffusion PDE stuff.
//#################################################################################################

//-------------------------------------------------------------------------------------------------
// Function computes flow components given a time.
//-------------------------------------------------------------------------------------------------
void Flow(double & vx, double & vy, const Configuration & conf, const double t)
{
    double max_vx = conf.asDouble("flow_model_max_vx");
    double max_vy = conf.asDouble("flow_model_max_vy");
    vx = -max_vx * std::sin(0.1 * t / conf.asDouble("integration_period") - M_PI);
    vy = -max_vy * std::sin(0.2 * t / conf.asDouble("integration_period") - M_PI);
}

//-------------------------------------------------------------------------------------------------
// Function initializes inverse matrix of implicit Euler time-integrator:
// B * x_{t+1} = x_{t}, where B = A^{-1} is the matrix returned by this function.
// The matrix must be inverted while iterating forward in time: x_{t+1} = A * x_{t}.
// Note, the matrix we generate here is acting on a subdomain.
//-------------------------------------------------------------------------------------------------
void InverseModelMatrix(DA_matrix_t & B, const Configuration & conf, const double t)
{
    using SF = Matrix<NELEMS_X,NELEMS_Y>;   // sub-field type

    const int Nx = NELEMS_X;    // short-hand aliases
    const int Ny = NELEMS_Y;

    const double D  = conf.asDouble("diffusion_coef");
    const double dx = conf.asDouble("dx");
    const double dy = conf.asDouble("dy");
    const double dt = conf.asDouble("dt");

    const double rho_x = D * dt / std::pow(dx,2);
    const double rho_y = D * dt / std::pow(dy,2);

    const double v0x = 2.0 * dx / dt;
    const double v0y = 2.0 * dy / dt;

    double vx = 0.0, vy = 0.0;
    Flow(vx, vy, conf, t);
    vx = vx / v0x;
    vy = vy / v0y;

    //std::cout << "rho_x: " << rho_x << ", rho_y: " << rho_y << ", v0x: " << v0x << ", v0y: " << v0y << endl;

    // Matrix B is supposed to be a sparse one. For now, since we do not have a fast utility
    // for sparse matrix inversion, we define B as a dense one with many zeros. The matrix
    // is assembled in exactly the same way on each iteration, and there is no need to set
    // unused entries to zero as we did it at the beginning.
#if 0
    // Process the internal subdomain points.
    for (int x = 1; x < Nx-1; ++x) {
    for (int y = 1; y < Ny-1; ++y) {
        int i = SF::sub2ind(x,y);
        B(i,i)                  = 1 + 2*(rho_x + rho_y);
        B(i,SF::sub2ind(x-1,y)) = - vx - rho_x;
        B(i,SF::sub2ind(x+1,y)) = + vx - rho_x;
        B(i,SF::sub2ind(x,y-1)) = - vy - rho_y;
        B(i,SF::sub2ind(x,y+1)) = + vy - rho_y;
    }}

    // Relatively expensive lambda function used to initialize matrix over the border points.
    auto SetEntry = [Nx, Ny, vx, vy, rho_x, rho_y, &B](int x, int y) {
        int i = SF::sub2ind(x,y);
        B(i,i) = 1 + 2*(rho_x + rho_y);

        if (x == 0) {
            B(i,SF::sub2ind(x  ,y)) = - 2*vx - rho_x;
            B(i,SF::sub2ind(x+1,y)) = + 2*vx - rho_x;
        } else if (x+1 == Nx) {
            B(i,SF::sub2ind(x-1,y)) = - 2*vx - rho_x;
            B(i,SF::sub2ind(x  ,y)) = + 2*vx - rho_x;
        } else {
            B(i,SF::sub2ind(x-1,y)) = - vx - rho_x;
            B(i,SF::sub2ind(x+1,y)) = + vx - rho_x;
        }

        if (y == 0) {
            B(i,SF::sub2ind(x,y  )) = - 2*vy - rho_y;
            B(i,SF::sub2ind(x,y+1)) = + 2*vy - rho_y;
        } else if (y+1 == Ny) {
            B(i,SF::sub2ind(x,y-1)) = - 2*vy - rho_y;
            B(i,SF::sub2ind(x,y  )) = + 2*vy - rho_y;
        } else {
            B(i,SF::sub2ind(x,y-1)) = - vy - rho_y;
            B(i,SF::sub2ind(x,y+1)) = + vy - rho_y;
        }
    };

    // Top and bottom boundaries.
    for (int x = 1, y =    0; x < Nx-1; ++x) SetEntry(x, y);
    for (int x = 1, y = Ny-1; x < Nx-1; ++x) SetEntry(x, y);

    // Left and right boundaries.
    for (int x =    0, y = 1; y < Ny-1; ++y) SetEntry(x, y);
    for (int x = Nx-1, y = 1; y < Ny-1; ++y) SetEntry(x, y);

    // Corner points.
    SetEntry(   0,    0);
    SetEntry(Nx-1,    0);
    SetEntry(   0, Ny-1);
    SetEntry(Nx-1, Ny-1);
#endif

    // TODO: this operation can be avoided in case of sparse matrix.
    FillMatrix(B, 0.0);

    // Process the internal subdomain points.
    for (int x = 0; x < Nx; ++x) {
    for (int y = 0; y < Ny; ++y) {
        int i = SF::sub2ind(x,y);

        if ((x == 0) || (x+1 == Nx) || (y == 0) || (y+1 == Ny)) {
            B(i,i) += 1 + 2*(rho_x + rho_y);

            if (x == 0) {
                B(i,SF::sub2ind(x,  y)) += - 2*vx - rho_x;
                B(i,SF::sub2ind(x+1,y)) += + 2*vx - rho_x;
            } else if (x+1 == Nx) {
                B(i,SF::sub2ind(x-1,y)) += - 2*vx - rho_x;
                B(i,SF::sub2ind(x,  y)) += + 2*vx - rho_x;
            } else {
                B(i,SF::sub2ind(x-1,y)) += - vx - rho_x;
                B(i,SF::sub2ind(x+1,y)) += + vx - rho_x;
            }

            if (y == 0) {
                B(i,SF::sub2ind(x,y  )) += - 2*vy - rho_y;
                B(i,SF::sub2ind(x,y+1)) += + 2*vy - rho_y;
            } else if (y+1 == Ny) {
                B(i,SF::sub2ind(x,y-1)) += - 2*vy - rho_y;
                B(i,SF::sub2ind(x,y  )) += + 2*vy - rho_y;
            } else {
                B(i,SF::sub2ind(x,y-1)) += - vy - rho_y;
                B(i,SF::sub2ind(x,y+1)) += + vy - rho_y;
            }
        } else {
            B(i,i) = 1 + 2*(rho_x + rho_y);
            B(i,SF::sub2ind(x-1,y)) = - vx - rho_x;
            B(i,SF::sub2ind(x+1,y)) = + vx - rho_x;
            B(i,SF::sub2ind(x,y-1)) = - vy - rho_y;
            B(i,SF::sub2ind(x,y+1)) = + vy - rho_y;
        }
    }}
}

//-------------------------------------------------------------------------------------------------
// Function collects remote neighboring boundaries for each subdomain skipping the global border.
//-------------------------------------------------------------------------------------------------
void GetRemoteBoundary(boundary_grid_t & boundary, const domain_t & domain)
{
    const point2d_t OR = Origin;            // short-hand alias
    const point2d_t SZ = SubDomGridSize;    // short-hand alias

    for (Direction dir : { Up, Down, Left, Right }) assert_true((0 <= dir) && (dir < 4));

    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        std::vector<double> * side = boundary[idx].side;
        if (idx[_Y_] < SZ[_Y_]-1) side[Up   ] = domain[idx + point2d_t{0,+1}].getBoundary(Down );
        if (idx[_Y_] > OR[_Y_])   side[Down ] = domain[idx + point2d_t{0,-1}].getBoundary(Up   );
        if (idx[_X_] < SZ[_X_]-1) side[Right] = domain[idx + point2d_t{+1,0}].getBoundary(Left );
        if (idx[_X_] > OR[_X_])   side[Left ] = domain[idx + point2d_t{-1,0}].getBoundary(Right);

        assert_true(side[Up   ].empty() || (static_cast<int>(side[Up   ].size()) == NELEMS_X));
        assert_true(side[Down ].empty() || (static_cast<int>(side[Down ].size()) == NELEMS_X));
        assert_true(side[Right].empty() || (static_cast<int>(side[Right].size()) == NELEMS_Y));
        assert_true(side[Left ].empty() || (static_cast<int>(side[Left ].size()) == NELEMS_Y));

if (idx[0]==1 && idx[1]==SZ[1]-1) {
    bool ok = true;
    const auto & subdomain = domain[idx].getLayer<L_100m>();
    double mean = 0;
    for (int x = 0; x < NELEMS_X; ++x) {
    for (int y = 0; y < NELEMS_Y; ++y) { mean += std::fabs(subdomain[{x,y}]); }}

    double up = 0, down=0, right=0, left=0;
    ok = (side[Down].empty() || side[Down].size() == (size_t)NELEMS_X); assert_true(ok);
    ok = (side[Up  ].empty() || side[Up  ].size() == (size_t)NELEMS_X); assert_true(ok);
    for (int x = 0; x < NELEMS_X; ++x) {
        if (!side[Down].empty()) {
            //ok = (side[Down][x] == subdomain[{x,0}]);           assert_true(ok);
            down += subdomain[{x,0}];
        }
        //if (!side[Up  ].empty())
        {
            //ok = (side[Up  ][x] == subdomain[{x,NELEMS_Y-1}]);  assert_true(ok);
            up   += subdomain[{x,NELEMS_Y-1}];
        }
    }
    ok = (side[Left ].empty() || side[Left ].size() == (size_t)NELEMS_Y); assert_true(ok);
    ok = (side[Right].empty() || side[Right].size() == (size_t)NELEMS_Y); assert_true(ok);
    for (int y = 0; y < NELEMS_Y; ++y) {
        if (!side[Left ].empty()) {
            //ok = (side[Left ][y] == subdomain[{0,y}]);          assert_true(ok);
            left  += subdomain[{0,y}];
        }
        if (!side[Right].empty()) {
            //ok = (side[Right][y] == subdomain[{NELEMS_X-1,y}]); assert_true(ok);
            right += subdomain[{NELEMS_X-1,y}];
        }
    }

    std::cout << "mean: " << (mean / (NELEMS_X * NELEMS_Y)) << '\n'
              << "left: " << left/NELEMS_Y << ", right: " << right/NELEMS_Y
              << ", down: " << down/NELEMS_X << ", up: " << up/NELEMS_X << endl << flush;
}

    });
}

//-------------------------------------------------------------------------------------------------
// Function implements the idea of Schwartz method where the boundary values of subdomain
// are get updated depending on flow direction.
//-------------------------------------------------------------------------------------------------
double SchwartzUpdate(const Configuration & conf,
                      boundary_grid_t & boundary, domain_t & domain, double physical_time)
{
    // Collect. remote neighboring boundaries for each subdomain.
    GetRemoteBoundary(boundary, domain);

    Grid<double,2> glo_diff(SubDomGridSize);    // accumulator of differences across the domain
    Grid<double,2> glo_mean(SubDomGridSize);    // accumulator of mean density across the domain

    // For each subdomain, update its boundary by the boundary values
    // of neighboring subdomain, if the density flow comes in.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        const std::vector<double> * side = boundary[idx].side;
        double                      vx = 0.0, vy = 0.0;
        double                      diff = 0.0, mean = 0.0;
        size_t                      diff_scale = 0;

        Flow(vx, vy, conf, physical_time);
        for (Direction dir : { Up, Down, Left, Right }) {
            const auto & remote = side[dir];
            if (remote.empty()) continue;               // skip external borders

            const auto myself = domain[idx].getBoundary(dir);
            assert_true(myself.size() == remote.size());

            for (auto a = remote.begin(), b = myself.begin(); b != myself.end(); ++a, ++b) {
                diff += std::fabs(*a - *b);
                ++diff_scale;
            }

            // Make a normal vector pointing outside the subdomain.
            int normal_x = (dir == Right) ? +1 : ((dir == Left) ? -1 : 0);
            int normal_y = (dir ==    Up) ? +1 : ((dir == Down) ? -1 : 0);

            // Update the boundary points if flow enters the subdomain.
            if (normal_x * vx + normal_y * vy < 0) {
                domain[idx].setBoundary(dir, remote);
            }
        }
        glo_diff[idx] = (diff /= static_cast<double>(std::max(diff_scale, size_t(1))));

        // Compute the average density field over the domain.
        const auto & subdomain = domain[idx].getLayer<L_100m>();
        for (int x = 0; x < NELEMS_X; ++x) {
        for (int y = 0; y < NELEMS_Y; ++y) { mean += std::fabs(subdomain[{x,y}]); }}
        glo_mean[idx] = (mean /= (NELEMS_X * NELEMS_Y));
    });

#if 0
{
for (int x = 0; x < NELEMS_X; ++x) {
for (int y = 0; y < NELEMS_Y; ++y) {
    std::cout << "x: " << x << ", y: " << y
        << ", diff: " << glo_diff[{x,y}] << ", mean: " <<  glo_mean[{x,y}] << endl;
}}
std::cout << endl << endl;
}
#endif

    double rel_error = (ReduceMean(glo_diff) / std::max(ReduceMean(glo_mean), TINY));
    //std::cout << "rel_error: " << rel_error << endl << flush;
    return rel_error;
}

//-------------------------------------------------------------------------------------------------
// Using model matrix A, the function integrates advection-diffusion equation forward in time
// (x_{t+1} = A * x_{t}) and records all the solutions - state fields. These fields are
// considered as the "true" state of the nature and the source of the "true" observations.
//-------------------------------------------------------------------------------------------------
void ComputeTrueFields(const Configuration & conf)
{
    std::cout << "Computing the observations, a.k.a 'true' density fields" << endl << flush;
    assert_true((Origin[_X_] == 0) && (Origin[_Y_] == 0)) << "origin is expected at (0,0)";
    //true_fields = np.zeros((conf.nx, conf.ny, conf.Nt))

    std::unique_ptr<gnuplot::Gnuplot> gp(new gnuplot::Gnuplot(nullptr, "-background yellow"));

    // Global grid structures defined for all the subdomains.
    // The difference between subdomains 'curr[i]' and 'state[i]' is that the latter
    // can handle multi-resolution case and communicate with the neighbour subdomains. We copy
    // 'state[i]' to 'curr[i]' on each iteration, do processing and copy it back.

    mat_grid_t                 B(SubDomGridSize);           // inverse model matrices
    field_grid_t               curr(SubDomGridSize);        // copy of the current subdomain states
    field_grid_t               next(SubDomGridSize);        // copy of the next subdomain states
    lu_grid_t                  lu(SubDomGridSize);          // objects for LU decomposition
    std::vector<unsigned char> image_buffer;                // temporary buffer
    domain_t                   state(SubDomGridSize);       // states of the density fields
    boundary_grid_t            boundary(SubDomGridSize);    // neighboring boundary vaklues

    // Clear remote boundary arrays of each subdomain.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        for (Direction dir : { Up, Down, Left, Right }) { boundary[idx].side[dir].clear(); }
    });
    // Generate the initial density field.
    InitialField(state, curr, conf);

    // Time integration forward in time.
    for (int t = 0, Nt = conf.asInt("Nt"); t < Nt; ++t) {
        //std::cout << '+' << flush;
        WriteImage(conf, curr, "field", t, image_buffer, true, gp.get());

        const double physical_time = t * conf.asDouble("dt");
        do {
WriteImage(conf, curr, "field", t, image_buffer, false, gp.get());
            //std::cout << '.' << flush;
            pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
                MatrixFromAllscale(curr[idx], state[idx].getLayer<L_100m>());

                InverseModelMatrix(B[idx], conf, physical_time);    // compute inverse model matrix
                lu[idx].Init(B[idx]);                               // decompose: B = L * U
                lu[idx].Solve(next[idx], curr[idx]);                // x_{t+1} = B^{-1} * x_{t}
                curr[idx] = next[idx];

                AllscaleFromMatrix(state[idx].getLayer<L_100m>(), next[idx]);
            });
        } while (SchwartzUpdate(conf, boundary,
                                state, physical_time) > conf.asDouble("schwartz_tol"));
    }
    std::cout << endl << endl << endl;
}

#if 0
{
            // 1) update boundaries
            for (Direction dir : { Up, Down, Left, Right }) {

                // skip global boarder
                if (dir == Up    && idx[0] == 0)         continue; // if direction == up and no neighbour to south
                if (dir == Down  && idx[0] == size[0]-1) continue;
                if (dir == Left  && idx[1] == 0)         continue;
                if (dir == Right && idx[1] == size[1]-1) continue;

                // obtain the local boundary
                auto local_boundary = cur.getBoundary(dir);

                // obtain the neighboring boundary
                auto remote_boundary =
                    (dir == Up)   ? A[idx + point2d_t{-1,0}].getBoundary(Down)  : // remote boundary is bottom strip of neighbour
                    (dir == Down) ? A[idx + point2d_t{ 1,0}].getBoundary(Up)    : // remote boundary is top of neighbour
                    (dir == Left) ? A[idx + point2d_t{0,-1}].getBoundary(Right) : // remote boundary is left of domain
                                    A[idx + point2d_t{0, 1}].getBoundary(Left);

                // compute local flow in domain to decide if flow in or out of domain
                if (dir == Down) ny = -1; // flow into domain from below
                if (dir == Left) nx = -1; // flow into domain from left
                double flow_boundary =  nx*flowu + ny*flowv;
                // TODO: scale the boundary vectors to the same resolution

                // compute updated boundary
                assert(local_boundary.size() == remote_boundary.size());
                if (flow_boundary < 0) {  // then flow into domain need to update boundary with neighbour value
                    for(size_t i = 0; i<local_boundary.size(); i++) {
                        // for now, we just take the average
                        // need to update this to account for flow direction (Fearghal)
                        local_boundary[i] = remote_boundary[i];
                    }
                }

//std::cout << "DIRECTION: "
//<< (dir == Up ? "Up" : (dir == Down ? "Down" : (dir == Left ? "Left" : "Right"))) << std::endl;
assert(CheckNoNan(res.getLayer<L_100m>()));

                // update boundary in result
                res.setBoundary(dir,local_boundary);

assert(CheckNoNan(res.getLayer<L_100m>()));
            }

}
#endif


} // end namespace app
} // end namespace allscale

int Amdados2DMain()
{
    std::cout << "***** Amdados2D application *****" << std::endl << std::endl << std::flush;
    using namespace ::amdados::app;
    using namespace ::amdados::app::utils;

    // Read the primary parameters.
    Configuration conf;
    conf.ReadConfigFile("amdados.conf");
    // Initialize the rest of parameters that can be deduced from the primary ones.
    InitDependentParams(conf);
    conf.PrintParameters();
    // Create and clear the output directory.
    CreateAndCleanOutputDir(conf);

    // Computing the observations, a.k.a 'true' density fields
    ComputeTrueFields(conf);

    return EXIT_SUCCESS;
}

//Compute(zero, size_global);
//std::cout << "active layer: " << A[{1,1}].getLayer<L_100m>() << std::endl;  // XXX ???



//long max_x = -1, max_y = -1;

        //// Get the remote neighboring boundaries for each subdomain, skipping the global border.
        //pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
            //for (Direction dir : { Up, Down, Left, Right }) {
                //assert_true((0 <= dir) && (dir < 4));
                //if ((dir == Up)    && (idx[0] == Origin[0]))           continue;
                //if ((dir == Down)  && (idx[0] == SubDomGridSize[0]-1)) continue;
                //if ((dir == Left)  && (idx[1] == Origin[1]))           continue;
                //if ((dir == Right) && (idx[1] == SubDomGridSize[1]-1)) continue;
                //glo_boundary[idx].side[dir] =
                    //(dir == Up)   ? glo_state[idx + point2d_t{-1,0}].getBoundary(Down)  :
                    //(dir == Down) ? glo_state[idx + point2d_t{ 1,0}].getBoundary(Up)    :
                    //(dir == Left) ? glo_state[idx + point2d_t{0,-1}].getBoundary(Right) :
                                    //glo_state[idx + point2d_t{0, 1}].getBoundary(Left);
            //}
//g_mutex.lock();
//max_x = std::max(max_x, idx[0]);
//max_y = std::max(max_y, idx[1]);
//g_mutex.unlock();
        //});

        ////auto & kkk = glo_state[{0,0}];
//auto & hhh = glo_boundary[{3,4}];
//for (Direction dir : { Up, Down, Left, Right }) {
    //std::cout << "dir = " << dir << endl;
    //std::cout << "size: " << hhh.side[dir].size() << endl << endl;
//}
//std::cout << "max_x: " << max_x << ", max_y: " << max_y << endl;
//auto & lll = glo_state[{-100,0}]; std::cout << &lll << endl;
//return;

