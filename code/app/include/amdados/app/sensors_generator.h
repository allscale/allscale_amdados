//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// This class is used as a namespace with virtual functions, which are never
// inlined or duplicated. The latter allows us to keep implementation in the
// header file included in both Allscale and MPI projects.
//=============================================================================
class SensorsGenerator
{
public:
// Destructor.
virtual ~SensorsGenerator() {}

//-----------------------------------------------------------------------------
// Function generates pseudo-randomly distributed sensor locations.
// N O T E, the sensor points are stored unordered and their coordinates are
// global, i.e. defined with respected to the whole domain's coordinate system.
//-----------------------------------------------------------------------------
virtual void MakeSensors(const Configuration & conf,
                         std::vector<point2d_t> & sensors) const
{
    MY_TIME_IT("Running sensors' generator ...")

    const size_t Sx = conf.asUInt("subdomain_x");
    const size_t Sy = conf.asUInt("subdomain_y");

    const size_t Nx = conf.asUInt("num_subdomains_x");
    const size_t Ny = conf.asUInt("num_subdomains_y");
    const size_t Nsubdom = Nx * Ny;                 // number of subdomains
    const size_t Ntotal = (Sx * Sy) * (Nx * Ny);    // total number of points

    // Compute the number of sensors.
    // TODO: here we define the number
    // of sensors as a fraction of the total number of nodal points. On the
    // other hand, we want not more than one sensor per subdomain. Should
    // we instead count the fraction in a number of subdomains having sensor?
    const double fraction =
            std::min(std::max(conf.asDouble("sensor_fraction"), 0.001), 0.75);
    const size_t Nsensors = std::min(std::max(static_cast<size_t>(
            std::floor(fraction * Ntotal + 0.5)), size_t(1)), Nsubdom);
    assert_true((1 <= Nsensors) && (Nsensors <= Nsubdom));
    MY_LOG(INFO) << "fraction of sensor points = " << fraction
                 << ", #sensors = " << Nsensors
                 << ", #subdomains = " << Nsubdom
                 << ", #nodal points = " << Ntotal;

    // Random permutation of flat subdomain indices.
    std::vector<size_t> perm(Nsubdom);
    for (size_t i = 0; i < Nsubdom; ++i) { perm[i] = i; }
    std::shuffle(perm.begin(), perm.end(),
            std::default_random_engine(static_cast<unsigned>(RandomSeed())));

    // Randomly select a sub-set of subdomains with a single sensor therein.
    std::srand(static_cast<unsigned>(RandomSeed()));
    sensors.resize(Nsensors);
    for (size_t i = 0; i < Nsensors; ++i) {
        // Get positions of a randomly selected subdomain.
        size_t xs = perm[i] % Nx;       // i = xs + ys * Nx
        size_t ys = perm[i] / Nx;

        // Get randomly selected point (sensor location) inside the subdomain.
        size_t x = static_cast<size_t>(std::rand());
        size_t y = static_cast<size_t>(std::rand());
        x = (Sx >= 3) ? (1 + (x % (Sx - 2))) : (x % Sx);
        y = (Sy >= 3) ? (1 + (y % (Sy - 2))) : (y % Sy);

        // Store sensor location in global coordinates.
        sensors[i] = Sub2Glo(point2d_t(x,y), point2d_t(xs,ys), size2d_t(Sx,Sy));
    }
}

};  // class SensorsGenerator

}   // namespace amdados

