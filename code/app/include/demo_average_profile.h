//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

using ::allscale::api::user::data::Grid;

#ifdef AMDADOS_DEBUGGING

point2d_t GetGridSize(const Configuration & conf);

/**
 * This is debugging/testing/visualization facility. Class accumulates some
 * entity over time in per-subdomain fashion. At the end of time progression,
 * the entities are aggregated across all the subdomains and the average
 * profile (as a function of time) is printed into a text file.
 */
class AverageProfile
{
public:
/**
 * Constructor.
 */
AverageProfile(const Configuration & conf) : m_accums(GetGridSize(conf))
{
    Clear(conf);
}

/**
 * Destructor. Virtual functions is a portable way to avoid inlining.
 */
virtual ~AverageProfile() {}

/**
 * Accumulate differences as the process progresses over time.
 */
virtual void Accumulate(const point2d_t & idx, double diff)
{
    m_accums[idx].push_back(diff);
}

/**
 * Print the aggregated (over all the subdomains) profile of the average entity
 * as a function of time (iteration number).
 */
virtual void PrintProfile(const Configuration & conf, const char * entity_name)
{
    const size_t len = m_accums[{0,0}].size();
    double_array_t profile(len);
    // Doing reduction sequentially.
    m_accums.forEach([&](const double_array_t & acc) {
        assert_true(acc.size() == len)
            << "difference profiles of all the subdomains must have "
            << "the same length" << std::endl;
        std::transform(profile.begin(), profile.end(), acc.begin(),
                       profile.begin(), std::plus<double>());
    });

    assert_true(entity_name != nullptr);
    const int Nx = conf.asInt("num_subdomains_x") * conf.asInt("subdomain_x");
    const int Ny = conf.asInt("num_subdomains_y") * conf.asInt("subdomain_y");
    std::stringstream filename;
    filename << conf.asString("output_dir") << PathSep << entity_name
             << "_Nx" << Nx << "_Ny" << Ny << ".txt";

    std::fstream f(filename.str(), std::ios::out | std::ios::trunc);
    assert_true(f.good()) << "failed to open file for writing: "
                          << filename.str();
    // Scale the profile to get average over all the subdomains.
    const point2d_t size = m_accums.size();
    const double scale = 1.0 / static_cast<double>(size.x * size.y);
    for (auto v : profile) f << (v * scale) << std::endl;
    f.flush();
}

/**
 * Clear this object, so it can be reused.
 */
virtual void Clear(const Configuration & conf)
{
    using ::allscale::api::user::algorithm::pfor;
    pfor(point2d_t(0,0), m_accums.size(), [&](const point2d_t & idx) {
        auto & acc = m_accums[idx];
        acc.clear();
        acc.reserve(conf.asInt("Nt"));
    });
}

private:
    Grid<double_array_t,2> m_accums;    // accumulator per subdomain

}; // class AverageProfile

#else   // !AMDADOS_DEBUGGING

/**
 * Stub AverageProfile class does nothing in production mode.
 */
class AverageProfile
{
public:
    AverageProfile(const Configuration &) {}
    virtual ~AverageProfile() {}
    virtual void Accumulate(const point2d_t &, double) {}
    virtual void PrintProfile(const Configuration &, const char *) {}
    virtual void Clear(const Configuration &) {}
};

#endif  // AMDADOS_DEBUGGING

} // namespace amdados


