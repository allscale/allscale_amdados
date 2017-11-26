//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {
namespace utils {

using ::amdados::app::utils::Configuration;

#ifdef AMDADOS_ENABLE_AVERAGE_PROFILE

//=================================================================================================
// This is debugging/testing/visualization facility.
// Class accumulates some entity over time in per-subdomain fashion. At the end of time
// progression, the entities are aggregated across all the subdomains and the average profile
// (as a function of time) is printed into a text file.
//=================================================================================================
class AverageProfile
{
public:
    // Constructor.
    AverageProfile(const Configuration & conf) : m_accums(SubDomGridSize)
    {
        Clear(conf);
    }

    // Accumulate differences as the process progresses over time.
    void Accumulate(const point2d_t & idx, double diff)
    {
        m_accums[idx].push_back(diff);
    }

    // Print the aggregated (over all the subdomains) profile of the average entity
    // as a function of time (iteration number).
    void PrintProfile(const Configuration & conf, const char * entity_name)
    {
        assert_true(Origin == point2d_t(0,0));
        const size_t len = m_accums[Origin].size();
        std::vector<double> profile(len);
        for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
        for (int y = 0; y < SubDomGridSize[_Y_]; ++y) {
            const auto & acc = m_accums[{x,y}];
            assert_true(acc.size() == len)
                << "difference profiles of all the subdomains must have the same length" << endl;
            std::transform(profile.begin(), profile.end(), acc.begin(),
                           profile.begin(), std::plus<double>());
        }}

        assert_true(entity_name != nullptr);
        std::fstream f(conf.asString("output_dir") + "/" + entity_name + ".txt",
                       std::ios::out | std::ios::trunc);
        assert_true(f.good());
        for (auto v : profile) f << v/(SubDomGridSize[_X_]*SubDomGridSize[_Y_]) << endl;
        f << flush;
    }

    // Clear this object, so it can be resused.
    void Clear(const Configuration & conf)
    {
        using ::allscale::api::user::algorithm::pfor;
        pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
            auto & acc = m_accums[idx];
            acc.clear();
            acc.reserve(conf.asInt("Nt"));
        });
    }

private:
    ::allscale::api::user::data::Grid<std::vector<double>,2> m_accums; // accumulator per subdomain

}; // class AverageProfile

#else // !AMDADOS_ENABLE_AVERAGE_PROFILE

//=================================================================================================
// Stub AverageProfile class does nothing.
//=================================================================================================
class AverageProfile
{
public:
    AverageProfile(const Configuration &) {
        std::cout << "AverageProfile is disabled" << std::endl;
    }
    void Accumulate(const point2d_t &, double) {}
    void PrintProfile(const Configuration &, const char *) {}
    void Clear(const Configuration &) {}
};

#endif // AMDADOS_ENABLE_AVERAGE_PROFILE

} // end namespace utils
} // end namespace app
} // end namespace amdados


