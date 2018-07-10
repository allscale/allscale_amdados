//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

namespace amdados {

//=============================================================================
// Class loads sensors and measurements. The rationale for this class is to
// define non-inline (virtual) functions in the header file.
//=============================================================================
class InputData
{
public:
//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
InputData() {}

//-----------------------------------------------------------------------------
// Destructor.
//-----------------------------------------------------------------------------
virtual ~InputData() {}

//-----------------------------------------------------------------------------
// Function reads sensor locations and all the observations associated with
// the subdomains attached to this process.
//-----------------------------------------------------------------------------
virtual void Load(const Configuration & conf, MpiGrid & grid) const
{
    LoadSensorLocations(conf, grid);
    LoadSensorMeasurements(conf, grid);
}

private:
//-----------------------------------------------------------------------------
// Function reads the file of sensor locations and places their coordinates
// into subdomains attached to this process.
//-----------------------------------------------------------------------------
virtual void LoadSensorLocations(const Configuration & conf, MpiGrid & grid)
const
{
    const size2d_t grid_size(conf.asInt("num_subdomains_x"),
                             conf.asInt("num_subdomains_y"));
    assert_true(grid_size == grid.getGridSize());

    const size2d_t cell_size(conf.asInt("subdomain_x"),
                             conf.asInt("subdomain_y"));

    // Pre-initialize sensors of subdomains attached to this process.
    grid.forAllLocal([&conf](SubDomain * sd) {
        sd->StartSensorsInitialization(conf);
    });

    // Read the file and add sensors to subdomains attached to this process.
    std::string filename = MakeFileName(conf, "sensors");
    CheckFileExists(conf, filename);
    std::fstream f(filename, std::ios::in);
    assert_true(f.good()) << "failed to open sensors' file: " << filename;
    for (point2d_t pt; f >> pt.x >> pt.y;) {
        point2d_t pos = Glo2Cell(pt, cell_size);
        assert_true((0 <= pos.x) && (pos.x < grid_size.x));
        assert_true((0 <= pos.y) && (pos.y < grid_size.y));
        SubDomain * sd = grid.getSubdomain(pos);
        if (sd != nullptr) {
            sd->m_sensors.push_back(Glo2Sub(pt, cell_size));
        }
    }

    // Post-initialize sensors of subdomains attached to this process.
    grid.forAllLocal([&conf](SubDomain * sd) {
        sd->FinalizeSensorsInitialization(conf);
    });
    MY_LOG(INFO) << "Sensor locations have been loaded";

//#ifdef AMDADOS_DEBUGGING   TODO: mpi reduce mean
//    size_t num_sensors = 0;
//    sensors.forEach([&](point_array_t & arr) { num_sensors += arr.size(); });
//    MY_LOG(INFO) << "Average number of sensors per subdomain: " <<
//            (double(num_sensors) / double(MpiGridSize[0] * MpiGridSize[1]));
//#endif
}

//-----------------------------------------------------------------------------
// Function reads the file of sensor measurements.
// N O T E, this function must be invoked after \sa LoadSensorLocations().
//-----------------------------------------------------------------------------
virtual void LoadSensorMeasurements(const Configuration & conf, MpiGrid & grid)
const
{
    MY_LOG(INFO) << "\n"
        << "--------------------------------------------------------\n"
        << "B E W A R E: after regeneration of sensors through the\n"
        << "scenario 'sensors' you have to compute new observations:\n"
        << "\t python3 python/ObservationsGenerator.py\n"
        << "because observations stick to sensor locations.\n"
        << "--------------------------------------------------------\n";

    const size2d_t grid_size(conf.asInt("num_subdomains_x"),
                             conf.asInt("num_subdomains_y"));
    assert_true(grid_size == grid.getGridSize());

    const size2d_t cell_size(conf.asInt("subdomain_x"),
                             conf.asInt("subdomain_y"));

    const int Nt = conf.asInt("Nt");

    // Pre-initialize observations of subdomains attached to this process.
    grid.forAllLocal([&conf](SubDomain * sd) {
        sd->StartObservationsInitialization(conf);
    });

    // Read observations for subdomains attached to this process.
    std::string filename = MakeFileName(conf, "analytic");
    CheckFileExists(conf, filename);
    std::fstream f(filename, std::ios::in);
    assert_true(f.good()) << "failed to open observations' file: " << filename;

    std::vector<int> counters((size_t)(grid_size.x * grid_size.y));

    // Read the file of all the observations by time slices.
    int last_timestamp = -1;
    while (1) {
        int t = 0, num = 0;
        point2d_t pt;
        float val = 0.0f;

        // Read the header of a new time-slice (timestamp and num. of records).
        if (!(f >> t >> num)) {
            break;
        }
        assert_true(t < Nt);
        assert_true(last_timestamp + 1 == t);

        // Reset all the counters upon arrival of a new time-slice.
        std::fill(counters.begin(), counters.end(), 0);

        // Read all the records of the time-slice.
        for (int i = 0; i < num; ++i) {
            // Read the global coordinates of a sensor and a measurement.
            if (!(f >> pt.x >> pt.y >> val)) {
                assert_true(0) << "failed to read observations' file";
            }

            // Get subdomain position on the grid.
            point2d_t pos = Glo2Cell(pt, cell_size);
            assert_true((0 <= pos.x) && (pos.x < grid_size.x));
            assert_true((0 <= pos.y) && (pos.y < grid_size.y));

            // Get the data matrix of the subdomain, and the counter.
            SubDomain * sd = grid.getSubdomain(pos);
            if (sd == nullptr) {
                continue;           // skip foreign subdomains
            }

            Matrix & obs = sd->m_observations;
            int    & cnt = counters[(size_t)grid.sub2ind(pos)];
            int      nsensors = static_cast<int>(sd->m_sensors.size());

            // Allocate the data matrix, if necessary.
            // It is fine if nsensors == 0.
            if (obs.Empty()) {
                obs.Resize(Nt, nsensors);
            }

            // Check that sensor coordinates appear in exactly the same order
            // as in array of sensor locations.
            assert_true(cnt < nsensors);
            assert_true(Glo2Sub(pt, cell_size) == sd->m_sensors[(size_t)cnt]);

            // Save the measurement in the data matrix.
            obs(t,cnt) = val;
            ++cnt;
        }
        last_timestamp = t;

        // Check that all the entries of data matrix had been set.
        grid.forAllLocal([&counters,&grid](SubDomain * sd) {
            assert_true(sd != nullptr);
            assert_true(counters[(size_t)grid.sub2ind(sd->m_pos)] ==
                                     static_cast<int>(sd->m_sensors.size()));
        });
    }
    assert_true(last_timestamp + 1 == Nt);

    // Post-initialize observations of subdomains attached to this process.
    grid.forAllLocal([&conf](SubDomain * sd) {
        sd->FinalizeObservationsInitialization(conf);
    });
    MY_LOG(INFO) << "Sensor observations have been loaded";
}

//-----------------------------------------------------------------------------
// Function converts global coordinates to the subdomain index on the grid.
//-----------------------------------------------------------------------------
static point2d_t Glo2Cell(const point2d_t & p, const size2d_t & cell_size)
{
    return point2d_t(p.x / cell_size.x,
                     p.y / cell_size.y);
}

//-----------------------------------------------------------------------------
// Function converts global coordinates to the local ones inside a subdomain.
//-----------------------------------------------------------------------------
static point2d_t Glo2Sub(const point2d_t & p, const size2d_t & cell_size)
{
    return point2d_t(p.x % cell_size.x,
                     p.y % cell_size.y);
}

};  // class InputData

}   // namespace amdados
