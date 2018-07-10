//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Variables and data associated with a subdomain.
//=============================================================================
class SubDomain
{
public:
    Matrix        m_curr_field;   // current subdomain state field
    Matrix        m_next_field;   // next subdomain state field

    KalmanFilter  m_Kalman;       // Kalman filter
    Matrix        m_B;            // inverse model matrix

    Matrix        m_P;            // process model covariance
    Matrix        m_Q;            // process noise covariance
    Matrix        m_H;            // observation matrix
    Matrix        m_R;            // observation noise covariance
    Vector        m_z;            // observation vector

    point_array_t m_sensors;      // sensor locations at the finest resolution
    Matrix        m_observations; // all the measurements at sensors

    LUdecomposition m_LU;         // used for state propagation without sensors

    size2d_t      m_size;         // size of this subdomain
    size2d_t      m_ex_size;      // size of extended subdomain
    size2d_t      m_grid_size;    // grid size in number of subdomains
    point2d_t     m_pos;          // subdomain's position on the grid
    size_t        m_Nt;           // number of time integration steps
    size_t        m_Nsubiter;     // number of sub-iterations

private:
    // Neighbour subdomain to this one.
    struct Neighbour {
        int  rank;      // rank of the neighbour's process
        int  pos;       // neighbour's flat index position on the grid
        int  dir;       // identifier of common boundary on neighbour side
        int  tag;       // identifies neighbour and corresponding boundary
        bool insider;   // true if neighbour subdomain is located inside domain

        Neighbour() : rank(-1), pos(-1), dir(-1), tag(-1), insider(false) {}
    };

    int            m_ready_stage;           // indicates initialization stage
    int            m_rank;                  // rank of my process
    long           m_flat_pos;              // my flat index position on grid
    Neighbour      m_neighbour[NSides];     // neighbours of this subdomain
    double_array_t m_send_boundary[NSides]; // boundaries sent to neighbours
    double_array_t m_recv_boundary[NSides]; // received from neighbours
    MPI_Request    m_send_request[NSides];  // statuses of send requests
    size_t         m_send_count;            // number of sent boundaries

public:
//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
SubDomain(const Configuration & conf, const MpiGrid & grid,
          int subdom_rank, point2d_t position)
    : m_curr_field()
    , m_next_field()
    , m_Kalman(), m_B()
    , m_P(), m_Q(), m_H(), m_R(), m_z()
    , m_sensors(), m_observations()
    , m_LU()
    , m_size(), m_ex_size(), m_grid_size(grid.getGridSize()), m_pos(position)
    , m_Nt(0), m_Nsubiter(0)
    , m_ready_stage(0), m_rank(subdom_rank)
    , m_flat_pos(grid.sub2ind(position)), m_neighbour()
    , m_send_boundary(), m_recv_boundary()
    , m_send_request(), m_send_count(0)
{
    // Extended subdomain has extra point layers on either side and these
    // points actually belong to the neighbour subdomains.
    m_ex_size.x = (m_size.x = conf.asInt("subdomain_x")) + 2;
    m_ex_size.y = (m_size.y = conf.asInt("subdomain_y")) + 2;
    assert_true((m_size.x >= 3) && (m_size.y >= 3));

    assert_true(std::max(m_size.x, m_size.y) <= 64) <<
        "recommended subdomain size: not greater than 16x16, maximum: 64x64";
    assert_true(subdom_rank == GetRank());

    m_send_boundary[Left ].resize((size_t)m_size.y);
    m_send_boundary[Right].resize((size_t)m_size.y);

    m_send_boundary[Down ].resize((size_t)m_size.x);
    m_send_boundary[Up   ].resize((size_t)m_size.x);

    m_recv_boundary[Left ].resize((size_t)m_size.y);
    m_recv_boundary[Right].resize((size_t)m_size.y);

    m_recv_boundary[Down ].resize((size_t)m_size.x);
    m_recv_boundary[Up   ].resize((size_t)m_size.x);

    // Neighbour subdomain sees my boundary in opposite way.
    Directions neighbour_boundary_dir[NSides];
    neighbour_boundary_dir[Left]  = Right;
    neighbour_boundary_dir[Right] = Left;
    neighbour_boundary_dir[Down]  = Up;
    neighbour_boundary_dir[Up]    = Down;

    // Initialize structure that keeps information about neighbour subdomains.
    // Note, neighbour subdomains outside the main domain get this process
    // rank, and that fact will be used while sending boundary values.
    for (int i = 0; i < NSides; ++i) {
        m_send_request[i] = 0;
        point2d_t nei_coords = m_pos.neighbour((Directions)i);
        Neighbour & nei = m_neighbour[i];
        if (grid.is_inside(nei_coords)) {
            nei.rank    = grid.getRank(nei_coords);
            nei.pos     = static_cast<int>(grid.sub2ind(nei_coords));
            nei.dir     = static_cast<int>(neighbour_boundary_dir[i]);
            nei.tag     = nei.dir + NSides * nei.pos;
            nei.insider = true;
        } else {
            nei = Neighbour();
        }
    }
}

//-----------------------------------------------------------------------------
// Destructor.
//-----------------------------------------------------------------------------
virtual ~SubDomain()
{
    assert_true(m_ready_stage == 4) << "expects stage 4";
}

//-----------------------------------------------------------------------------
// Returns flat index of a point inside extended (!) subdomain.
//-----------------------------------------------------------------------------
inline long sub2ind_ex(long x, long y) const
{
    return base_sub2ind(x, y, m_ex_size.x, m_ex_size.y);
}

//-----------------------------------------------------------------------------
// Converts flat index inside extended (!) subdomain into point coordinates.
//-----------------------------------------------------------------------------
inline point2d_t ind2sub_ex(long idx) const
{
    return base_ind2sub(idx, m_ex_size.x, m_ex_size.y);
}

//-----------------------------------------------------------------------------
// Function maps a subdomain local coordinates to global ones.
//-----------------------------------------------------------------------------
inline point2d_t sub2glo(const point2d_t & local_point) const
{
    const long x = local_point.x;
    const long y = local_point.y;

    assert_true((0 <= x) && (x < m_size.x) &&
                (0 <= y) && (y < m_size.y));

    return point2d_t(x + m_pos.x * m_size.x,
                     y + m_pos.y * m_size.y);
}

//-----------------------------------------------------------------------------
// Call this function right before loading sensor locations from file.
//-----------------------------------------------------------------------------
virtual void StartSensorsInitialization(const Configuration &)
{
    assert_true(m_ready_stage == 0) << "expects initialization stage 0";
    m_sensors.clear();
    ++m_ready_stage;
}

//-----------------------------------------------------------------------------
// Call this function right after loading sensor locations from file.
//-----------------------------------------------------------------------------
virtual void FinalizeSensorsInitialization(const Configuration & conf)
{
    assert_true(m_ready_stage == 1) << "expects initialization stage 1";
    assert_true(m_observations.Size() == 0) << "no observations at this stage";

    // Mind the extended subdomain: one extra point layer on either side.
    const int sub_problem_size = (int)(m_ex_size.x * m_ex_size.y);
    const int Nsensors = static_cast<int>(m_sensors.size());

    m_Nt = conf.asUInt("Nt");
    m_Nsubiter = conf.asUInt("num_sub_iter");

    m_curr_field.Resize((int)m_ex_size.x, (int)m_ex_size.y);
    m_next_field.Resize((int)m_ex_size.x, (int)m_ex_size.y);
    m_B.Resize(sub_problem_size, sub_problem_size);

    // Note, we initialize the Kalman filter matrices only in the
    // presence of sensor(s), otherwise they are useless.
    if (Nsensors > 0) {
        m_P.Resize(sub_problem_size, sub_problem_size);
        m_Q.Resize(sub_problem_size, sub_problem_size);
        m_H.Resize(Nsensors, sub_problem_size);
        m_R.Resize(Nsensors, Nsensors);
        m_z.Resize(Nsensors);

        // Compute the matrix of observations H. Recall, sensor coordinates
        // were defined for normal subdomain but we operate on extended one,
        // for this reason x+1 and y+1.
        for (size_t k = 0; k < m_sensors.size(); ++k) {
            long x = m_sensors[k].x;    // subdomain-local coordinates
            long y = m_sensors[k].y;
            m_H(static_cast<int>(k),
                static_cast<int>(sub2ind_ex(x + 1, y + 1))) = 1.0;
        }
    }
    ++m_ready_stage;
}

//-----------------------------------------------------------------------------
// Call this function right before loading sensor observations from file.
//-----------------------------------------------------------------------------
virtual void StartObservationsInitialization(const Configuration &)
{
    assert_true(m_ready_stage == 2) << "expects initialization stage 2";
    assert_true(m_observations.Size() == 0) << "no observations at this stage";
    m_observations.Clear();
    ++m_ready_stage;
}

//-----------------------------------------------------------------------------
// Call this function right after loading sensor observations from file.
//-----------------------------------------------------------------------------
virtual void FinalizeObservationsInitialization(const Configuration &)
{
    assert_true(m_ready_stage == 3) << "expects initialization stage 3";
    ++m_ready_stage;
}

//-----------------------------------------------------------------------------
// Function sends the boundary values of this subdomain to its immediate
// neighbours using non-blocking MPI_Isend() routine. The boundary values
// will be copied directly to those subdomains that share the same process
// with this one.
// For example, 5x5 subdomain (7x7 extended subdomain). Internal points 'i'
// belong to subdomain, external points 'e' actually belong to neighbour
// subdomains but we copy their values to extended variant of this subdomain.
// When we send, say, the right boundary point values, they are copied
// to the left column of extended subdomain of the right neighbour:
//
// my extended              my subdomain,        right neighbour sees my
// subdomain:               right boundary:      right boundary on its left:
//
// e e e e e e e            e e e e e e e        e e e e e e e
// e i i i i i e            e i i i i r e        r i i i i i e
// e i i i i i e            e i i i i r e        r i i i i i e
// e i i i i i e            e i i i i r e  --->  r i i i i i e
// e i i i i i e            e i i i i r e        r i i i i i e
// e i i i i i e            e i i i i r e        r i i i i i e
// e e e e e e e            e e e e e e e        e e e e e e e.
//-----------------------------------------------------------------------------
virtual void SendBoundariesToNeighbours(const MpiGrid & grid, long timestamp)
{
    (void)timestamp;

    // Recall, these are normal (not extended) subdomain sizes.
    const int Sx = static_cast<int>(m_size.x);
    const int Sy = static_cast<int>(m_size.y);

//if (tt == 1 && m_pos.x == 0 && m_pos.y == 0) {
//    std::cout << "!!!(0,0)!!!" << std::endl;
//    Print();
//}

    // Get boundaries. N O T E: curr_field occupies extended subdomain.
    for (int x = 0; x < Sx; ++x) {
        m_send_boundary[Down][(size_t)x] = m_curr_field(x + 1, 1);
        m_send_boundary[Up  ][(size_t)x] = m_curr_field(x + 1, Sy);
    }
    for (int y = 0; y < Sy; ++y) {
        m_send_boundary[Left ][(size_t)y] = m_curr_field(1,  y + 1);
        m_send_boundary[Right][(size_t)y] = m_curr_field(Sx, y + 1);
    }

    // Send boundaries to all neighbour subdomains in non-blocking fashion.
    // If a neighbour subdomain is attached to the same process and located
    // inside the domain, then we copy directly.
    m_send_count = 0;
    for (int b = 0; b < NSides; ++b) {
        const Neighbour & nei = m_neighbour[b];
        if (!nei.insider)
            continue;
        if (m_rank != nei.rank) {
            MPI_CHECK(MPI_Isend(static_cast<void*>(m_send_boundary[b].data()),
                                static_cast<int>(m_send_boundary[b].size()),
                                MPI_DOUBLE, nei.rank, nei.tag,
                                MPI_COMM_WORLD,
                                &(m_send_request[m_send_count])));
            ++m_send_count;
        } else {
            SubDomain * sd = grid.getSubdomain(nei.pos);
            sd->m_recv_boundary[nei.dir] = m_send_boundary[b];
        }
    }
}

//-----------------------------------------------------------------------------
// Function receives boundary values of neighbour subdomains sent to this one.
// The function does not return until this subdomain has received the boundary
// values from its neighbour (blocking MPI receive routine).
//-----------------------------------------------------------------------------
virtual void ReceiveBoundariesFromNeighbours(long timestamp)
{
    (void)timestamp;

    // Note, these are normal (not extended) subdomain sizes.
    const int Sx = static_cast<int>(m_size.x);
    const int Sy = static_cast<int>(m_size.y);
    const int tag_base = static_cast<int>(m_flat_pos * NSides);

    auto Receive = [tag_base,this](Directions dir) -> double_array_t & {
        int nei_rank = m_neighbour[dir].rank;
        if (nei_rank != m_rank) {
            int tag = tag_base + (int)dir;
            MPI_CHECK(MPI_Recv(static_cast<void*>(m_recv_boundary[dir].data()),
                               static_cast<int>(m_recv_boundary[dir].size()),
                               MPI_DOUBLE, nei_rank, tag, MPI_COMM_WORLD,
                               MPI_STATUS_IGNORE));
        }
        return m_recv_boundary[dir];
    };

    Matrix & field = m_curr_field;  // short-hand alias
//if (tt == 1 && m_pos.x == 0 && m_pos.y == 1) {
//    Print();
//}
    // Set up left-most points from the right boundary of the left neighbour.
    if (m_pos.x > 0) {
        double_array_t & boundary = Receive(Left);
        for (int y = 0; y < Sy; ++y) field(0, y+1) = boundary[(size_t)y];
    } else {
        for (int y = 0; y < Sy; ++y) field(0, y+1) = field(2, y+1);
    }
//if (tt == 1 && m_pos.x == 0 && m_pos.y == 1) {
//    Print();
//}
    // Set up right-most points from the left boundary of the right neighbour.
    if (m_pos.x + 1 < m_grid_size.x) {
        double_array_t & boundary = Receive(Right);
        for (int y = 0; y < Sy; ++y) field(Sx+1, y+1) = boundary[(size_t)y];
    } else {
        for (int y = 0; y < Sy; ++y) field(Sx+1, y+1) = field(Sx-1, y+1);
    }
//if (tt == 1 && m_pos.x == 0 && m_pos.y == 1) {
//    Print();
//}
    // Set up bottom-most points from the top boundary of the bottom neighbour.
    if (m_pos.y > 0) {
        double_array_t & boundary = Receive(Down);
        for (int x = 0; x < Sx; ++x) field(x+1, 0) = boundary[(size_t)x];
    } else {
        for (int x = 0; x < Sx; ++x) field(x+1, 0) = field(x+1, 2);
    }
//if (tt == 1 && m_pos.x == 0 && m_pos.y == 1) {
//    Print();
//}
    // Set up top-most points from the bottom boundary of the top neighbour.
    if (m_pos.y + 1 < m_grid_size.y) {
        double_array_t & boundary = Receive(Up);
        for (int x = 0; x < Sx; ++x) field(x+1, Sy+1) = boundary[(size_t)x];
    } else {
        for (int x = 0; x < Sx; ++x) field(x+1, Sy+1) = field(x+1, Sy-1);
    }
//if (tt == 1 && m_pos.x == 0 && m_pos.y == 1) {
//    Print();
//}
//if (tt == 1 && m_pos.x == 0 && m_pos.y == 0) {
//    std::cout << "!!!(0,0)!!!" << std::endl;
//    Print();
//}
    // The corner points are not used in finite-difference scheme applied to
    // internal subdomain points, however, Kalman filter uses the entire
    // subdomain field (as currently implemented but could be avoided).
    // For the latter reason, we have to assign some reasonable values to the
    // corner points. Note, since we operate on extended subdomains, their
    // corners belong to unreachable diagonal neighbours.
    field(0,0)       = (field(1,0)     + field(1,1)   + field(0,1)    ) / 3.0;
    field(0,Sy+1)    = (field(1,Sy+1)  + field(1,Sy)  + field(0,Sy)   ) / 3.0;
    field(Sx+1,0)    = (field(Sx,0)    + field(Sx,1)  + field(Sx+1,1) ) / 3.0;
    field(Sx+1,Sy+1) = (field(Sx,Sy+1) + field(Sx,Sy) + field(Sx+1,Sy)) / 3.0;
}

//-----------------------------------------------------------------------------
// Function waits until all neighbour subdomains confirmed that they have
// received the boundary values sent by this subdomain.
//-----------------------------------------------------------------------------
virtual void WaitForExchangeCompletion()
{
    MPI_Status status[NSides];
    MPI_CHECK(MPI_Waitall((int)m_send_count, m_send_request, status));
}

//-----------------------------------------------------------------------------
// Debugging function.
//-----------------------------------------------------------------------------
void Print() const
{
    for (int x = 0; x < m_ex_size.x; ++x) {
        for (int y = 0; y < m_ex_size.y; ++y) {
            std::cout << "  " << std::setw(6) << m_curr_field(x, y);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

};  // class SubDomain

}   // namespace amdados

