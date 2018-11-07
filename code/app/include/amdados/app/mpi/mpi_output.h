//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Class writes out the current solution field into a file.
//=============================================================================
class MpiOutputWriter
{
public:
//-----------------------------------------------------------------------------
// Constructor open the output file for writing.
//-----------------------------------------------------------------------------
MpiOutputWriter(const Configuration & conf)
    : m_file(MPI_FILE_NULL)
    , m_buffer()
    , m_base(0)
    , m_info(MPI_INFO_NULL)
{
    MPI_Info_create(&m_info);
    /* Disables ROMIO's data-sieving */
    MPI_Info_set(m_info, "romio_ds_read", "disable");  
    MPI_Info_set(m_info, "romio_ds_write", "disable");

    std::string filename = MakeFileName(conf, "field");
    if (MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                      MPI_MODE_CREATE | MPI_MODE_WRONLY,
                      m_info, &m_file) != MPI_SUCCESS) {
        MY_LOG(ERROR) << "failed to open file " << filename << " for writing";
        PrintErrorAndExit(nullptr);
    }
    MPI_CHECK(MPI_File_set_size(m_file, 0));
    MPI_CHECK(MPI_File_get_position(m_file, &m_base));
}

//-----------------------------------------------------------------------------
// Destructor closes the output file.
//-----------------------------------------------------------------------------
virtual ~MpiOutputWriter()
{
    // Collective close operation.
    if (MPI_File_close(&m_file) != MPI_SUCCESS) {
        MY_LOG(ERROR) << "MPI_File_close() failed";
        PrintErrorAndExit(nullptr);
    }
    MPI_Info_free(&m_info);
}

//-----------------------------------------------------------------------------
// Function appends the current solution field to the file.
//-----------------------------------------------------------------------------
virtual void Write(long timestamp, const MpiGrid & grid)
{
    // Allocate temporary buffer once at the beginning. Although we recompute
    // the buffer size every time, this should be a little overhead.
    // Recall, field occupies extended subdomain.
    {
        size_t buffer_size = 0;
        grid.forAllLocal([&](const SubDomain * sd) {
            buffer_size += 4 * static_cast<size_t>(sd->m_size.x * sd->m_size.y);
            assert_true((sd->m_curr_field.NRows() == sd->m_size.x + 2) &&
                        (sd->m_curr_field.NCols() == sd->m_size.y + 2));
        });
        if (m_buffer.empty()) m_buffer.resize(buffer_size);
        assert_true(m_buffer.size() == buffer_size);
    }

    // Fill up 4-column binary buffer: (time, abscissa, ordinate, value).
    // Recall, field occupies extended subdomain so we take field(x+1,y+1)
    // rather than field(x,y).
    size_t idx = 0;
    grid.forAllLocal([&idx,timestamp,this](const SubDomain * sd) {
        const Matrix & field = sd->m_curr_field;
        for (int x = 0; x < sd->m_size.x; ++x) {
        for (int y = 0; y < sd->m_size.y; ++y) {
            point2d_t glo = sd->sub2glo(point2d_t(x,y)); // global coordinates
            m_buffer[idx + 0] = static_cast<float>(timestamp);
            m_buffer[idx + 1] = static_cast<float>(glo.x);
            m_buffer[idx + 2] = static_cast<float>(glo.y);
            m_buffer[idx + 3] = static_cast<float>(field(x + 1, y + 1));
            idx += 4;
        }}
    });
    assert_true(m_buffer.size() == idx);
    assert_true(m_buffer.size() < (unsigned)std::numeric_limits<int>::max());

    // Collectively write into the output file.
//    if (MPI_File_write_all(m_file, m_buffer.data(),
//                           static_cast<int>(m_buffer.size()),
//                           MPI_FLOAT, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
//        MY_LOG(ERROR) << "failed to write in the output file";
//        PrintErrorAndExit(nullptr);
//    }

    MPI_Offset len = static_cast<MPI_Offset>(m_buffer.size() * sizeof(float));
    MPI_Offset offset = 0;
    MPI_CHECK(MPI_Exscan(&len, &offset, 1,
                         MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD));

    if (MPI_File_write_at_all(m_file, m_base + offset,
                              m_buffer.data(),
                              static_cast<int>(m_buffer.size()),
                              MPI_FLOAT, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        MY_LOG(ERROR) << "failed to write in the output file";
        PrintErrorAndExit(nullptr);
    }

    MPI_Offset global_len = 0;
    MPI_CHECK(MPI_Allreduce(&len, &global_len, 1,
                            MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD));
    m_base += global_len;
}

private:
    MPI_File      m_file;           ///< output file
    float_array_t m_buffer;         ///< temporary buffer
    MPI_Offset    m_base;           ///< position in the file
    MPI_Info      m_info;			///< additional IO settings

};  // class MpiOutputWriter

}   // namespace amdados

