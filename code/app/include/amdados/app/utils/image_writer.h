#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

namespace amdados {
namespace app {
namespace utils {

//=============================================================================
// Class (member function) converts a state field into image and
// saves it into "pgm" file.
//=============================================================================
class ImageWriter
{
private:
    std::string                m_output_dir;        // output directory
    std::vector<unsigned char> m_image_buffer;      // temporary image buffer

public:
//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
explicit ImageWriter(const std::string & output_dir) : m_output_dir(output_dir)
{
    if (m_output_dir.empty()) m_output_dir.append(".");
}

//-----------------------------------------------------------------------------
// Function writes the whole property field into a file in binary, gray-scaled
// PGM format, where all the values are adjusted to [0..255] interval.
// A T T E N T I O N: the function must be used outside pfor(..) loop.
// \param  title        file title (extension will be "pgm").
// \param  field        global grid of subdomains to be saved as an image.
// \param  time_index   discrete timestamp.
// \param  write_image  true, if saving into file was intended.
// \return  reference to image that can be plotted, for example, by Gnuplot.
//-----------------------------------------------------------------------------
const std::vector<unsigned char> & Write(const char * title,
        const ::allscale::api::user::data::Grid<Matrix,2> & field,
        int time_index, bool write_image)
{
	using ::allscale::api::user::algorithm::pfor;
	
    const double TINY = numeric_limits<double>::min() /
               std::pow(numeric_limits<double>::epsilon(),3);
    const bool flipY = false;

        const int ImageSize = GLOBAL_NELEMS_X * GLOBAL_NELEMS_Y;
        if (m_image_buffer.capacity() < ImageSize) {
            m_image_buffer.reserve(ImageSize);
        }
        m_image_buffer.resize(ImageSize);

    // Compute the minimum and maximum values of the property field.
    double lower = 0.0, upper = 0.0;
    {
        ::allscale::api::user::data::Grid<double,2> minvalues(SubDomGridSize),
                                                    maxvalues(SubDomGridSize);
        pfor(Origin, SubDomGridSize,
            [&](const point2d_t & idx) {
                minvalues[idx] = *(std::min_element(field[idx].begin(),
                                                    field[idx].end()));
                maxvalues[idx] = *(std::max_element(field[idx].begin(),
                                                    field[idx].end()));
            });
        lower = ReduceAbsMin(minvalues);
        upper = ReduceAbsMax(maxvalues);
        //assert_true(upper - lower > TINY);
    }

    // Convert the property field into one-byte-per-pixel representation.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        const auto & subfield = field[idx];
        for (int y = 0; y < NELEMS_Y; ++y) { int yg = Sub2GloY(idx, y);
        for (int x = 0; x < NELEMS_X; ++x) { int xg = Sub2GloX(idx, x);
            int pos = xg +
                GLOBAL_NELEMS_X * (flipY ? yg : (GLOBAL_NELEMS_Y - 1 - yg));
            if (!(static_cast<unsigned>(pos) <
                  static_cast<unsigned>(ImageSize))) assert_true(0);
            m_image_buffer[pos] = static_cast<unsigned char>(
                    Round(255.0 * (subfield(x,y) - lower) /
                            std::max(upper - lower, TINY)) );
        }}
    });

    // Write a file in binary, gray-scaled PGM format.
    if (write_image) {
        std::stringstream filename;
        filename << m_output_dir << "/" << title << std::setfill('0')
                 << std::setw(5) << time_index << ".pgm";
        std::ofstream f(filename.str(), std::ios_base::out    |
                                        std::ios_base::binary |
                                        std::ios_base::trunc);
        assert_true(f.good()) << "failed to open file for writing: "
                              << filename.str() << endl;
        f << "P5\n" << GLOBAL_NELEMS_X << " "
                    << GLOBAL_NELEMS_Y << "\n" << int(255) << "\n";
        f.write(reinterpret_cast<const char*>(m_image_buffer.data()),
                m_image_buffer.size());
        f << std::flush;
    }
    return m_image_buffer;
}

}; // ImageWriter

} // end namespace utils
} // end namespace app
} // end namespace allscale

