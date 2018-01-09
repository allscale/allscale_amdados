//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

#include "demo_average_profile.h"
#include "demo_gnuplot.h"

namespace amdados {

using ::allscale::api::user::data::Grid;

#ifdef AMDADOS_DEBUGGING

/**
 * Class (member function) converts a state field into image and
 * saves it into "pgm" file.
 */
class Visualizer
{
private:
    std::string      m_output_dir;  // output directory
    ubyte_array_t    m_image;       // image buffer
    gnuplot::Gnuplot m_gnuplot;     // image plotter

public:
/**
 * Constructor.
 */
explicit Visualizer(const std::string & output_dir)
: m_output_dir(output_dir), m_image(), m_gnuplot()
{
    if (m_output_dir.empty()) m_output_dir.append(".");
}

/**
 * Destructor. Virtual functions is a portable way to avoid inlining.
 */
virtual ~Visualizer() {}

/**
 * Function generates and optionally writes the whole property field
 * into a file in binary, gray-scaled PGM format, where all the values
 * are adjusted to [0..255] interval.
 * N O T E: when image is being created, its content is implicitly transposed,
 *          i.e. abscissa in faster than ordinate as opposed to the situation
 *          where we address the matrix elements.
 * N O T E: it does not make sense to call this function inside pfor(..) loop.
 * \param  title        file title (extension will be "pgm").
 * \param  field        global grid of subdomains to be saved as an image.
 * \param  time_index   discrete timestamp.
 * \param  write_image  true, if saving into file was intended.
 * \return  reference to image that can be plotted, for example, by Gnuplot.
 */
virtual const ubyte_array_t & Write(const char * title,
            const Grid<Matrix,2> & field, int time_index, bool write_image)
{
	using ::allscale::api::user::algorithm::pfor;
    const bool flipY = false;

    point2d_t glo = GlobalSize(field);
    m_image.resize(glo.x * glo.y);

    // Compute the minimum and maximum values of the property field.
    double lower = 0.0, upper = 0.0;
    {
        Grid<double,2> minvalues(field.size()), maxvalues(field.size());
        pfor(point2d_t(0,0), field.size(), [&](const point2d_t & idx) {
            minvalues[idx] = *(std::min_element(field[idx].begin(),
                                                field[idx].end()));
            maxvalues[idx] = *(std::max_element(field[idx].begin(),
                                                field[idx].end()));
        });
        // Doing reduction sequentially.
        minvalues.forEach([&](const double & v) { if (lower > v) lower = v; });
        maxvalues.forEach([&](const double & v) { if (upper < v) upper = v; });
    }

    // Convert the property field into one-byte-per-pixel representation.
    pfor(point2d_t(0,0), field.size(), [&](const point2d_t & idx) {
        const auto & subfield = field[idx];
        for (int y = 0; y < SUBDOMAIN_Y; ++y) { int yg = Sub2GloY(idx, y);
        for (int x = 0; x < SUBDOMAIN_X; ++x) { int xg = Sub2GloX(idx, x);
            int pos = xg + glo.x * (flipY ? yg : (glo.y - 1 - yg));
            if (!(static_cast<size_t>(pos) < m_image.size()))
                assert_true(0);
            m_image[pos] = static_cast<unsigned char>(
                    Round(255.0 * (subfield(x,y) - lower) /
                            std::max(upper - lower, TINY)) );
        }}
    });

    // Write a file in binary, gray-scaled PGM format.
    if (write_image) {
        std::stringstream filename;
        filename << m_output_dir << PathSep << title;
        if (time_index >= 0) {
            filename << std::setfill('0') << std::setw(5) << time_index;
        }
        filename << ".pgm";
        std::fstream f(filename.str(), std::ios_base::out    |
                                       std::ios_base::binary |
                                       std::ios_base::trunc);
        assert_true(f.good()) << "failed to open file for writing: "
                              << filename.str() << std::endl;
        f << "P5\n" << glo.x << " " << glo.y << "\n" << int(255) << "\n";
        f.write(reinterpret_cast<const char*>(m_image.data()),
                static_cast<std::streamsize>(m_image.size()));
        f.flush();
    }
    return m_image;
}

/**
 * Function writes into file the image that shows the sensor locations
 * in the whole domain.
 * \param  H  grid of observation matrices.
 */
virtual void WriteImageOfSensors(const Grid<Matrix,2> & H)
{
    using ::allscale::api::user::algorithm::pfor;
    Grid<Matrix,2> field(H.size());
    pfor(point2d_t(0,0), H.size(), [&](const point2d_t & idx) {
        field[idx].Resize(SUBDOMAIN_X, SUBDOMAIN_Y);

        Matrix       & subfield = field[idx];
        const Matrix & h = H[idx];

        // Everything is set to black except border points painted in dark gray.
        Fill(subfield, 0.0);
        for (int x = 0; x < SUBDOMAIN_X; ++x) {
            subfield(x,0) = subfield(x,SUBDOMAIN_Y - 1) = 128.0;
        }
        for (int y = 0; y < SUBDOMAIN_Y; ++y) {
            subfield(0,y) = subfield(SUBDOMAIN_X - 1,y) = 128.0;
        }

        // Sensor locations are depicted in white.
        for (int r = 0; r < h.NRows(); ++r) {
        for (int c = 0; c < h.NCols(); ++c) {
            if (h(r,c) > 0.01) {                    // H is a 0/1 matrix
                div_t res = div(c, SUBDOMAIN_Y);    // 1D index to 2D point
                int x = res.quot;
                int y = res.rem;
                assert_true(Sub2Ind(x,y) == c);     // check index conversion
                subfield(x,y) = 255.0;
            }
        }}
    });

    const int Nx = field.size()[0] * SUBDOMAIN_X;
    const int Ny = field.size()[1] * SUBDOMAIN_Y;
    std::stringstream title;
    title << "sensors" << "_Nx" << Nx << "_Ny" << Ny;
    this->Write(title.str().c_str(), field, -1, true);
}

/**
 * Function generates and (optionally) plots image of a field, whose entries
 * are represented as the grayscaled values adjusted to [0..255] interval.
 * Optionally, the image can be written on disk in PGM format.
 * \param  field        2D field to be visualized as a grayscaled image.
 * \param  title        name of the field (will be used to create a file name).
 * \param  timestamp    current iteration index (discrete time).
 * \param  write_image  if true, then image will be saved on the dhard disk.
 * \param  plot_image   if true, image will be plotted in a graphical window.
 */
virtual void ShowImage(const Grid<Matrix,2> & field,
                       const char * title,
                       int timestamp,
                       bool write_image = true,
                       bool plot_image = true)
{
    if (!write_image && !plot_image) return;
    this->Write(title, field, timestamp, write_image);
    if (plot_image) {
        const size_t LEN = 63;
        char title[LEN + 1];
        snprintf(title, LEN, "frame %05d", timestamp);
        const point2d_t glo = GlobalSize(field);
        m_gnuplot.PlotGrayImage(m_image.data(), glo.x, glo.y, title, true);
    }
}

}; // Visualizer

#else   // !AMDADOS_DEBUGGING

/**
 * Stub class does nothing in production mode.
 */
class Visualizer
{
    ubyte_array_t empty;
public:
    explicit Visualizer(const std::string &) {}
    virtual ~Visualizer() {}
    virtual const ubyte_array_t & Write(const char *,
                    const Grid<Matrix,2> &, int, bool) { return empty; }
    virtual void WriteImageOfSensors(const Grid<Matrix,2> &) {}
    virtual void ShowImage(const Grid<Matrix,2> &, const char *,
                    int, bool w = true, bool p = true) { (void)w; (void)p; }
};

#endif  // AMDADOS_DEBUGGING

} // namespace amdados

