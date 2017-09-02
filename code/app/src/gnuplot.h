//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#undef AMDADOS_ENABLE_GNUPLOT
// Comment the line below to disable Gnuplot.
#define AMDADOS_ENABLE_GNUPLOT

#ifdef AMDADOS_ENABLE_GNUPLOT

#include <cstdlib>
#include <cstdio>
#include <unistd.h>            // for access(), mkstemp()

namespace amdados {
namespace app {
namespace gnuplot {
namespace details {
namespace {

// Linux or Mac OS.
#if defined(__APPLE__) || defined(__linux__) || defined(__unix__)

// On Mac OS ( https://apple.stackexchange.com/questions/103814/cant-plot-with-gnuplot-on-my-mac ):
// 1) Install XQuartz (via brew);
// 2) brew uninstall gnuplot
// 3) brew install gnuplot --with-x11
// Now gnuplot supports the x11 terminal.

const char   PathDelimiter[] = ":";
const char   GnuplotName[] = "gnuplot";
const char * GnuplotPath[] = {"/usr/local/bin/", "/usr/bin/"};
const char   StdTerminal[] = "x11";

// Function checks if file exists and has executable permission.
bool CheckGnuplotExists(const std::string & filename) {
    return (access(filename.c_str(), X_OK) == 0);
}

// System dependent way to open pipe.
FILE * PipeOpen(const char * name) {
    return popen(name, "w");
}

// System dependent way to close pipe.
bool PipeClose(FILE * p) {
    return !(pclose(p) == -1);
}

// Function checks existence of 'DISPLAY' environmental variable (and hence the display).
bool CheckDisplayExists() {
    if (getenv("DISPLAY") != NULL) return true;
    std::cout << "WARNING: 'DISPLAY' environment variable does not exist" << std::endl;
    return false;
}

// Function searches (1) user specified path; (2) standard folders;
// (3) PATH environmental variable for Gnuplot.
std::string GetProgramPath(const char * user_specified)
{
    std::string filename;

    // First look in standard folders. Start from user supplied path to Gnuplot.
    for (int i = -1; i < static_cast<int>(sizeof(GnuplotPath)/sizeof(GnuplotPath[0])); ++i) {
        if (i < 0) {                            // first check the user specified location
            if (user_specified != nullptr) {
                filename = std::string(user_specified) + "/" +  GnuplotName;
                if (CheckGnuplotExists(filename))
                    return filename;
            }
        } else {
            filename = std::string(GnuplotPath[i]) + "/" +  GnuplotName;
            if (CheckGnuplotExists(filename))
                return filename;
        }
    }
    filename.clear();

    // Second look in PATH environmental variable.
    char * path = getenv("PATH");
    if (path == nullptr) {
        std::cout << "PATH environment variable is not set" << std::endl;
        return std::string();
    }

    std::list<std::string> ls;
    std::string            pathStr = path;

    // Split the path string into a list of strings.
    for (size_t i = 0; i < pathStr.size();) {
        i = pathStr.find_first_not_of(" \t\r\n", i);         // skip leading white-spaces
        if (i == std::string::npos)
            break;                                           // nothing left
        size_t e = pathStr.find_first_of(PathDelimiter, i);  // end of the token
        if (e == std::string::npos) {
            ls.push_back(pathStr.substr(i));
            break;
        } else {
            ls.push_back(pathStr.substr(i, e - i));
        }
        i = e + 1;
    }

    // Scan list for Gnuplot program files.
    for (std::list<std::string>::const_iterator i = ls.begin(); i != ls.end(); ++i) {
        filename = (*i) + "/" +  GnuplotName;
        if (CheckGnuplotExists(filename))
            return filename;
    }

    std::cout << "WARNING: cannot find Gnuplot in PATH, standard or user-specified folder(s)"
              << std::endl;
    return std::string();
}

#endif // checking for the OS

} // anonymous namespace
} // namespace details

//=================================================================================================
/** C++ wrapper to Gnuplot command interface. */
//=================================================================================================
class Gnuplot
{
private:
    FILE *            mPipe;    ///< pointer to the stream that can be used to write to the pipe
    bool              mValid;   ///< validation of Gnuplot session
    std::vector<char> mBuffer;  ///< command buffer

private:
    Gnuplot(const Gnuplot &);
    Gnuplot & operator=(const Gnuplot &);

//-------------------------------------------------------------------------------------------------
/** Function closes connection to Gnuplot process, clears this object and returns false. */
//-------------------------------------------------------------------------------------------------
bool Clear()
{
    if (mPipe != nullptr) {
        if (details::PipeClose(mPipe)) {
            std::cout << "Gnuplot session was closed normally" << std::endl;
        } else {
            std::cout << "WARNING: failed to close Gnuplot session" << std::endl;
        }
    }
    mPipe = NULL;
    mValid = false;
    return false;
}

public:
//-------------------------------------------------------------------------------------------------
/** Constructor sets the optional alternative (non-standard) Gnuplot path.
 *  \note For windows: use path with slash '/' not backslash '\' .
 *  \param  path  user-specified path to Gnuplot.
 *  \param  args  Gnuplot command-line arguments. */
//-------------------------------------------------------------------------------------------------
explicit Gnuplot(const char * path = nullptr, const char * args = nullptr)
: mPipe(nullptr), mValid(false), mBuffer()
{
    std::string filename = details::GetProgramPath(path);
    Clear();
    if (!filename.empty() && details::CheckDisplayExists()) {
        if (args != nullptr) (filename += " ") += args;

        // Pipe is used to open the Gnuplot application and communicate with it.
        mPipe = details::PipeOpen(filename.c_str());
        if (mPipe) {
            mValid = true;
            SetStdTerminal();  // set terminal type
            return;
        }
    }
    std::cout << "WARNING: failed to open Gnuplot session, plotting is unavailable" << std::endl;
    Clear();
}

//-------------------------------------------------------------------------------------------------
/** Destructor closes connection to Gnuplot process. */
//-------------------------------------------------------------------------------------------------
virtual ~Gnuplot()
{
    Clear();
}

//-------------------------------------------------------------------------------------------------
/** Function sends a command to an active Gnuplot session. */
//-------------------------------------------------------------------------------------------------
Gnuplot & Command(const char * command)
{
    if (mValid && (command != nullptr)) {
        std::fprintf(mPipe, "%s\n", command);
        std::fflush(mPipe);
    }
    return (*this);
}

//-------------------------------------------------------------------------------------------------
// Function repeats the last 'plot' or 'splot' command. This can be useful for viewing a plot
// with different set options or when generating the same plot for several devices,
//-------------------------------------------------------------------------------------------------
Gnuplot & Replot()
{
    if (mValid) {
        std::fputs("replot\n", mPipe);
        std::fflush(mPipe);
    }
    return (*this);
}

//-------------------------------------------------------------------------------------------------
/** Function resets a Gnuplot session and sets all variables to default. */
//-------------------------------------------------------------------------------------------------
Gnuplot & ResetAll()
{
    if (mValid) {
        std::fputs("reset\nclear\n", mPipe);
        std::fflush(mPipe);
        SetStdTerminal();
    }
    return (*this);
}

//-------------------------------------------------------------------------------------------------
/** Function sets the standard (screen) terminal. */
//-------------------------------------------------------------------------------------------------
Gnuplot & SetStdTerminal()
{
    if (mValid) {
        std::fprintf(mPipe, "set output\nset terminal %s\n", details::StdTerminal);
        std::fflush(mPipe);
    }
    return (*this);
}

//-------------------------------------------------------------------------------------------------
/** Function saves a Gnuplot session to a postscript file, filename without extension.
 *  Default name is "gnuplot_output". */
//-------------------------------------------------------------------------------------------------
Gnuplot & SetPostscriptTerminal(const char * filename)
{
    if (mValid) {
        if (filename == nullptr) filename = "gnuplot_output";
        fprintf(mPipe, "set terminal postscript color\nset output \"%s.ps\"\n", filename);
        fflush(mPipe);
    }
    return (*this);
}

//-------------------------------------------------------------------------------------------------
// Function plots greyscaled image (Gnuplot 4.2+).
// It assumed that a pixel value at (x,y) is addressed as follows: image[x*height + y],
// i.e., ordinate changes faster then abscissa. The image is implicitly transposed upon
// drawing, so that the "image(x,y)" looks as expected.
//-------------------------------------------------------------------------------------------------
Gnuplot & PlotGrayImage(const unsigned char * image, int width, int height,
                        const std::string & title, bool flipY = false)
{
    if (!mValid) return (*this);
    assert_true((image != nullptr) && (width > 0) && (height > 0));

    const size_t expected_size = 1024 + 4*width*height; // header size + 4 bytes per pixel
    if (mBuffer.size() < expected_size) mBuffer.resize(expected_size, '\0');

    while (1) {
        long   N = static_cast<long>(mBuffer.size());
        char * p = mBuffer.data();
        std::fill(mBuffer.begin(), mBuffer.end(), 0);   // not actually necessary but just in case

        // Write the header.
        long n = snprintf(p, static_cast<size_t>(N),
            "\t\t\t\t\t\t\t"                        // sometimes pipe "swallows" the first symbol
            /*"set terminal x11\n"*/
            "unset key\n"
            "set title \"%s\"\n"
            "set xrange [0:%d] noreverse nowriteback\n"
            "set yrange [0:%d] noreverse nowriteback\n"
            "set palette gray\n"
            "unset colorbox\n"
            "set tics out\n"
            "set autoscale noextend\n"      // keepfix
            "unset logscale\n"
            "plot '-' matrix with image pixels\n", title.c_str(), (width-1), (height-1));
        // Return upon error, otherwise proceed if buffer is not overwhelmed.
        if (n < 0) {
            std::cout << "WARNING: snprintf() failed to create a Gnuplot command" << std::endl;
            return (*this);
        } else if (n < N) {
            p += n;             // advance pointer right beyond the header
            N -= n;             // reduce effective size

            // Copy image into string buffer in textual format. I use 'n+9 < N' guard to exclude
            // any possibility to overrun the buffer (actually 'n+5' would be enough).
            n = 0;
            for (int y = 0; y < height; ++y) {
                const int yy = flipY ? (height-1-y) : y;
                for (int x = 0; x < width; ++x) {
                    if (n+9 < N) {  // convert unsigned char into up to 3-char string
                        div_t res = div(static_cast<int>(image[x + width * yy]), 100);
                        if (res.quot > 0) {             // pixel value >= 100
                            p[n++] = '0' + res.quot;
                            res = div(res.rem, 10);
                            p[n++] = '0' + res.quot;
                        } else {                        // pixel value < 100
                            res = div(res.rem, 10);
                            if (res.quot > 0) p[n++] = '0' + res.quot;
                        }
                        p[n++] = '0' + res.rem;
                        p[n++] = ' ';
                    } else {
                        n = N;      // signals insufficient buffer size
                    }
                }
                p[n++] = '\n';
            }
            // Enclosing footer.
            if (n+9 < N) {
                p[n++] = 'e';
                p[n++] = '\n';
                p[n++] = 'e';
                p[n++] = '\n';
                p[n++] = '\0';
            } else {
                n = N;              // signals insufficient buffer size
            }
        }

        // FILE * f = fopen("cmd.txt","wt"); fprintf(f,"%s",mBuffer.data()); fclose(f); exit(0);

        // Break if command and image have been written normally.
        if (n < N) break;
        // Otherwise reallocate buffer and repeat.
        mBuffer.resize(std::max(2*mBuffer.size(), size_t(1024)), '\0');
    }

    // Send the command along with image.
    Command(mBuffer.data());
    return (*this);
}

}; // class Gnuplot

} // namespace gnuplot
} // namespace app
} // namespace amdados

#else   // !AMDADOS_ENABLE_GNUPLOT

namespace amdados {
namespace app {
namespace gnuplot {

//=================================================================================================
/** Stub class does nothing. */
//=================================================================================================
class Gnuplot
{
public:
    explicit Gnuplot(const char * a = nullptr, const char * b = nullptr) { (void)a; (void)b; }
    virtual ~Gnuplot() {}
    Gnuplot & Command(const char *) { return *this; }
    Gnuplot & Replot() { return *this; }
    Gnuplot & ResetAll() { return *this; }
    Gnuplot & SetStdTerminal() { return *this; }
    Gnuplot & SetPostscriptTerminal(const char *) { return *this; }
    Gnuplot & PlotGrayImage(const unsigned char *, int, int,
                        const std::string &, bool flipY = false) { (void)flipY; return *this; }
};

} // namespace gnuplot
} // namespace app
} // namespace amdados

#endif  // AMDADOS_ENABLE_GNUPLOT

