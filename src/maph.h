
// TODO:
//
//     - More JSON inputs
//         * subsampling step size
//         * background color:  0, NaN, or other
//     - Parse data from other apps -- Apple Activity?
//     - Look into Strava API to automatically pull updated activities
//     - Upload build artifacts?
//     - Documentation
//     - Error checking
//     - Benchmark
//     - Refactor
//         * try/throw/catch?
//

//======================================================================

// Standard

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <sstream>

#if !(defined(_WIN32) || defined(_WIN64))
	#include <glob.h>
#endif

// mkdir
#if defined(_WIN32) || defined(_WIN64)
	#include <direct.h>
#else
	#include <sys/stat.h>
#endif

//======================================================================

// Third party

#include <json.hpp>
using json = nlohmann::json;

#include <pugixml.hpp>

//======================================================================

// Custom

#include <colormapper.h>
#include <constants.h>

//======================================================================

const std::string me = "maph";

enum Kernel {cylinder, pyramid, cone, gaussian};

// explicitly defined for backwards-compatability
enum Sampling {pointwise = 0, autosample = 1, linear = 2};

// Also canoe, e-bike ride, etc. but those are not supported here
// (mainly because I've never done them personally!)
enum Type {hike, nordicSki, ride, run, walk};

// TODO:
// - speed
// - heartrate
// - distance (or time, but not like strava #TimeMap)
enum MapType {elevation, gradient, heat};

class Settings
{
	public:

		int nx, ny, nxy;

		bool fit, fitx, fity, fitnx, fitny;
		double minx, maxx, miny, maxy, mins, maxs, margin;
		double degPerPix;

		// Kernel radius (pixels)
		int r;

		// Kernel width
		int tr1;

		// Coefficient (spikiness) of Gaussian kernel.
		double gka;

		// Variance of Gaussian kernel.  Note r, gka, and gkb are
		// not independent.
		double gkb;

		Sampling sampling;

		Kernel kernel;

		// Linear sampling stepsize in units of kernel radii.  TODO:
		// optionally load from JSON
		double step = 0.25;

		std::string fname;
		std::string fgpx;
		int verb;

		bool stat;

		irwincolor::Map c;
		bool allCmaps;
		std::string cmapFile;
		std::vector<std::string> mapNames;

		// Filter activity type?
		bool filterType;
		Type type;

		MapType mapType;

		// Filter by start/end date?
		bool filterAfter, filterBefore;
		std::tm aftertm, beforetm;
		time_t after, before;

		// Bytes per pixel:  RGBA
		const int bpp = 4;

		const std::string dlm = "_-_", imgext = ".png";
};

class Data
{
	public:
		std::vector<std::string> gpxs;      // filenames
		std::vector<double> lats;           // lattitudes
		std::vector<double> lons;           // longitudes
		std::vector<double> eles;           // elevations
		std::vector<double> scas;           // scalars for non-heat map type
		std::vector<unsigned int> iEndSeg;  // index of each segment's end
		int ntrkptsum;                      // number of trackpoints
};

class Transformation
{
	public:

		bool enable;
		double x, y, rot;
		std::vector<double> mat;
		std::string outdir;
};

template<class T> std::ostream& operator<<(std::ostream& stream, const std::vector<T>& values)
{
	// cout (or other stream) << a vector
	stream << "{";
	copy(begin(values), end(values), std::ostream_iterator<T>(stream, ", "));
	stream << '}';
	return stream;
}

std::string getKernelName(Kernel k)
{
	// If you like it then you shoulda put a string on it
	if (k == gaussian)
		return "gaussian";
	else if (k == cylinder)
		return "cylinder";
	else if (k == pyramid)
		return "pyramid";
	else
		return "cone";
}

std::string getSamplingName(Sampling s)
{
	if (s == pointwise)
		return "pointwise";
	else if (s == linear)
		return "linear";
	else
		return "autosample";
}

std::string getMapTypeName(MapType t)
{
	if (t == elevation)
		return "elevation";
	else if (t == gradient)
		return "gradient";
	else
		return "heat";
}

std::string getTypeName(Type t)
{
	if (t == hike)
		return "Hike";
	if (t == nordicSki)
		return "Nordic Ski";
	else if (t == ride)
		return "Ride";
	else if (t == walk)
		return "Walk";
	else
		return "Run";
}

int getTypeInt(Type t)
{
	// Is this a strava convention?  Or GPX standard?
	if (t == hike)
		return 4;
	if (t == nordicSki)
		return 7;
	else if (t == ride)
		return 1;
	else if (t == walk)
		return 10;
	else // run
		return 9;
}

void getFiles(std::string &pattern, std::vector<std::string> &fileList)
{
	//std::cout << "starting getFiles..." << std::endl;
#if defined(_WIN32) || defined(_WIN64)

	std::replace(pattern.begin(), pattern.end(), '/', '\\');

	// Wildcard globbing only works on Windows after the last
	// directory separator and only with a single wildcard.
	std::string dir = pattern.substr(0, pattern.find_last_of("\\"));

	// There are some fancy Windows API methods for calling dir, but
	// they don't work with wildcards!
	std::string ftmp = std::tmpnam(nullptr);

	// In cmd, ftmp is a full path, e.g.
	// "C:\Users\jirwi\AppData\Local\Temp\s26c.0".  In git bash, it
	// is just "\sgjg.", so we need to prepend a leading ".".
	if (ftmp.substr(0, 1) == "\\") ftmp = "." + ftmp;
	//std::cout << "ftmp = " << ftmp << std::endl;

	std::string cmd = (std::string) "dir /b " + pattern + " > " + ftmp;
	system(cmd.c_str());

	std::ifstream ifs(ftmp);
	std::string line;
	while (std::getline(ifs, line))
	{
		// Dir lists the filename but not its path, so we need to
		// add it back.
		std::string file = dir + slash + line;
		//std::replace(file.begin(), file.end(), '\\', '/');
		//std::cout << "file = " << file << "\n";
		fileList.push_back(file);
	}

	ifs.close();
	cmd = (std::string) "del " + ftmp;
	system(cmd.c_str());

#else

	// Declare glob_t for storing the results of globbing
	glob_t globbuf;

	// Glob. GLOB_TILDE tells the globber to expand "~" in the pattern to the home directory
	glob(pattern.c_str(), GLOB_TILDE, NULL, &globbuf);

	for (int i = 0; i < globbuf.gl_pathc; ++i)
		fileList.push_back(globbuf.gl_pathv[i]);

	// Free the globbuf structure
	if (globbuf.gl_pathc > 0)
		globfree(&globbuf);

#endif
	//std::cout << "done getFiles..." << std::endl;
}

void printBounds(const Settings& s)
{
	std::cout << "Bounds:\n\t"
			<< s.minx << " <= x <= " << s.maxx
			<< "\n\t"
			<< s.miny << " <= y <= " << s.maxy
			<< std::endl;
}

std::vector<uint8_t> colorPixels(const Settings& s,
		const std::vector<unsigned int>& img,
		const std::vector<unsigned int>& idx)
{
	std::vector<uint8_t> pix(s.bpp * s.nxy);
	double x, x0;
	x0 = 0.0;

	// Number of non-zeros
	unsigned int nnonzero = 0;
	for (int i = 0; i < s.nxy; i++)
	{
		if (img[i] != 0) nnonzero++;
	}

	int j = 0, n;
	for (int i = 0; i < s.nxy; i++)
	{
		// Map img[idx[i]] to x in range [0, 1].
		if (img[idx[i]] == 0)
		{
			if (s.mapType != heat)
			{
				// Map 0 to NaN color.  Otherwise evenly distribute.
				x = 2.0;
			}
			else
			{
				// Leave 0 as is.
				x = 0.0;
			}
		}
		else if (i > 0 && img[idx[i]] == img[idx[i-1]])
		{
			j++;
			x = x0;
		}
		else
		{
			j++;
			x = (double) j / nnonzero;
		}
		x0 = x;

		// Map x to RGB triple.
		std::vector<uint8_t> rgb = s.c.map(x);

		//std::cout << "rgb = " << rgb[0] << " " << rgb[1] << " " << rgb[2] << "\n";

		// Save RGB triple to pixel.
		int ib = idx[i] * s.bpp;
		pix[ib + 0] = rgb[0];
		pix[ib + 1] = rgb[1];
		pix[ib + 2] = rgb[2];
		pix[ib + 3] = irwincolor::twofivefive;  // alpha
	}
	return pix;
}

void addKernel(const Settings& s, double lat, double lon,
		std::vector<unsigned int>& img, std::vector<unsigned int>& kernel,
		double scalar)
{
	// Add a single heat kernel to the image at coordinate (lat, lon).

	int ix0 = floor(s.nx * (lon - s.minx) / (s.maxx - s.minx));
	if (ix0 < -s.r || ix0 > s.nx + s.r) return;
	int iy0 = floor(s.ny * (lat - s.miny) / (s.maxy - s.miny));
	if (iy0 < -s.r || iy0 > s.ny + s.r) return;

	for (int dy = -s.r; dy <= s.r; dy++)
	{
		int iy = iy0 + dy;
		if (0 <= iy && iy < s.ny)
		{
			int dyr = dy + s.r;
			int iyk = s.tr1 * (s.tr1 - dyr - 1);
			int iyp = s.nx  * (s.ny  - iy  - 1);
			for (int dx = -s.r; dx <= s.r; dx++)
			{
				int ix = ix0 + dx;
				if (0 <= ix && ix < s.nx)
				{
					int dxr = dx + s.r;
					int ik = iyk + dxr;
					int ip = iyp + ix ;

					if (s.mapType != heat)
					{
						// Just set, don't add.  This should
						// probably be the default for non-heat
						// maps.  Could use max() instead for
						// non-cylinder kernels.
						if (kernel[ik] > 0)
							img[ip] = floor(kernel[ik] * scalar) + 1;
					}
					else // heat
					{
						img[ip] += kernel[ik];
					}

				}
			}
		}
	}
}

int usage()
{
	std::cout << "Usage:" << std::endl;
	std::cout << "\t" << me << " input.json [-t dir xshift yshift rot]" << std::endl;
	return ERR_CMD_ARG;
}

int loadArgs(int argc, char* argv[], Settings& s, Transformation& t, std::string& fjson)
{
	t.enable = false;
	s.stat = false;

	if (argc < 2)
	{
		return usage();
	}

	int i = 0;
	while (i < argc - 1)
	{
		try
		{
			std::string str = argv[++i];
			if (str == "-t")
			{
				// Transformation args
				t.enable = true;
				t.outdir = argv[++i];
				std::cout << "Transformation directory = \"" << t.outdir << "\".\n";
				t.x   = std::stod(argv[++i]);
				t.y   = std::stod(argv[++i]);
				t.rot = std::stod(argv[++i]);
				std::cout << "X translation = " << t.x   << "\n";
				std::cout << "Y translation = " << t.y   << "\n";
				std::cout << "Rotation      = " << t.rot << "\n";
				std::cout << std::endl;
			}
			else if (str == "-s")
			{
				// Summary statistics
				s.stat = true;
			}
			else
				fjson = str;
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << std::endl;
			return usage();
		}
	}
	return 0;
}

int loadJson(std::string& fjson, json& inj)
{
	std::cout << "Loading JSON input file \"" << fjson << "\" ..." << std::endl;
	std::ifstream ifs(fjson, std::ifstream::in);

	try
	{
		ifs >> inj;
		ifs.close();  // close to release lock on file
	}
	catch (const std::exception& e)
	{
		std::cout << "\nError:  cannot load JSON input file \"" << fjson << "\"." << std::endl;
		std::cout << e.what() << std::endl;
		return irwincolor::ERR_JSON;
	}

	return 0;
}


std::tm str2tm(std::string str)
{
	// Convert string str to time tm

	// TODO:  exceptions

	//std::cout << "str = " << str << std::endl;

	std::tm tm = {};
	std::istringstream ss(str);

	// ISO 8601:  YYYY-MM-DD"T"hh:mm::ss"Z"
	ss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%S");
	//std::cout << std::put_time(&tm, "%c") << "\n\n";

	return tm;
}

int loadSettings(Settings& s, json& inj, std::string& fjson)
{
	// Initial defaults
	s.nx = 1280;
	s.ny = 720;
	s.verb = 0;
	s.sampling = autosample;
	s.r = 10;
	s.minx = 0.0;
	s.miny = 0.0;
	s.maxx = 0.0;
	s.maxy = 0.0;
	s.margin = 0.0;
	s.fitx = false;
	s.fity = false;
	s.fitnx = false;
	s.fitny = false;
	s.fname = fjson.substr(0, fjson.find_last_of("."));
	s.fgpx = "*.gpx";
	s.cmapFile = "";
	std::string cmapName = "";
	s.allCmaps = false;
	s.gka = 10.0;
	s.kernel = cone;
	s.filterType = false;
	s.filterAfter = false;
	s.after = 0;
	s.filterBefore = false;
	s.before = 0;
	s.mapType = heat;

	bool bminx = false, bminy = false, bmaxx = false, bmaxy = false;

	// JSON keys
	const std::string sizexId = "Image size x";
	const std::string sizeyId = "Image size y";
	const std::string verbId = "Verbosity";
	const std::string sampleId = "Sampling";
	const std::string radiusId = "Kernel radius";
	const std::string minxId = "Min x";
	const std::string minyId = "Min y";
	const std::string maxxId = "Max x";
	const std::string maxyId = "Max y";
	const std::string marginId = "Margin";
	const std::string fitxId = "Fit x";
	const std::string fityId = "Fit y";
	const std::string fitnxId = "Fit nx";
	const std::string fitnyId = "Fit ny";
	const std::string fnameId = "File name prefix";
	const std::string fgpxId = "GPX files";
	const std::string cmapFileId = "Colormap file";
	const std::string cmapNameId = "Colormap name";
	const std::string allCmapsId = "All colormaps";
	const std::string imapId = "Color map";
	const std::string invertMapId = "Invert color map";
	const std::string gaussianAmpId = "Gaussian amplitude";
	const std::string kernelId = "Kernel";
	const std::string typeId = "Type";
	const std::string afterId = "After";
	const std::string beforeId = "Before";
	const std::string mapTypeId = "Map type";

	for (json::iterator it = inj.begin(); it != inj.end(); it++)
	{
		//std::cout << "\nkey = " << it.key();
		if (it.key() == sizexId && it.value().is_number())
		{
			s.nx = it.value();
		}
		else if (it.key() == sizeyId && it.value().is_number())
		{
			s.ny = it.value();
		}
		else if (it.key() == verbId && it.value().is_number())
		{
			s.verb = it.value();
		}
		else if (it.key() == sampleId
				&& (it.value().is_number() || it.value().is_string()))
		{
			if (it.value().is_number())
				s.sampling = it.value();
			else
			{
				if (it.value() == getSamplingName(pointwise))
					s.sampling = pointwise;
				else if (it.value() == getSamplingName(linear))
					s.sampling = linear;
				else if (it.value() == getSamplingName(autosample))
					s.sampling = autosample;
				else
					std::cout << "\nWarning:  unknown sampling value "
							<< it.value() << ", defaulting to \""
							<< getSamplingName(autosample) << "\"\n";
			}
		}
		else if (it.key() == radiusId && it.value().is_number())
		{
			s.r = it.value();
		}
		else if (it.key() == minxId && it.value().is_number())
		{
			s.minx = it.value();
			bminx = true;
		}
		else if (it.key() == minyId && it.value().is_number())
		{
			s.miny = it.value();
			bminy = true;
		}
		else if (it.key() == maxxId && it.value().is_number())
		{
			s.maxx = it.value();
			bmaxx = true;
		}
		else if (it.key() == maxyId && it.value().is_number())
		{
			s.maxy = it.value();
			bmaxy = true;
		}
		else if (it.key() == marginId && it.value().is_number())
		{
			s.margin = it.value();
		}
		else if (it.key() == fitxId && it.value().is_boolean())
		{
			s.fitx = it.value();
		}
		else if (it.key() == fityId && it.value().is_boolean())
		{
			s.fity = it.value();
		}
		else if (it.key() == fitnxId && it.value().is_boolean())
		{
			s.fitnx = it.value();
		}
		else if (it.key() == fitnyId && it.value().is_boolean())
		{
			s.fitny = it.value();
		}
		else if (it.key() == fnameId && it.value().is_string())
		{
			s.fname = it.value();
		}
		else if (it.key() == fgpxId && it.value().is_string())
		{
			s.fgpx = it.value();
		}
		else if (it.key() == cmapFileId && it.value().is_string())
		{
			s.cmapFile = it.value();
		}
		else if (it.key() == cmapNameId && it.value().is_string())
		{
			cmapName = it.value();
		}
		else if (it.key() == imapId && it.value().is_number())
		{
			s.c.imap = it.value();
		}
		else if (it.key() == invertMapId && it.value().is_boolean())
		{
			s.c.inv = it.value();
		}
		else if (it.key() == allCmapsId && it.value().is_boolean())
		{
			s.allCmaps = it.value();
		}
		else if (it.key() == gaussianAmpId && it.value().is_number())
		{
			s.gka = it.value();
		}
		else if (it.key() == kernelId && it.value().is_string())
		{
			if (it.value() == getKernelName(pyramid))
				s.kernel = pyramid;
			else if (it.value() == getKernelName(cylinder))
				s.kernel = cylinder;
			else if (it.value() == getKernelName(gaussian))
				s.kernel = gaussian;
			else if (it.value() == getKernelName(cone))
				s.kernel = cone;
			else
				std::cout << "\nWarning:  unknown kernel value "
						<< it.value() << ", defaulting to \""
						<< getKernelName(cone) << "\"\n";
		}
		else if (it.key() == typeId && it.value().is_string())
		{
			s.filterType = true;
			if (it.value() == getTypeName(hike))
				s.type = hike;
			else if (it.value() == getTypeName(nordicSki))
				s.type = nordicSki;
			else if (it.value() == getTypeName(ride))
				s.type = ride;
			else if (it.value() == getTypeName(run))
				s.type = run;
			else if (it.value() == getTypeName(walk))
				s.type = walk;
			else
			{
				s.filterType = false;
				std::cout << "\nWarning:  unknown type value "
						<< it.value() << ", defaulting to all activity types.\n";
			}
		}
		else if (it.key() == mapTypeId && it.value().is_string())
		{
			if (it.value() == getMapTypeName(elevation))
				s.mapType = elevation;
			else if (it.value() == getMapTypeName(gradient))
				s.mapType = gradient;
			else if (it.value() == getMapTypeName(heat))
				s.mapType = heat;
			else
			{
				std::cout << "\nWarning:  unknown map type value "
						<< it.value() << ", defaulting to heat map.\n";
			}
		}
		else if (it.key() == afterId && it.value().is_string())
		{
			s.filterAfter = true;
			s.aftertm = str2tm((std::string) it.value());
			s.after = mktime(&s.aftertm);
		}
		else if (it.key() == beforeId && it.value().is_string())
		{
			s.filterBefore = true;
			s.beforetm = str2tm((std::string) it.value());
			s.before = mktime(&s.beforetm);
		}
		else
		{
			std::cout << "\nWarning:  unknown JSON key \"" << it.key()
					<< "\" with value \"" << it.value() << "\"\n";
		}
	}

	s.fit = !(bminx && bminy && bmaxx && bmaxy);

	if (s.mapType != heat)
	{
		// Linear sampling, cylinder kernel only for non-heat maps
		s.sampling = linear;
		s.kernel = cylinder;
	}

	// Echo inputs
	std::cout << "\n";
	std::cout << "Image size" << " = " << s.nx << ", " << s.ny << "\n";
	std::cout << fnameId << " = \"" << s.fname << "\"\n";
	std::cout << fgpxId << " = \"" << s.fgpx << "\"\n";
	std::cout << cmapFileId << " = \"" << s.cmapFile << "\"\n";
	std::cout << cmapNameId << " = \"" << cmapName << "\"\n";
	std::cout << allCmapsId << " = \"" << s.allCmaps << "\"\n";
	std::cout << imapId << " = " << s.c.imap << "\n";
	std::cout << invertMapId << " = " << s.c.inv << "\n";
	std::cout << verbId << " = " << s.verb << "\n";
	std::cout << sampleId << " = " << getSamplingName(s.sampling) << "\n";
	std::cout << radiusId << " = " << s.r << "\n";
	std::cout << gaussianAmpId << " = " << s.gka << "\n";
	std::cout << kernelId << " = " << getKernelName(s.kernel) << "\n";
	std::cout << marginId << " = " << s.margin << "\n";
	std::cout << fitxId << " = " << s.fitx << "\n";
	std::cout << fityId << " = " << s.fity << "\n";
	std::cout << fitnxId << " = " << s.fitnx << "\n";
	std::cout << fitnyId << " = " << s.fitny << "\n";
	std::cout << mapTypeId << " = " << getMapTypeName(s.mapType) << "\n";

	if (s.filterType) std::cout << typeId << " = " << getTypeName(s.type) << "\n";

	if (s.filterAfter) std::cout << afterId << " = " << std::put_time(&s.aftertm, "%c") << "\n";
	if (s.filterBefore) std::cout << beforeId << " = " << std::put_time(&s.beforetm, "%c") << "\n";

	//std::cout << "s.fit = " << s.fit << "\n";
	if (!s.fit)
		printBounds(s);

	std::cout << std::endl;

	if (s.cmapFile != "" && (cmapName != "" || s.allCmaps))
	{
		int io;
		if (s.allCmaps)
		{
			io = irwincolor::loadMapNames(s.cmapFile, s.mapNames);
			//std::cout << "mapNames = " << s.mapNames << std::endl;
		}
		else
			io = s.c.load(s.cmapFile, cmapName);
		if (io != 0) return io;
	}
	return 0;
}

double distance(double lat1i, double lon1i, double lat2i, double lon2i,
		double ele1i, double ele2i)
{
	// Use the haversine formula to get the distance in miles
	// between two lat/lon (degrees) / elevation (meters) points on
	// the Earth

	//double r = 3958.7613;   // mean radius (mi)
	double r = 3949.9028;   // polar radius (mi)
	double mipm = 1609.34;  // 1 mile in meters (misnomer)

	// Fudge factors are for matching 2939617387.gpx to strava.
	// It's not exact.  Strava may apply smoothing to GPS points
	// before summing distances.

	double fudge = 1.0;
	//double fudge = 0.9910412; // haversine w/o ele
	//double fudge = 0.9902064;  // haversine w/ ele
	//double fudge = 0.9873515;  // WGS-84 w/ ele
	//double fudge = 0.9881791;  // WGS-84 w/o ele.  Closest so far.

	// WGS-84
	double a = 6378137.0   ;  // equatorial radius     (m)
	double b = 6356752.3142;  // polar semi-minor axis (m)

	//std::cout << "pi = " << pi << "\n";
	//std::cout << "lats/lons = " << lat1i << " " << lon1i
	//                     << " " << lat2i << " " << lon2i
	//                     << " " << ele1i << " " << ele2i << "\n";

	// Radians
	double lat1 = lat1i * pi / 180.0;
	double lon1 = lon1i * pi / 180.0;
	double lat2 = lat2i * pi / 180.0;
	double lon2 = lon2i * pi / 180.0;

	//// WGS-84

	// TODO:  function
	double n1 = pow(a, 2) / sqrt(pow(a, 2) * pow(cos(lat1), 2) + pow(b, 2) * pow(sin(lat1), 2));
	double x1 = (              n1 + ele1i) * cos(lat1) * cos(lon1);
	double y1 = (              n1 + ele1i) * cos(lat1) * sin(lon1);
	double z1 = (pow(b/a, 2) * n1 + ele1i) * sin(lat1);

	double n2 = pow(a, 2) / sqrt(pow(a, 2) * pow(cos(lat2), 2) + pow(b, 2) * pow(sin(lat2), 2));
	double x2 = (              n2 + ele2i) * cos(lat2) * cos(lon2);
	double y2 = (              n2 + ele2i) * cos(lat2) * sin(lon2);
	double z2 = (pow(b/a, 2) * n2 + ele2i) * sin(lat2);

	double d = sqrt(pow(x2-x1, 2) + pow(y2-y1, 2) + pow(z2-z1, 2));
	//d = sqrt(pow(d, 2) - pow(ele2i - ele1i, 2));
	d /= mipm;

	//if (d > 1.0)
	//{
	//	std::cout << "d = " << d << "\n";
	//	std::cout << "lats/lons = " << lat1i << " " << lon1i
	//	                     << " " << lat2i << " " << lon2i << "\n";
	//}

	d *= fudge;

	return d;
}

int loadGpxs(Settings&s, Transformation& t, Data& d)
{
	d.ntrkptsum = 0;
	double lat, lon, ele, ele0;
	for (int ig = 0; ig < d.gpxs.size(); ig++)
	{
		int ntrkpt = 0;
		double dist = 0.0;
		bool errEle = false;

		if (s.verb > 0 || s.stat) std::cout << "GPX file = " << d.gpxs[ig] << std::endl;

		pugi::xml_document doc;
		if (int io = irwincolor::loadXml(d.gpxs[ig], doc) != 0)
			return io;

		std::string xquery;
		pugi::xpath_node trkpt0, trk;
		try
		{
			if (s.filterType)
			{
				xquery = "/gpx/trk";
				trk = doc.select_node(xquery.c_str());
				if (!trk) throw irwincolor::xPathException;

				int typel = std::stoi(trk.node().child_value("type"));
				//std::cout << "typel = " << typel << "\n";

				if (typel != getTypeInt(s.type)) continue;
			}

			if (s.filterAfter || s.filterBefore)
			{
				xquery = "/gpx/metadata";
				pugi::xpath_node meta = doc.select_node(xquery.c_str());
				if (!meta) throw irwincolor::xPathException;

				std::string time = meta.node().child_value("time");
				//std::cout << "time = " << time << "\n";

				std::tm t = str2tm(time);
				time_t tt = mktime(&t);

				if (s.filterBefore && s.before < tt     ) continue;
				if (s.filterAfter  &&       tt < s.after) continue;

			}

			xquery = "/gpx/trk/trkseg/trkpt";
			trkpt0 = doc.select_node(xquery.c_str());
			if (!trkpt0) throw irwincolor::xPathException;

			for (pugi::xml_node trkpt = trkpt0.node(); trkpt;
					trkpt = trkpt.next_sibling("trkpt"))
			{
				ntrkpt++;
				lat = std::stod(trkpt.attribute("lat").value());
				lon = std::stod(trkpt.attribute("lon").value());
				//std::cout << "lat, lon = " << lat << ", " << lon << "\n";

				ele0 = ele;
				try
				{
					ele = std::stod(trkpt.child_value("ele"));
					//std::cout << "ele = " << ele << "\n";
				}
				catch (const std::exception& e)
				{
					ele = 0.0;
					errEle = true;
				}

				d.lats.push_back(lat);
				d.lons.push_back(lon);
				d.eles.push_back(ele);

				if (s.stat)
				{
					if (ntrkpt > 1)
					{
						// TODO:  this doesn't match strava exactly.
						// Need to account for elevation (careful,
						// GPX is in meters) and maybe a better (not
						// just spherical) geophysical model.

						dist += distance(d.lats[d.ntrkptsum + ntrkpt - 2],
						                 d.lons[d.ntrkptsum + ntrkpt - 2],
						                 d.lats[d.ntrkptsum + ntrkpt - 1],
						                 d.lons[d.ntrkptsum + ntrkpt - 1],
						                 ele0, ele);
					}
				}

			}  // trackpoint loop

			// Save scalars
			if (s.mapType == elevation)
			{
				for (int i = 0; i < ntrkpt; i++)
					d.scas.push_back(d.eles[d.ntrkptsum + i]);
			}
			else if (s.mapType == gradient)
			{
				// TODO:  gradient units?  Not that it matters at
				// all for the colormap, but I think these are
				// meters per mile.

				// First point, a little WET
				double dist = distance(
					d.lats[d.ntrkptsum],
					d.lons[d.ntrkptsum],
					d.lats[d.ntrkptsum + 1],
					d.lons[d.ntrkptsum + 1],
					d.eles[d.ntrkptsum],
					d.eles[d.ntrkptsum + 1]);

				double s0;
				double s = (d.eles[1] - d.eles[0]) / dist;
				d.scas.push_back(s);

				for (int i = 1; i < ntrkpt - 1; i++)
				{
					s0 = s;
					dist = distance(
						d.lats[d.ntrkptsum + i],
						d.lons[d.ntrkptsum + i],
						d.lats[d.ntrkptsum + i + 1],
						d.lons[d.ntrkptsum + i + 1],
						d.eles[d.ntrkptsum + i],
						d.eles[d.ntrkptsum + i + 1]);

					s = (d.eles[d.ntrkptsum + i + 1] - d.eles[d.ntrkptsum + i]) / dist;

					// Average of preceeding and following line segments
					d.scas.push_back(0.5 * (s0 + s));

				}

				// Final point
				d.scas.push_back(s);
			}

			if (s.stat) std::cout << "Distance = " << dist << " mi\n";
			if (s.verb > 0) std::cout << "\tNumber of track points = " << ntrkpt << std::endl;

			if (t.enable)
			{
				std::string fogpx = t.outdir + slash
						+ d.gpxs[ig].substr(d.gpxs[ig].find_last_of(slash) + 1);

				if (s.verb > 0) std::cout << "fogpx = " << fogpx << std::endl;

				std::ofstream ofgpx(fogpx);
				ofgpx << std::setprecision(16);
				ofgpx << "<gpx>\n";
				ofgpx << "\t<trk>\n";
				ofgpx << "\t\t<trkseg>\n";

				// Shuffle.  Could initialize with iota() here.
				std::vector<unsigned int> id(ntrkpt);
				for (unsigned int i = 0; i < ntrkpt; i++)
					id[i] = i;
				std::random_shuffle(id.begin(), id.end());
				//std::cout << "id = " << id << std::endl;

				for (int it = 0; it < ntrkpt; it++)
				{
					unsigned int j = id[it];

					// Translate
					d.lons[j] += t.x;
					d.lats[j] += t.y;

					// Rotate
					lon = t.mat[0] * d.lons[j] + t.mat[1] * d.lats[j];
					lat = t.mat[2] * d.lons[j] + t.mat[3] * d.lats[j];

					ofgpx << "\t\t\t<trkpt lat=\"" << lat
					                << "\" lon=\"" << lon << "\"/>\n";
				}

				d.lats.clear();
				d.lats.shrink_to_fit();
				d.lons.clear();
				d.lons.shrink_to_fit();

				ofgpx << "\t\t</trkseg>\n";
				ofgpx << "\t</trk>\n";
				ofgpx << "</gpx>\n";

				ofgpx.close();
			}

			d.ntrkptsum += ntrkpt;
			d.iEndSeg.push_back(d.ntrkptsum - 1);

			if (errEle)
			{
				std::cout << "\nWarning:  cannot parse elevations in GPX file \"" << d.gpxs[ig]
						<< "\"." << std::endl;
			}

		}
		catch (const std::exception& e)
		{
			// Ignore empty GPX files with a warning.
			std::cout << "\nWarning:  cannot parse GPX file \"" << d.gpxs[ig]
					<< "\"." << std::endl;
			std::cout << e.what() << std::endl;
			//return ERR_XML_PARSE;
		}

	}  // GPX file loop

	//std::cout << "d.lats = " << d.lats << std::endl;
	//std::cout << "d.lons = " << d.lons << std::endl;
	//std::cout << "d.eles = " << d.eles << std::endl;
	//std::cout << "d.iEndSeg = " << d.iEndSeg << std::endl;

	std::cout << "\nTotal number of track points = " << d.ntrkptsum << std::endl;

	return 0;
}

std::vector<unsigned int> getKernel(Settings& s)
{
	//std::cout << "getting kernel..." << std::endl;
	s.tr1 = 2 * s.r + 1;
	std::vector<unsigned int> kernel(s.tr1 * s.tr1);

	s.gkb = s.r / sqrt(-log(1.0 / s.gka));
	double gkb2 = s.gkb * s.gkb;
	//std::cout << "gkb, gkb2 = " << s.gkb << " " << gkb2 << "\n";

	int ix, iy, ik;
	unsigned int inc, kmax;

	kmax = 0;
	for (int dy = -s.r; dy <= s.r; dy++)
	{
		iy = dy + s.r;
		int iyk = s.tr1 * (s.tr1 - iy - 1);
		for (int dx = -s.r; dx <= s.r; dx++)
		{
			ix = dx + s.r;

			if (s.kernel == pyramid)
				inc = std::max(0, s.r - abs(dx) - abs(dy));
			else if (s.kernel == gaussian)
				inc = s.gka * exp(-(dx*dx + dy*dy) / gkb2);
			else if (s.kernel == cylinder)
				inc = sqrt(dx*dx + dy*dy) <= s.r ? 1: 0;
			else // cone
				inc = std::max(0, s.r - (int) sqrt(dx*dx + dy*dy));

			ik = iyk + ix;
			kernel[ik] = inc;
			kmax = std::max(kmax, inc);
		}
	}

	unsigned int ndigits = ceil(log10(kmax + 1));
	if (s.verb > 0)
	{
		std::cout << "kernel =\n";// << kernel << std::endl;
		for (int i = 0; i < s.tr1; i++)
		{
			for (int j = i * s.tr1; j < (i+1) * s.tr1; j++)
				printf(("%" + std::to_string(ndigits+1) + "d").c_str(), kernel[j]);
			std::cout << "\n";
		}
		std::cout << std::endl;
	}

	return kernel;
}

double getAvgLength(Settings& s, Data& d)
{
	double lensum = 0;
	unsigned int iseg = 0;
	for (int i = 0; i < d.ntrkptsum - 1; i++)
	{
		if (i == d.iEndSeg[iseg])
		{
			iseg++;
		}
		else
		{
			double len = sqrt(pow(d.lats[i+1] - d.lats[i], 2)
			                + pow(d.lons[i+1] - d.lons[i], 2));
			lensum += len;
			//std::cout << "len = " << len << "\n";
		}
	}

	double lenavg = lensum / d.ntrkptsum;
	double lenavgpix = lenavg / s.degPerPix;

	//std::cout << "lensum    = " << lensum    << std::endl;
	//std::cout << "lenavg    = " << lenavg    << std::endl;
	//std::cout << "lenavgpix / s.r = " << lenavgpix / s.r << std::endl;

	std::cout << "Average track length = " << lenavgpix
			<< " pixels" << std::endl;

	return lenavgpix;
}

void printBar(int np, const std::string pstr)
{
	std::cout << "|";
	for (int i = 0; i < np; i++)
		std::cout << pstr;
	std::cout << "|" << std::endl;
	std::cout << "|";
}

void printProgress(int ip, const std::string pstr)
{
	static int ip0 = -1; // this won't work with multiple progress bars
	if (ip > ip0)
	{
		std::cout << pstr << std::flush;
		ip0 = ip;
	}
}

std::vector<unsigned int> convolute(Settings& s, Data& d,
		std::vector<unsigned int>& kernel)
{
	s.nxy = s.nx * s.ny;
	std::vector<unsigned int> img(s.nxy);

	// Initialize
	std::fill(img.begin(), img.end(), 0);

	// Number of characters in progress bar
	int np = 70;
	const std::string pstr = "=";

	// Kernel scalar for map types other than heat
	double scalar = 1.0;
	double scaleMax = 1.e9;

	if (s.sampling == pointwise)
	{
		std::cout << "Convoluting pointwise..." << std::endl;
		printBar(np, pstr);
		for (int i = 0; i < d.ntrkptsum; i++)
		{
			printProgress(np * i / d.ntrkptsum, pstr);
			addKernel(s, d.lats[i], d.lons[i], img, kernel, scalar);
		}
	}
	else
	{
		std::cout << "Convoluting linearly..." << std::endl;
		printBar(np, pstr);
		unsigned int iseg = 0;
		double dd = 0.0;
		double dd0 = 0.0;

		int i = 0;
		while (i < d.ntrkptsum)
		{
			double dpix;
			printProgress(np * i / d.ntrkptsum, pstr);
			if (i == d.iEndSeg[iseg])
			{
				// Final point
				iseg++;
				i++;
				dd = 0.0;
			}
			else
			{
				//std::cout << "sampling dd = " << dd << std::endl;

				// Linear subsampling

				double x0 = d.lons[i];
				double y0 = d.lats[i];
				double x1 = d.lons[i + 1];
				double y1 = d.lats[i + 1];

				double s0, s1;
				if (s.mapType != heat)
				{
					s0 = d.scas[i];
					s1 = d.scas[i + 1];
				}

				// Distance between trackpoint i and i+1 in units of
				// pixels.  Already calculated?
				dpix = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)) / s.degPerPix;

				while (dd < 1.0)
				{
					//std::cout << "dd = " << dd << std::endl;

					double lat = y0 + dd * (y1 - y0);
					double lon = x0 + dd * (x1 - x0);

					if (s.mapType != heat)
					{
						double sca = s0 + dd * (s1 - s0);
						scalar = scaleMax * (sca - s.mins) / (s.maxs - s.mins);
						//std::cout << "scalar = " << scalar << "\n";
					}

					addKernel(s, lat, lon, img, kernel, scalar);

					dd0 = dd;
					if (dpix == 0.0)
					{
						i++;
						break;
					}
					else
					{
						// Move a step of a fraction of the kernel radius
						dd += s.step * s.r / dpix;
					}

				}
			}
			//std::cout << "d.iEndSeg[iseg] = " << d.iEndSeg[iseg] << "\n";

			// Offset next kernel one step size away from last one
			double dpix0 = dpix;
			if (dpix0 != 0.0)
			{
				while (dd >= 1.0 && i < d.iEndSeg[iseg] && i < d.ntrkptsum)
				{
					i++;
					//std::cout << "i = " << i << std::endl;

					double x0 = d.lons[i];
					double y0 = d.lats[i];
					double x1 = d.lons[i + 1];
					double y1 = d.lats[i + 1];

					if (x0 != x1 || y0 != y1)
					{
						dpix0 = dpix;
						dpix = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)) / s.degPerPix;

						// Residual distance leftover from the last
						// line, as a fraction of last line's length
						double ddres = 1.0 - dd0;

						// Residual in pixels
						double pixres = ddres * dpix0;

						// Offset from new starting point in pixels
						double pixstart = s.step * s.r - pixres;

						// Dimensionless offset (and hypothetical
						// previous point in case of overflow)
						dd = pixstart / dpix;
						dd0 = dd - s.step * s.r / dpix;

						//std::cout << "dd = " << dd << std::endl;
					}
				}
			}

			// Just in case
			dd = std::max(dd, 0.0);

		}
	}
	std::cout << "|" << std::endl;

	//std::cout << img << std::endl;
	return img;
}

void setTransformation(Transformation& t)
{
	if (t.enable)
	{
		t.mat.resize(4);

		// Get rotation matrix
		t.rot *= pi / 180.0;  // degrees to radians
		t.mat[0] =  std::cos(t.rot);
		t.mat[1] = -std::sin(t.rot);
		t.mat[2] =  std::sin(t.rot);
		t.mat[3] =  std::cos(t.rot);

		// Make output directory
		#if defined(_WIN32) || defined(_WIN64)
			_mkdir(t.outdir.c_str());
		#else
			mkdir(t.outdir.c_str(), 0777);
		#endif
	}
}

void setFitting(Settings& s, Data& d)
{
	// Get the image bounds depending on which input options were
	// specified, and from the zoom level, set the sampling
	// strategy.

	if (s.fit)
	{
		s.minx = *min_element(begin(d.lons), end(d.lons));
		s.maxx = *max_element(begin(d.lons), end(d.lons));
		s.miny = *min_element(begin(d.lats), end(d.lats));
		s.maxy = *max_element(begin(d.lats), end(d.lats));

		if (s.mapType != heat)
		{
			s.mins = *min_element(begin(d.scas), end(d.scas));
			s.maxs = *max_element(begin(d.scas), end(d.scas));
		}
	}

	if (s.fitnx)
	{
		double avgy = 0.5 * (s.miny + s.maxy);
		double difx = s.maxx - s.minx;
		double dify = s.maxy - s.miny;
		s.nx = s.ny * cos(avgy * pi / 180.0) * difx / dify;
	}

	if (s.fitny)
	{
		double avgy = 0.5 * (s.miny + s.maxy);
		double difx = s.maxx - s.minx;
		double dify = s.maxy - s.miny;
		s.ny = s.nx / cos(avgy * pi / 180.0) * dify / difx;
	}

	if (s.fitx)
	{
		double avgx = 0.5 * (s.minx + s.maxx);
		double avgy = 0.5 * (s.miny + s.maxy);
		double dify = s.maxy - s.miny;
		double difx = dify / cos(avgy * pi / 180.0) * s.nx / s.ny;
		s.minx = avgx - 0.5 * difx;
		s.maxx = avgx + 0.5 * difx;
	}

	if (s.fity)
	{
		double avgy = 0.5 * (s.miny + s.maxy);
		double difx = s.maxx - s.minx;
		double dify = cos(avgy * pi / 180.0) * s.ny / s.nx * difx;
		s.miny = avgy - 0.5 * dify;
		s.maxy = avgy + 0.5 * dify;
	}

	double difx = s.maxx - s.minx;
	double dify = s.maxy - s.miny;
	s.minx -= s.margin * difx;
	s.maxx += s.margin * difx;
	s.miny -= s.margin * dify;
	s.maxy += s.margin * dify;

	if (s.fitnx || s.fitny)
	{
		std::cout << "Resizing image to " << s.nx << ", " << s.ny << std::endl;
	}

	if (s.fit || s.fitx || s.fity || s.margin != 0.0)
		printBounds(s);

	if (s.mapType != heat)
		std::cout << "\t" << s.mins << " <= scalar <= " << s.maxs << std::endl;

	s.degPerPix
			= sqrt(pow(s.maxx - s.minx, 2) + pow(s.maxy - s.miny, 2))
			/ sqrt(pow(s.nx           , 2) + pow(s.ny           , 2));

	if (s.sampling == autosample)
	{
		double lenavgpix = getAvgLength(s, d);
		if (lenavgpix / s.r < 0.5)  // TODO:  this limit depends on kernel type
			s.sampling = pointwise;
		else
			s.sampling = linear;
	}
}

int maph(int argc, char* argv[])
{
	int io = 0;

	std::cout << std::setprecision(16);

	Settings s;
	Transformation t;
	Data d;
	std::string fjson;

	io = loadArgs(argc, argv, s, t, fjson);
	if (io != 0)
		return io;

	json inj;
	io = loadJson(fjson, inj);
	if (io != 0)
		return io;

	io = loadSettings(s, inj, fjson);
	if (io != 0)
		return io;

	getFiles(s.fgpx, d.gpxs);

	//std::cout << "GPX files = " << d.gpxs << std::endl;
	std::cout << "Number of GPX files = " << d.gpxs.size() << std::endl;

	setTransformation(t);

	io = loadGpxs(s, t, d);
	if (io != 0)
		return io;

	if (t.enable)
		return 0;

	setFitting(s, d);
	auto kernel = getKernel(s);
	auto img = convolute(s, d, kernel);
	//return 0;  // for benchmarking kernel convolution

	// Sort for histogram coloring.
	auto idx = irwincolor::sortidx(img);

	if (s.allCmaps)
	{
		for (int i = 0; i < s.mapNames.size(); i++)
		{
			if (io = s.c.load(s.cmapFile, s.mapNames[i]) != 0)
				return io;
			auto pixels = colorPixels(s, img, idx);
			if (io = irwincolor::savePng(pixels, s.nx, s.ny,
					s.fname + s.dlm + s.mapNames[i] + s.imgext)
					!= 0)
				return io;
		}
	}
	else
	{
		auto pixels = colorPixels(s, img, idx);
		io = irwincolor::savePng(pixels, s.nx, s.ny, s.fname + s.imgext);
	}

	return io;
}

