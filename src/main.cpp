
// TODO:
//
//     - Cone kernel
//     - Gaussian kernel?
//     - More colormaps
//     - Delete lawnmowing activity
//     - Transformation
//     - Shuffling
//     - More JSON inputs
//     - Refactor
//

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <regex>
//#include <experimental/filesystem>
#include <glob.h>

#include <json.hpp>
using json = nlohmann::json;

#include <pugixml.hpp>

#include <lodepng.h>

const int twofivefive = 255;

// Error return codes
const int ERR_CMD_ARG   = -1;
const int ERR_JSON      = -2;
const int ERR_XML_OPEN  = -3;
const int ERR_XML_PARSE = -4;
const int ERR_PNG       = -5;

class XPathException: public std::exception
{
  virtual const char* what() const throw()
  {
    return "No match for XPath expression.";
  }
} xPathException;

template<class T> std::ostream& operator<<(std::ostream& stream, const std::vector<T>& values)
{
	// cout (or other stream) << a vector
	stream << "{";
	copy(begin(values), end(values), std::ostream_iterator<T>(stream, ", "));
	stream << '}';
	return stream;
}

template <typename T> std::vector<unsigned int> sortidx(const std::vector<T>& v)
{
	// Instead of sorting a vector v, sort an index array that gives
	// the ascending sorted values v[idx[0: ]].

	// initialize original index locations
	std::vector<unsigned int> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indices based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](unsigned int i1, unsigned int i2) {return v[i1] < v[i2];});

	return idx;
}

class Settings
{
	public:

		int nx, ny, nxy;

		bool fit;
		double minx, maxx, miny, maxy;

		std::string fname;
		std::string fgpx;
		int verb;

		// Bytes per pixel:  RGBA
		const int bpp = 4;

		const std::string imgext = ".png";
};

void printBounds(const Settings& s)
{
	std::cout << "Bounds:\n\t"
			<< s.minx << " <= x <= " << s.maxx
			<< "\n\t"
			<< s.miny << " <= y <= " << s.maxy
			<< std::endl;
}

template <typename T> T loadJsonOrDefault(std::string id, T dflt, json j)
{
	// Return a value from a JSON object j by its key id, or return
	// a default dflt if it cannot be parsed.
	T val;
	try
	{
		// TODO:  see https://stackoverflow.com/questions/16337610/how-to-know-if-a-type-is-a-specialization-of-stdvector

		val = j[id];
		//val = j[id].get<std::vector>;
	}
	catch (const std::exception& e)
	{
		std::cout << "Input did not contain properly formatted \"" << id << "\".  Using default." << std::endl;
		val = dflt;
	}
	return val;
}

int save_png(const std::vector<uint8_t>& b, int nx, int ny, std::string f)
{
	std::cout << "Writing file \"" << f << "\"..." << std::endl;
	std::vector<uint8_t> imageBuffer;
	unsigned error = lodepng::encode(imageBuffer, b, nx, ny);
	if (error)
	{
		std::cout << "\nError: PNG encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
		return ERR_PNG;
	}
	lodepng::save_file(imageBuffer, f);
	return 0;
}

// Linux only?
void getFiles(const std::string &pattern, std::vector<std::string> &fileList)
{
	//Declare glob_t for storing the results of globbing
	glob_t globbuf;

	// Glob. GLOB_TILDE tells the globber to expand "~" in the pattern to the home directory
	glob(pattern.c_str(), GLOB_TILDE, NULL, &globbuf);

	for (int i = 0; i < globbuf.gl_pathc; ++i)
		fileList.push_back(globbuf.gl_pathv[i]);

	// Free the globbuf structure
	if (globbuf.gl_pathc > 0)
		globfree(&globbuf);
}

std::vector<uint8_t> map(double y)
{
	// TODO:  inputs
	int imap = 6;
	bool inv = false;
	double nanr = 0.0, nang = 0.0, nanb = 0.0;

	// Almost, we need to floor result
	const float twofivesix = 255.999;

	// Map y in range [0, 1] (clipped if outside) to an RGB triple
	// in range [0, 255].

	std::vector<uint8_t> rgb(3);

	// Map out-of-range values to NaN color.
	if (y > 1.0)
	{
		rgb[0] = twofivesix * nanr;
		rgb[1] = twofivesix * nang;
		rgb[2] = twofivesix * nanb;
		return rgb;
	}

	y = std::min(std::max(y, 0.0), 1.0);

	if (inv)
		y = 1 - y;

	//// Note, using a function like sqrt here can skew the map
	//// towards the 255 end of the spectrum, while squaring can shift
	//// towards 0.
	//
	// TODO:  make an option for the exponent here.  It's useful for
	// things like black-to-red where the bright reds are
	// imperceptibly different from one another and a larger
	// exponent (e.g. 3 or 4) is helpful, while for HSV 1 is OK.
	//
	//y = sqrt(y);
	//y = y * y * y * y;
	if (0 <= imap && imap < 3)
	{
		// Black-to-(red|green|blue) map
		//            0    1     2
		rgb[imap] = y * twofivesix;
	}
	else if (3 <= imap && imap < 6)
	{
		// (cyan|magenta|yellow)-to-white map
		//    3     4       5
		rgb = {twofivefive, twofivefive, twofivefive};
		rgb[imap - 3] = y * twofivesix;
	}
	else
	{
		// HSV
		uint8_t c = twofivefive;
		double hp = 6.0 * y;
		uint8_t xu = twofivesix * (1.0 - abs(fmod(hp, 2.0) - 1.0));
		if (hp < 1.0)
			rgb = {c, xu, 0};
		else if (hp < 2.0)
			rgb = {xu, c, 0};
		else if (hp < 3.0)
			rgb = {0, c, xu};
		else if (hp < 4.0)
			rgb = {0, xu, c};
		else if (hp < 5.0)
			rgb = {xu, 0, c};
		else
			rgb = {c, 0, xu};
	}

	//std::cout << "y = " << y << "\n";
	//std::cout << "rgb = " << rgb << "\n";

	return rgb;
}

void colorPixels(const Settings& s, const std::vector<unsigned int>& img, const std::vector<unsigned int>& idx, std::vector<uint8_t>& b)
{
	double x, x0;
	x0 = 0.0;

	// Get min, max, number of non-zeros
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
			// Map 0 to NaN color.  Otherwise evenly distribute.
			// TODO:  add option to use 0 color without affecting histogram.
			x = 2.0;
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
		std::vector<uint8_t> rgb = map(x);

		//std::cout << "rgb = " << rgb[0] << " " << rgb[1] << " " << rgb[2] << "\n";

		// Save RGB triple to pixel.
		int ib = idx[i] * s.bpp;
		b[ib + 0] = rgb[0];
		b[ib + 1] = rgb[1];
		b[ib + 2] = rgb[2];
		b[ib + 3] = twofivefive;     // alpha
	}
}

int maph(int argc, char* argv[])
{

	std::cout << std::setprecision(16);
	std::string fjson;

	int i = 0;
	while (i < argc - 1)
	{
		i++;
		std::string str = argv[i];
		if (str == "-t")
		{
			// TODO:  transformation args
		}
		else
			fjson = str;
	}

	std::ifstream ifs(fjson);
	json inj;

	try
	{
		ifs >> inj;
		ifs.close();  // close to release lock on file
	}
	catch (const std::exception& e)
	{
		std::cout << "\nError:  cannot load JSON input file \"" << fjson << "\"." << std::endl;
		std::cout << e.what() << std::endl;
		return ERR_JSON;
	}

	Settings s;

	const std::string sizexId = "Image size x";
	s.nx = loadJsonOrDefault(sizexId, 1280, inj);

	const std::string sizeyId = "Image size y";
	s.ny = loadJsonOrDefault(sizeyId, 720, inj);

	const std::string verbId = "Verbosity";
	s.verb = loadJsonOrDefault(verbId, 0, inj);

	const std::string minxId = "Min x";
	s.minx = loadJsonOrDefault(minxId, 0.0, inj);

	const std::string minyId = "Min y";
	s.miny = loadJsonOrDefault(minyId, 0.0, inj);

	const std::string maxxId = "Max x";
	s.maxx = loadJsonOrDefault(maxxId, 0.0, inj);

	const std::string maxyId = "Max y";
	s.maxy = loadJsonOrDefault(maxyId, 0.0, inj);

	s.fit = s.minx == s.maxx && s.miny == s.maxy;

	const std::string fnameId = "File name prefix";
	const std::string fnameDflt = fjson;
	s.fname = loadJsonOrDefault(fnameId, fnameDflt, inj);

	const std::string fgpxId = "GPX files";
	const std::string fgpxDflt = "*.gpx";
	s.fgpx = loadJsonOrDefault(fgpxId, fgpxDflt, inj);

	// Echo inputs
	std::cout << "Image size" << " = " << s.nx << ", " << s.ny << "\n";
	std::cout << fnameId << " = \"" << s.fname << "\"\n";
	std::cout << fgpxId << " = \"" << s.fgpx << "\"\n";
	std::cout << verbId << " = " << s.verb << "\n";

	//std::cout << "s.fit = " << s.fit << "\n";
	if (!s.fit)
		printBounds(s);

	std::cout << std::endl;

	std::vector<std::string> gpxs;
	getFiles(s.fgpx, gpxs);

	//std::cout << "GPX files = " << gpxs << std::endl;
	std::cout << "Number of GPX files = " << gpxs.size() << std::endl;

	s.nxy = s.nx * s.ny;
	std::vector<unsigned int> img(s.nxy);
	std::vector<uint8_t> b(s.bpp * s.nxy);

	// Initialize
	int ip;
	for (ip = 0; ip < s.nxy; ip++)
		img[ip] = 0;

	std::vector<double> lats, lons;
	int ntrkptsum = 0;
	double lat, lon;
	for (int ig = 0; ig < gpxs.size(); ig++)
	{
		int ntrkpt = 0;

		if (s.verb > 0) std::cout << "GPX file = " << gpxs[ig] << std::endl;

		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(gpxs[ig].c_str());
		if (!result)
		{
			std::cout << "\nError:  cannot open or parse GPX file \"" << gpxs[ig] << "\"." << std::endl;
			return ERR_XML_OPEN;
		}

		std::string xquery;
		pugi::xpath_node trkpt0;
		try
		{
			xquery = "/gpx/trk/trkseg/trkpt";
			trkpt0 = doc.select_node(xquery.c_str());
			if (!trkpt0) throw xPathException;

			for (pugi::xml_node trkpt = trkpt0.node(); trkpt; trkpt = trkpt.next_sibling("trkpt"))
			{
				ntrkpt++;
				lat = std::stod(trkpt.attribute("lat").value());
				lon = std::stod(trkpt.attribute("lon").value());
				//std::cout << "lat, lon = " << lat << ", " << lon << "\n";

				lats.push_back(lat);
				lons.push_back(lon);
			}

			if (s.verb > 0) std::cout << "\tNumber of track points = " << ntrkpt << std::endl;

			ntrkptsum += ntrkpt;

		}
		catch (const std::exception& e)
		{
			// Ignore empty GPX files with a warning.
			std::cout << "\nWarning:  cannot parse GPX file \"" << gpxs[ig] << "\"." << std::endl;
			std::cout << e.what() << std::endl;
			//return ERR_XML_PARSE;
		}
	}

	//std::cout << "lats = " << lats << std::endl;
	//std::cout << "lons = " << lons << std::endl;

	std::cout << "\nTotal number of track points = " << ntrkptsum << std::endl;

	if (s.fit)
	{
		s.minx = *min_element(begin(lons), end(lons));
		s.maxx = *max_element(begin(lons), end(lons));
		s.miny = *min_element(begin(lats), end(lats));
		s.maxy = *max_element(begin(lats), end(lats));

		printBounds(s);
	}

	// Kernel radius (pixels)
	int r = 10;

	int ix, iy, ix0, iy0;
	unsigned int inc;
	for (i = 0; i < ntrkptsum; i++)
	{
		lat = lats[i];
		lon = lons[i];
		ix0 = floor(s.nx * (lon - s.minx) / (s.maxx - s.minx));
		iy0 = floor(s.ny * (lat - s.miny) / (s.maxy - s.miny));

		// Pyramid kernel
		for (int dx = -r; dx <= r; dx++)
		{
			ix = ix0 + dx;
			if (0 <= ix && ix < s.nx)
			{
				for (int dy = -r; dy <= r; dy++)
				{
					iy = iy0 + dy;
					if (0 <= iy && iy < s.ny)
					{
						ip = s.nx * (s.ny - iy - 1) + ix;
						inc = std::max(0, r - abs(dx) - abs(dy));
						img[ip] = img[ip] + inc;
					}
				}
			}
		}
	}
	//std::cout << img << std::endl;
	//return 0;  // for benchmarking kernel convolution

	// Sort for histogram coloring.
	auto idx = sortidx(img);

	colorPixels(s, img, idx, b);
	int io = save_png(b, s.nx, s.ny, s.fname + s.imgext);
	return io;
}

int main(int argc, char* argv[])
{
	std::string me = "maph";
	std::cout << std::endl;
	std::cout << "Starting " << me << std::endl;
	std::cout << std::endl;

	if (argc < 2)
	{
		// TODO:  move to maph()
		std::cout << "Usage:" << std::endl;
		std::cout << "\t" << me << " input.json" << std::endl;
		return ERR_CMD_ARG;
	}

	int io = maph(argc, argv);

	std::cout << "Done " << me << std::endl;
	std::cout << std::endl;
	return io;
}

