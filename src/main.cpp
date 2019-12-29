#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <regex>
#include <experimental/filesystem>
#include <glob.h>

#include <json.hpp>
using json = nlohmann::json;

#include <pugixml.hpp>

#include <lodepng.h>

const double pi = 4.0 * atan(1.0);

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

class Settings
{
	public:

		int nx, ny, nxy;

		bool fit;
		double minx, maxx, miny, maxy;

		std::string fname;
		std::string fgpx;

		// Bytes per pixel:  RGBA
		const int bpp = 4;

		const std::string dlm = "_", imgext = ".png";
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
	const int sizexDflt = 1280;
	s.nx = loadJsonOrDefault(sizexId, sizexDflt, inj);

	const std::string sizeyId = "Image size y";
	const int sizeyDflt = 720;
	s.ny = loadJsonOrDefault(sizeyId, sizeyDflt, inj);

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
	std::cout << "Image size" << " = " << s.nx << " " << s.ny << "\n";
	std::cout << fnameId << " = \"" << s.fname << "\"\n";
	std::cout << fgpxId << " = \"" << s.fgpx << "\"\n";

	//std::cout << "s.fit = " << s.fit << "\n";
	if (!s.fit)
		printBounds(s);

	std::cout << std::endl;

	std::vector<std::string> gpxs;
	getFiles(s.fgpx, gpxs);

	//std::cout << "GPX files = " << gpxs << std::endl;
	std::cout << "Number of GPX files = " << gpxs.size() << std::endl;

	std::vector<unsigned char> img;
	s.nxy = s.nx * s.ny;
	img.resize(s.nxy * s.bpp);

	// Initialize
	int ip = 0;
	int ix, iy, ib;
	for (ix = 0; ix < s.nx; ix++)
		for (iy = 0; iy < s.ny; iy++)
		{
			ip++;
			ib = ip * s.bpp;
			img[ib + 0] = 0;    // R
			img[ib + 1] = 0;    // G
			img[ib + 2] = 0;    // B
			img[ib + 3] = 255;  // alpha
		}

	std::vector<double> lats, lons;
	int ntrkptsum = 0;
	double lat, lon;
	for (int ig = 0; ig < gpxs.size(); ig++)
	{
		int ntrkpt = 0;

		// TODO:  add verbosity option
		std::cout << "GPX file = " << gpxs[ig] << std::endl;

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

			std::cout << "\tNumber of track points = " << ntrkpt << std::endl;
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

	// Line drawing
	for (i = 0; i < ntrkptsum; i++)
	{
		lat = lats[i];
		lon = lons[i];
		ix = floor(s.nx * (lon - s.minx) / (s.maxx - s.minx));
		iy = floor(s.ny * (lat - s.miny) / (s.maxy - s.miny));

		//std::cout << "ix, iy = " << ix << ", " << iy << std::endl;

		if (0 <= ix && ix < s.nx && 0 <= iy && iy < s.ny)
		{
			ip = s.nx * (s.ny - iy - 1) + ix;
			ib = ip * s.bpp;

			img[ib + 0] = 255;    // R
			img[ib + 1] = 255;    // G
			img[ib + 2] = 255;    // B
		}
	}

	int io = save_png(img, s.nx, s.ny, "test.png");

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
		std::cout << "Usage:" << std::endl;
		std::cout << "\t" << me << " input.json" << std::endl;
		return ERR_CMD_ARG;
	}

	int io = maph(argc, argv);

	std::cout << "Done " << me << std::endl;
	std::cout << std::endl;
	return io;
}

