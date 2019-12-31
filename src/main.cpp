
// TODO:
//
//     - Precompute kernel for optimization
//     - Gaussian kernel?
//     - Move kernel along line between trackpoints?
//     - Shuffling
//     - Add test suite comparing expected PNG hash
//     - More JSON inputs
//     - Refactor
//

//======================================================================

// Standard

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
//#include <regex>
//#include <experimental/filesystem>
#include <glob.h>

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

#include <lodepng.h>

//======================================================================

// Custom

#include <colormap.h>
#include <xml_helper.h>
#include <constants.h>
#include <sort.h>

//======================================================================

const std::string me = "maph";

template<class T> std::ostream& operator<<(std::ostream& stream, const std::vector<T>& values)
{
	// cout (or other stream) << a vector
	stream << "{";
	copy(begin(values), end(values), std::ostream_iterator<T>(stream, ", "));
	stream << '}';
	return stream;
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

class Settings
{
	public:

	int nx, ny, nxy;

	bool fit;
	double minx, maxx, miny, maxy;

	std::string fname;
	std::string fgpx;
	int verb;

	ColorMap c;

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
		std::vector<uint8_t> rgb = s.c.map(x);

		//std::cout << "rgb = " << rgb[0] << " " << rgb[1] << " " << rgb[2] << "\n";

		// Save RGB triple to pixel.
		int ib = idx[i] * s.bpp;
		b[ib + 0] = rgb[0];
		b[ib + 1] = rgb[1];
		b[ib + 2] = rgb[2];
		b[ib + 3] = twofivefive;     // alpha
	}
}

int usage()
{
	std::cout << "Usage:" << std::endl;
	std::cout << "\t" << me << " input.json [-t dir xshift yshift rot]" << std::endl;
	return ERR_CMD_ARG;
}

int maph(int argc, char* argv[])
{
	std::cout << std::setprecision(16);
	std::string fjson;

	std::string transdir;
	bool trans = false;
	double transx, transy, transd;

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
				trans = true;
				transdir = argv[++i];
				std::cout << "Transformation directory = \"" << transdir << "\".\n";
				transx = std::stod(argv[++i]);
				transy = std::stod(argv[++i]);
				transd = std::stod(argv[++i]);
				std::cout << "X translation = " << transx << "\n";
				std::cout << "Y translation = " << transy << "\n";
				std::cout << "Rotation      = " << transd << "\n";
				std::cout << std::endl;
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

	const std::string cmapFileId = "Colormap file";
	const std::string cmapFileDflt = "";
	std::string cmapFile = loadJsonOrDefault(cmapFileId, cmapFileDflt, inj);

	const std::string cmapNameId = "Colormap name";
	const std::string cmapNameDflt = "";
	std::string cmapName = loadJsonOrDefault(cmapNameId, cmapNameDflt, inj);

	const std::string imapId = "Color map";
	s.c.imap = loadJsonOrDefault(imapId, -1, inj);

	const std::string invertMapId = "Invert color map";
	s.c.inv = loadJsonOrDefault(invertMapId, false, inj);

	// Echo inputs
	std::cout << "Image size" << " = " << s.nx << ", " << s.ny << "\n";
	std::cout << fnameId << " = \"" << s.fname << "\"\n";
	std::cout << fgpxId << " = \"" << s.fgpx << "\"\n";
	std::cout << cmapFileId << " = \"" << cmapFile << "\"\n";
	std::cout << cmapNameId << " = \"" << cmapName << "\"\n";
	std::cout << imapId << " = " << s.c.imap << "\n";
	std::cout << invertMapId << " = " << s.c.inv << "\n";
	std::cout << verbId << " = " << s.verb << "\n";

	//std::cout << "s.fit = " << s.fit << "\n";
	if (!s.fit)
		printBounds(s);

	std::cout << std::endl;

	s.c.paraView = cmapFile != "" && cmapName != "";
	if (s.c.paraView)
	{
		int io = s.c.load(cmapFile, cmapName);
		if (io != 0) return io;
	}

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

	std::vector<double> rmat(4);
	if (trans)
	{
		// Get rotation matrix
		transd *= pi / 180.0;  // degrees to radians
		rmat[0] =  std::cos(transd);
		rmat[1] = -std::sin(transd);
		rmat[2] =  std::sin(transd);
		rmat[3] =  std::cos(transd);

		// Make output directory
		#if defined(_WIN32) || defined(_WIN64)
			_mkdir(transdir.c_str());
		#else
			mkdir(transdir.c_str(), 0777);
		#endif

	}

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

			if (trans)
			{
				// TODO:  Windows
				std::string fogpx = transdir + "/" + gpxs[ig].substr(gpxs[ig].find_last_of("/") + 1);

				if (s.verb > 0) std::cout << "fogpx = " << fogpx << std::endl;

				std::ofstream ofgpx(fogpx);
				ofgpx << std::setprecision(16);
				ofgpx << "<gpx>\n";
				ofgpx << "\t<trk>\n";
				ofgpx << "\t\t<trkseg>\n";

				// Shuffle
				std::vector<unsigned int> id(ntrkpt);
				for (unsigned int i = 0; i < ntrkpt; i++)
					id[i] = i;
				std::random_shuffle(id.begin(), id.end());
				//std::cout << "id = " << id << std::endl;

				for (int it = 0; it < ntrkpt; it++)
				{
					unsigned int j = id[it];

					// Translate
					lons[j] += transx;
					lats[j] += transy;

					// Rotate
					lon = rmat[0] * lons[j] + rmat[1] * lats[j];
					lat = rmat[2] * lons[j] + rmat[3] * lats[j];

					ofgpx << "\t\t\t<trkpt lat=\"" << lat
					                << "\" lon=\"" << lon << "\"/>\n";
				}

				lats.clear();
				lats.shrink_to_fit();
				lons.clear();
				lons.shrink_to_fit();

				ofgpx << "\t\t</trkseg>\n";
				ofgpx << "\t</trk>\n";
				ofgpx << "</gpx>\n";

				ofgpx.close();
			}

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

	if (trans)
		return 0;

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
						//// Pyramid kernel
						//inc = std::max(0, r - abs(dx) - abs(dy));

						// Cone kernel.  A bit slower but better looking.
						inc = std::max(0, r - (int) sqrt(dx*dx + dy*dy));

						ip = s.nx * (s.ny - iy - 1) + ix;
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
	std::cout << std::endl;
	std::cout << "Starting " << me << std::endl;
	std::cout << std::endl;

	int io = maph(argc, argv);

	std::cout << "Done " << me << std::endl;
	std::cout << std::endl;
	return io;
}

