
// TODO:
//
//     - Refactor
//     - Gaussian kernel
//     - Optionally use all colormaps from a colormap file
//     - More JSON inputs
//         * 'Fit nx' and 'Fit ny' options to change the number of
//           pixels for a consistent aspect instead of changing
//           lats/lons
//         * margin option to expand bounds (regardless of fitting options)
//         * subsampling step size
//         * kernel type
//         * background color:  0, NaN, or other
//     - Error checking
//

//======================================================================

// Standard

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>

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

#include <lodepng.h>

//======================================================================

// Custom

#include <colormap.h>
#include <xml_helper.h>
#include <constants.h>
#include <sort.h>

//======================================================================

const std::string me = "maph";

class Settings
{
	public:

		int nx, ny, nxy;

		bool fit, fitx, fity;
		double minx, maxx, miny, maxy;

		// Kernel radius (pixels)
		int r;

		// Kernel width
		int tr1;

		// isample == 0:  no subsampling, pointwise only
		// isample == 1:  automatically choose sampling based on avg length
		// isample == 2:  linear subsampling within every segment
		int isample;

		std::string fname;
		std::string fgpx;
		int verb;

		ColorMap c;

		// Bytes per pixel:  RGBA
		const int bpp = 4;

		const std::string imgext = ".png";
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

int savePng(const std::vector<uint8_t>& b, int nx, int ny, std::string f)
{
	std::cout << "Writing file \"" << f << "\" ..." << std::endl;
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

void getFiles(std::string &pattern, std::vector<std::string> &fileList)
{
#if defined(_WIN32) || defined(_WIN64)

	std::replace(pattern.begin(), pattern.end(), '/', '\\');

	// Wildcard globbing only works on Windows after the last
	// directory separator and only with a single wildcard.
	std::string dir = pattern.substr(0, pattern.find_last_of("\\"));

	// There are some fancy Windows API methods for calling dir, but
	// they don't work with wildcards!
	std::string ftmp = std::tmpnam(nullptr);
	std::string cmd = (std::string) "dir /b " + pattern + " > " + ftmp;
	system(cmd.c_str());

	std::ifstream ifs(ftmp);
	std::string line;
	while (std::getline(ifs, line))
	{
		// Dir lists the filename but not its path, so we need to
		// add it back.
		fileList.push_back(dir + "\\" + line);
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
			//// Map 0 to NaN color.  Otherwise evenly distribute.
			//x = 2.0;

			// Leave 0 as is.
			x = 0.0;
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
		pix[ib + 3] = twofivefive;  // alpha
	}
	return pix;
}

void addKernel(const Settings& s, double lat, double lon,
		std::vector<unsigned int>& img, std::vector<unsigned int>& kernel)
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
					img[ip] += kernel[ik];
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

int loadArgs(int argc, char* argv[], /* Settings& s, */ Transformation& t, std::string& fjson)
{
	t.enable = false;

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
	std::ifstream ifs(fjson);

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

	return 0;
}

int loadSettings(Settings& s, json& inj, std::string& fjson)
{
	// Initial defaults
	s.nx = 1280;
	s.ny = 720;
	s.verb = 0;
	s.isample = 1;
	s.r = 10;
	s.minx = 0.0;
	s.miny = 0.0;
	s.maxx = 0.0;
	s.maxy = 0.0;
	s.fitx = false;
	s.fity = false;
	s.fname = fjson.substr(0, fjson.find_last_of("."));
	s.fgpx = "*.gpx";
	std::string cmapFile = "";
	std::string cmapName = "";
	s.c.imap = -1;
	s.c.inv = false;

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
	const std::string fitxId = "Fit x";
	const std::string fityId = "Fit y";
	const std::string fnameId = "File name prefix";
	const std::string fgpxId = "GPX files";
	const std::string cmapFileId = "Colormap file";
	const std::string cmapNameId = "Colormap name";
	const std::string imapId = "Color map";
	const std::string invertMapId = "Invert color map";

	for (json::iterator it = inj.begin(); it != inj.end(); it++)
	{
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
		else if (it.key() == sampleId && it.value().is_number())
		{
			s.isample = it.value();
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
		else if (it.key() == fitxId && it.value().is_boolean())
		{
			s.fitx = it.value();
		}
		else if (it.key() == fityId && it.value().is_boolean())
		{
			s.fity = it.value();
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
			cmapFile = it.value();
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
		else
		{
			std::cout << "\nWarning:  unknown JSON key \"" << it.key()
					<< "\" with value \"" << it.value() << "\"\n";
		}
	}

	s.fit = !(bminx && bminy && bmaxx && bmaxy);

	// Echo inputs
	std::cout << "\n";
	std::cout << "Image size" << " = " << s.nx << ", " << s.ny << "\n";
	std::cout << fnameId << " = \"" << s.fname << "\"\n";
	std::cout << fgpxId << " = \"" << s.fgpx << "\"\n";
	std::cout << cmapFileId << " = \"" << cmapFile << "\"\n";
	std::cout << cmapNameId << " = \"" << cmapName << "\"\n";
	std::cout << imapId << " = " << s.c.imap << "\n";
	std::cout << invertMapId << " = " << s.c.inv << "\n";
	std::cout << verbId << " = " << s.verb << "\n";
	std::cout << sampleId << " = " << s.isample << "\n";
	std::cout << radiusId << " = " << s.r << "\n";
	std::cout << fitxId << " = " << s.fitx << "\n";
	std::cout << fityId << " = " << s.fity << "\n";

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
	return 0;
}

int loadGpxs(Settings&s, Transformation& t, std::vector<std::string>& gpxs,
		std::vector<double>& lats, std::vector<double>& lons,
		std::vector<unsigned int>& iEndSeg, int& ntrkptsum)
// TODO:  encapsulate some of these args.  c.f. similar functions.
{
	ntrkptsum = 0;
	double lat, lon;
	for (int ig = 0; ig < gpxs.size(); ig++)
	{
		int ntrkpt = 0;

		if (s.verb > 0) std::cout << "GPX file = " << gpxs[ig] << std::endl;

		// TODO:  make these next few lines a function in xml_helper.
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(gpxs[ig].c_str());
		if (!result)
		{
			std::cout << "\nError:  cannot open or parse GPX file \""
					<< gpxs[ig] << "\"." << std::endl;
			return ERR_XML_OPEN;
		}

		std::string xquery;
		pugi::xpath_node trkpt0;
		try
		{
			xquery = "/gpx/trk/trkseg/trkpt";
			trkpt0 = doc.select_node(xquery.c_str());
			if (!trkpt0) throw xPathException;

			for (pugi::xml_node trkpt = trkpt0.node(); trkpt;
					trkpt = trkpt.next_sibling("trkpt"))
			{
				ntrkpt++;
				lat = std::stod(trkpt.attribute("lat").value());
				lon = std::stod(trkpt.attribute("lon").value());
				//std::cout << "lat, lon = " << lat << ", " << lon << "\n";

				lats.push_back(lat);
				lons.push_back(lon);
			}

			if (s.verb > 0) std::cout << "\tNumber of track points = " << ntrkpt << std::endl;

			if (t.enable)
			{
				std::string fogpx = t.outdir + slash
						+ gpxs[ig].substr(gpxs[ig].find_last_of(slash) + 1);

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
					lons[j] += t.x;
					lats[j] += t.y;

					// Rotate
					lon = t.mat[0] * lons[j] + t.mat[1] * lats[j];
					lat = t.mat[2] * lons[j] + t.mat[3] * lats[j];

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
			iEndSeg.push_back(ntrkptsum - 1);

		}
		catch (const std::exception& e)
		{
			// Ignore empty GPX files with a warning.
			std::cout << "\nWarning:  cannot parse GPX file \"" << gpxs[ig]
					<< "\"." << std::endl;
			std::cout << e.what() << std::endl;
			//return ERR_XML_PARSE;
		}
	}

	//std::cout << "lats = " << lats << std::endl;
	//std::cout << "lons = " << lons << std::endl;
	//std::cout << "iEndSeg = " << iEndSeg << std::endl;

	std::cout << "\nTotal number of track points = " << ntrkptsum << std::endl;

	return 0;
}

std::vector<unsigned int> getKernel(Settings& s)
{
	//std::cout << "setting kernel..." << std::endl;

	s.tr1 = 2 * s.r + 1;
	std::vector<unsigned int> kernel(s.tr1 * s.tr1);
	int ix, iy, ik;
	unsigned int inc;
	for (int dy = -s.r; dy <= s.r; dy++)
	{
		iy = dy + s.r;
		int iyk = s.tr1 * (s.tr1 - iy - 1);
		for (int dx = -s.r; dx <= s.r; dx++)
		{
			ix = dx + s.r;

			//// Pyramid kernel
			//inc = std::max(0, s.r - abs(dx) - abs(dy));

			// Cone kernel, better looking
			inc = std::max(0, s.r - (int) sqrt(dx*dx + dy*dy));

			ik = iyk + ix;
			kernel[ik] = inc;
		}
	}
	return kernel;
}

double getAvgLength(std::vector<double>& lats, std::vector<double>& lons,
		std::vector<unsigned int> iEndSeg, int& ntrkptsum, double& degPerPix)
{
	double lensum = 0;
	unsigned int iseg = 0;
	for (int i = 0; i < ntrkptsum - 1; i++)
	{
		if (i == iEndSeg[iseg])
		{
			iseg++;
		}
		else
		{
			double len = sqrt(pow(lats[i+1] - lats[i], 2)
			                + pow(lons[i+1] - lons[i], 2));
			lensum += len;
			//std::cout << "len = " << len << "\n";
		}
	}

	double lenavg = lensum / ntrkptsum;
	double lenavgpix = lenavg / degPerPix;

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

std::vector<unsigned int> convolute(Settings& s, std::vector<double>& lats,
		std::vector<double>& lons, std::vector<unsigned int>& iEndSeg,
		int& ntrkptsum, double& degPerPix, std::vector<unsigned int>& kernel)
{
	s.nxy = s.nx * s.ny;
	std::vector<unsigned int> img(s.nxy);

	// Initialize
	std::fill(img.begin(), img.end(), 0);

	// Number of characters in progress bar
	int np = 70;
	const std::string pstr = "=";

	if (s.isample == 0)
	{
		std::cout << "Convoluting pointwise..." << std::endl;
		printBar(np, pstr);
		for (int i = 0; i < ntrkptsum; i++)
		{
			printProgress(np * i / ntrkptsum, pstr);
			addKernel(s, lats[i], lons[i], img, kernel);
		}
	}
	else // if (s.isample == 2)
	{
		std::cout << "Convoluting linearly..." << std::endl;
		printBar(np, pstr);
		unsigned int iseg = 0;
		for (int i = 0; i < ntrkptsum; i++)
		{
			printProgress(np * i / ntrkptsum, pstr);
			if (i == iEndSeg[iseg])
			{
				// Final point
				addKernel(s, lats[i], lons[i], img, kernel);
				iseg++;
			}
			else
			{
				// Linear subsampling

				double x0 = lons[i];
				double y0 = lats[i];
				double x1 = lons[i + 1];
				double y1 = lats[i + 1];

				double dpix = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)) / degPerPix;
				int n = std::max((int) (4 * (dpix / s.r)), 1);

				for (int j = 0; j < n; j++)
				{
					double lat = y0 + ((double) j / n) * (y1 - y0);
					double lon = x0 + ((double) j / n) * (x1 - x0);
					addKernel(s, lat, lon, img, kernel);
				}
			}
		}
	}
	std::cout << "|" << std::endl;

	//std::cout << img << std::endl;
	return img;
}

int maph(int argc, char* argv[])
{
	int io = 0;

	std::cout << std::setprecision(16);

	Settings s;
	Transformation t;
	std::string fjson;

	io = loadArgs(argc, argv, /* s, */ t, fjson);
	if (io != 0)
		return io;

	json inj;
	io = loadJson(fjson, inj);
	if (io != 0)
		return io;

	io = loadSettings(s, inj, fjson);
	if (io != 0)
		return io;

	std::vector<std::string> gpxs;
	getFiles(s.fgpx, gpxs);

	//std::cout << "GPX files = " << gpxs << std::endl;
	std::cout << "Number of GPX files = " << gpxs.size() << std::endl;

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

	std::vector<double> lats, lons;
	std::vector<unsigned int> iEndSeg;
	int ntrkptsum;
	io = loadGpxs(s, t, gpxs, lats, lons, iEndSeg, ntrkptsum);
	if (io != 0)
		return io;

	if (t.enable)
		return 0;

	if (s.fit)
	{
		s.minx = *min_element(begin(lons), end(lons));
		s.maxx = *max_element(begin(lons), end(lons));
		s.miny = *min_element(begin(lats), end(lats));
		s.maxy = *max_element(begin(lats), end(lats));
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

	if (s.fit || s.fitx || s.fity)
		printBounds(s);

	double degPerPix
			= sqrt(pow(s.maxx - s.minx, 2) + pow(s.maxy - s.miny, 2))
			/ sqrt(pow(s.nx           , 2) + pow(s.ny           , 2));

	if (s.isample == 1)
	{
		double lenavgpix = getAvgLength(lats, lons, iEndSeg, ntrkptsum, degPerPix);
		if (lenavgpix / s.r < 0.5)
			s.isample = 0;
		else
			s.isample = 2;
	}

	auto kernel = getKernel(s);
	auto img = convolute(s, lats, lons, iEndSeg, ntrkptsum, degPerPix, kernel);
	//return 0;  // for benchmarking kernel convolution

	// Sort for histogram coloring.
	auto idx = sortidx(img);

	auto pixels = colorPixels(s, img, idx);
	io = savePng(pixels, s.nx, s.ny, s.fname + s.imgext);
	return io;
}

