
// TODO:
//
//     - Refactor
//     - Gaussian kernel?
//     - More JSON inputs
//         * subsampling step size
//         * kernel type
//         * background color:  0, NaN, or other
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

class Settings
{
	public:

		int nx, ny, nxy;

		bool fit, fity;
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
	int iy0 = floor(s.ny * (lat - s.miny) / (s.maxy - s.miny));

	for (int dx = -s.r; dx <= s.r; dx++)
	{
		int ix = ix0 + dx;
		if (0 <= ix && ix < s.nx)
		{
			int dxr = dx + s.r;
			for (int dy = -s.r; dy <= s.r; dy++)
			{
				int iy = iy0 + dy;
				if (0 <= iy && iy < s.ny)
				{
					int dyr = dy + s.r;
					int ik = s.tr1 * (s.tr1 - dyr - 1) + dxr;
					int ip = s.nx  * (s.ny  - iy  - 1) + ix ;
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
	std::cout << "\nLoading JSON input file \"" << fjson << "\" ..." << std::endl;
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
	const std::string sizexId = "Image size x";
	s.nx = loadJsonOrDefault(sizexId, 1280, inj);

	const std::string sizeyId = "Image size y";
	s.ny = loadJsonOrDefault(sizeyId, 720, inj);

	const std::string verbId = "Verbosity";
	s.verb = loadJsonOrDefault(verbId, 0, inj);

	const std::string sampleId = "Sampling";
	s.isample = loadJsonOrDefault(sampleId, 1, inj);

	const std::string radiusId = "Kernel radius";
	s.r = loadJsonOrDefault(radiusId, 10, inj);

	const std::string minxId = "Min x";
	s.minx = loadJsonOrDefault(minxId, 0.0, inj);

	const std::string minyId = "Min y";
	s.miny = loadJsonOrDefault(minyId, 0.0, inj);

	const std::string maxxId = "Max x";
	s.maxx = loadJsonOrDefault(maxxId, 0.0, inj);

	const std::string maxyId = "Max y";
	s.maxy = loadJsonOrDefault(maxyId, 0.0, inj);

	const std::string fityId = "Fit y";
	s.fity = loadJsonOrDefault(fityId, false, inj);

	s.fit = s.minx == s.maxx && s.miny == s.maxy;

	const std::string fnameId = "File name prefix";
	const std::string fnameDflt = fjson.substr(0, fjson.find_last_of("."));
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
	std::cout << sampleId << " = " << s.isample << "\n";
	std::cout << radiusId << " = " << s.r << "\n";
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
	for (int dx = -s.r; dx <= s.r; dx++)
	{
		ix = dx + s.r;
		for (int dy = -s.r; dy <= s.r; dy++)
		{
			iy = dy + s.r;

			//// Pyramid kernel
			//inc = std::max(0, s.r - abs(dx) - abs(dy));

			// Cone kernel, better looking
			inc = std::max(0, s.r - (int) sqrt(dx*dx + dy*dy));

			ik = s.tr1 * (s.tr1 - iy - 1) + ix;
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

std::vector<unsigned int> convolute(Settings& s, std::vector<double>& lats,
		std::vector<double>& lons, std::vector<unsigned int>& iEndSeg,
		int& ntrkptsum, double& degPerPix, std::vector<unsigned int>& kernel)
{
	s.nxy = s.nx * s.ny;
	std::vector<unsigned int> img(s.nxy);

	// Initialize
	std::fill(img.begin(), img.end(), 0);

	if (s.isample == 0)
	{
		std::cout << "Convoluting pointwise..." << std::endl;
		for (int i = 0; i < ntrkptsum; i++)
			addKernel(s, lats[i], lons[i], img, kernel);
	}
	else // if (s.isample == 2)
	{
		std::cout << "Convoluting linearly..." << std::endl;
		unsigned int iseg = 0;
		for (int i = 0; i < ntrkptsum; i++)
		{
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

	if (s.fity)
	{
		double avgy = 0.5 * (s.miny + s.maxy);
		double difx = s.maxx - s.minx;
		double dify = cos(avgy * pi / 180.0) * s.ny / s.nx * difx;
		s.miny = avgy - 0.5 * dify;
		s.maxy = avgy + 0.5 * dify;
	}

	if (s.fit || s.fity)
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
	io = save_png(pixels, s.nx, s.ny, s.fname + s.imgext);
	return io;
}
