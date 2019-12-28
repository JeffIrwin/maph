#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h>

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
		std::string fname;

		// Bytes per pixel:  RGBA
		const int bpp = 4;

		const std::string dlm = "_", imgext = ".png";
};

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
		std::cout << "PNG encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
		return ERR_PNG;
	}
	lodepng::save_file(imageBuffer, f);
	return 0;
}

int maph(int argc, char* argv[])
{

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

	const std::string fnameId = "File name prefix";
	const std::string fnameDflt = fjson;
	s.fname = loadJsonOrDefault(fnameId, fnameDflt, inj);

	// Echo inputs
	std::cout << "Image size" << " = " << s.nx << " " << s.ny << "\n";
	std::cout << fnameId << " = \"" << s.fname << "\"\n";

	std::vector<unsigned char> img;
	s.nxy = s.nx * s.ny;
	img.resize(s.nxy * s.bpp);

	// Initialize
	int ip = 0;
	for (int ix = 0; ix < s.nx; ix++)
		for (int iy = 0; iy < s.ny; iy++)
		{
			ip++;
			int ib = ip * s.bpp;
			img[ib + 0] = 0;    // R
			img[ib + 1] = 0;    // G
			img[ib + 2] = 0;    // B
			img[ib + 3] = 255;  // alpha
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

