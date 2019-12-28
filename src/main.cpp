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

template<class T> std::ostream& operator<<(std::ostream& stream, const std::vector<T>& values)
{
	// cout (or other stream) << a vector
	stream << "{";
	copy(begin(values), end(values), std::ostream_iterator<T>(stream, ", "));
	stream << '}';
	return stream;
}

//int main(int argc, char* argv[])
int main()
{
	std::cout << "\nStarting maph" << std::endl;

	// Number of pixels

	//unsigned int nx = 420;
	//unsigned int ny = 300;

	//unsigned int nx = 1280;
	//unsigned int ny = 720;

	//unsigned int nx = 1920;
	//unsigned int ny = 1080;

	unsigned int nx = 3840;
	unsigned int ny = 2160;

	const unsigned int nbpp = 4;  // bytes per pixel

	std::vector<unsigned char> img;
	img.resize(nx * ny * nbpp);

	// Initialize
	int i = 0;
	for (int ix = 0; ix < nx; ix++)
		for (int iy = 0; iy < ny; iy++)
		{
			i++;
			int ib = i * nbpp;
			img[ib + 0] = 0;    // R
			img[ib + 1] = 0;    // G
			img[ib + 2] = 0;    // B
			img[ib + 3] = 255;  // alpha
		}

	std::string f = "test.png";
	std::cout << "Encoding and writing image \"" << f << "\"..." << std::endl;
	std::vector<unsigned char> buffer;
	unsigned error = lodepng::encode(buffer, &img[0], nx, ny);
	if (error)
	{
		std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
	}
	lodepng::save_file(buffer, f);

	std::cout << "\nDone!" << std::endl;
	return 0;
}
