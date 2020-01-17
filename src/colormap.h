
#ifndef INCLUDE_COLORMAP_H_
#define INCLUDE_COLORMAP_H_

//======================================================================

#include <vector>
#include <string>
#include <math.h>  // fmod

//======================================================================

#include <xml_helper.h>
#include <constants.h>
#include <sort.h>

//======================================================================

const unsigned char twofivefive = 255;

std::vector<float> rgb2xyz(const std::vector<float>& rgb)
{
	// https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
	return {
			0.4124f * rgb[0] + 0.3576f * rgb[1] + 0.1805f * rgb[2],
			0.2126f * rgb[0] + 0.7152f * rgb[1] + 0.0722f * rgb[2],
			0.0193f * rgb[0] + 0.1192f * rgb[1] + 0.9505f * rgb[2]
	};
}

std::vector<float> xyz2rgb(const std::vector<float>& xyz)
{
	std::vector<float> rgb(3);
	rgb[0] = 3.2406255f * xyz[0] - 1.5372080f * xyz[1] - 0.4986286f * xyz[2];
	rgb[1] =-0.9689307f * xyz[0] + 1.8757561f * xyz[1] + 0.0415175f * xyz[2];
	rgb[2] = 0.0557101f * xyz[0] - 0.2040211f * xyz[1] + 1.0569959f * xyz[2];
	rgb[0] = std::min(std::max(rgb[0], 0.0f), 1.0f);
	rgb[1] = std::min(std::max(rgb[1], 0.0f), 1.0f);
	rgb[2] = std::min(std::max(rgb[2], 0.0f), 1.0f);
	return rgb;
}

float labFunc(const float& x)
{
	// https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
	return x > 0.008856 ? pow(x, 1.0/3.0) : 7.787 * x + 16.0 / 116.0;
}

float labFuncInv(const float& f)
{
	return f > 0.206893 ? pow(f, 3) : (f - 16.0 / 116.0) / 7.787;
}

std::vector<float> xyz2lab(const std::vector<float>& xyz)
{
	// https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
	const float xn = 95.047, yn = 100.0, zn = 108.883;
	std::vector<float> lab(3);
	float fx = labFunc(xyz[0] / xn);
	float fy = labFunc(xyz[1] / yn);
	float fz = labFunc(xyz[2] / zn);
	lab[0] = 116.0 * (fy - 16.0 / 116.0);
	lab[1] = 500.0 * (fx - fy);
	lab[2] = 200.0 * (fy - fz);
	return lab;
}

std::vector<float> lab2xyz(const std::vector<float>& lab)
{
	std::vector<float> xyz(3);
	float fy = (lab[0] + 16.0) / 116.0;
	float fx = lab[1] / 500.0 + fy;
	float fz = fy - lab[2] / 200.0;
	const float xn = 95.047, yn = 100.0, zn = 108.883;
	xyz[0] = labFuncInv(fx) * xn;
	xyz[1] = labFuncInv(fy) * yn;
	xyz[2] = labFuncInv(fz) * zn;
	//std::cout << "lab2xyz: xyz = " << xyz << std::endl;
	return xyz;
}

std::vector<float> rgb2lab(const std::vector<float>& rgb)
{
	return xyz2lab(rgb2xyz(rgb));
}

std::vector<float> lab2msh(const std::vector<float>& lab)
{
	// https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf
	std::vector<float> msh(3);
	msh[0] = sqrt(lab[0] * lab[0] + lab[1] * lab[1] + lab[2] * lab[2]);
	msh[1] = acos(lab[0] / msh[0]);
	msh[2] = atan2(lab[2], lab[1]);
	return msh;
}

std::vector<float> rgb2msh(const std::vector<float>& rgb)
{
	return lab2msh(rgb2lab(rgb));
}

std::vector<float> lab2rgb(const std::vector<float>& lab)
{
	return xyz2rgb(lab2xyz(lab));
}

std::vector<float> msh2lab(const std::vector<float>& msh)
{
	return {
			msh[0] * cos(msh[1]),
			msh[0] * sin(msh[1]) * cos(msh[2]),
			msh[0] * sin(msh[1]) * sin(msh[2])
	};
}

std::vector<float> rgb2hsv(const std::vector<float>& rgb)
{
	std::vector<float> hsv(3);
	float m = std::min(std::min((double) rgb[0], (double) rgb[1]), (double) rgb[2]);
	float c, hp;
	if (rgb[0] > rgb[1] && rgb[0] > rgb[2])
	{
		c = rgb[0] - m;
		hsv[2] = rgb[0];
		hp = fmod((rgb[1] - rgb[2]) / c, 6.0);
	}
	else if (rgb[1] > rgb[2] && rgb[1] > rgb[0])
	{
		c = rgb[1] - m;
		hsv[2] = rgb[1];
		hp = (rgb[2] - rgb[0]) / c + 2.0;
	}
	else
	{
		c = rgb[2] - m;
		hsv[2] = rgb[2];
		hp = c <= 0.0 ? 0.0 : (rgb[0] - rgb[1]) / c + 4.0;
	}
	hsv[0] = hp * 60.0;
	hsv[1] = hsv[2] <= 0.0 ? 0.0 : c / hsv[2];
	return hsv;
}

std::vector<float> hsv2rgb(const std::vector<float>& hsv)
{
	std::vector<float> rgb(3);
	float c = std::min((double) hsv[1] * hsv[2], 1.0);
	float hp = hsv[0] / 60.0;
	float xu = std::min((double) c * (1.0 - abs(fmod(hp, 2.0) - 1.0)), 1.0);
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

	float m = hsv[2] - c;
	rgb[0] += m;
	rgb[1] += m;
	rgb[2] += m;

	rgb[0] = std::min(std::max((double) rgb[0], 0.0), 1.0);
	rgb[1] = std::min(std::max((double) rgb[1], 0.0), 1.0);
	rgb[2] = std::min(std::max((double) rgb[2], 0.0), 1.0);

	return rgb;
}

std::vector<float> msh2rgb(const std::vector<float>& msh)
{
	auto rgb = lab2rgb(msh2lab(msh));
	rgb[0] = std::min(std::max((double) rgb[0], 0.0), 1.0);
	rgb[1] = std::min(std::max((double) rgb[1], 0.0), 1.0);
	rgb[2] = std::min(std::max((double) rgb[2], 0.0), 1.0);
	return rgb;
	//return lab2rgb(msh2lab(msh));
}

int loadColorMapNames(std::string file, std::vector<std::string>& names)
{
	names.clear();

	std::string ext = file.substr(file.find_last_of(".") + 1);
	if (ext == "xml")
	{
		pugi::xml_document doc;
		if (int io = loadXml(file, doc) != 0)
			return io;

		std::string xquery;
		pugi::xpath_node x0;
		try
		{
			xquery = "/ColorMaps/ColorMap";
			x0 = doc.select_node(xquery.c_str());
			if (!x0) throw xPathException;
			for (pugi::xml_node x = x0.node(); x; x = x.next_sibling("ColorMap"))
				names.push_back(x.attribute("name").value());
		}
		catch (const std::exception& e)
		{
			std::cout << "\nError:  cannot parse colormap file \"" << file << "\"." << std::endl;
			std::cout << e.what() << std::endl;
			return ERR_XML_PARSE;
		}
	}
	else if (ext == "json")
	{
		std::ifstream ifs(file);
		json inj;

		try
		{
			ifs >> inj;
			ifs.close();  // close to release lock on file
		}
		catch (const std::exception& e)
		{
			std::cout << "\nError:  cannot load JSON colormap file \"" << file << "\"." << std::endl;
			std::cout << e.what() << std::endl;
			return ERR_JSON;
		}

		//std::cout << "Input JSON =\n" << std::setw(4) << inj << "\n" << std::endl;

		unsigned int i = 0;
		while (i < inj.size())
		{
			names.push_back(inj[i]["Name"]);
			i++;
		}
	}
	else
	{
		std::cout << "\nError:  unknown colormap file extension \"." << ext << "\"." << std::endl;
		return ERR_FILETYPE;
	}
	return 0;
}

class ColorMap
{
	public:

	std::vector<double> x, r, g, b;
	std::string space;
	bool inv, paraView;
	double nanr, nang, nanb;
	int imap;

	int load(std::string file, std::string name)
	{
		x.clear();
		r.clear();
		g.clear();
		b.clear();
		nanr = 0;
		nang = 0;
		nanb = 0;
		space = "RGB";

		std::string ext = file.substr(file.find_last_of(".") + 1);
		if (ext == "xml")
		{
			pugi::xml_document doc;
			if (int io = loadXml(file, doc) != 0)
				return io;

			std::string xquery;
			pugi::xpath_node point0;
			try
			{
				xquery = "/ColorMaps/ColorMap[@name='" + name + "']/Point";
				point0 = doc.select_node(xquery.c_str());
				if (!point0) throw xPathException;
				//std::cout << "point0 = " << point0.node().first_attribute().value() << "\n";
				space = point0.node().parent().attribute("space").value();
				for (pugi::xml_node point = point0.node(); point; point = point.next_sibling("Point"))
				{
					//std::cout << "b = " << point.attribute("b").value() << std::endl;
					x.push_back(std::stod(point.attribute("x").value()));
					r.push_back(std::stod(point.attribute("r").value()));
					g.push_back(std::stod(point.attribute("g").value()));
					b.push_back(std::stod(point.attribute("b").value()));
				}
			}
			catch (const std::exception& e)
			{
				std::cout << "\nError:  cannot parse colormap \"" << name << "\"." << std::endl;
				std::cout << e.what() << std::endl;
				return ERR_XML_PARSE;
			}

			try
			{
				xquery = "/ColorMaps/ColorMap[@name='" + name + "']/NaN";
				point0 = doc.select_node(xquery.c_str());
				if (!point0) throw xPathException;
				pugi::xml_node nan = point0.node();
				nanr = std::stod(nan.attribute("r").value());
				nang = std::stod(nan.attribute("g").value());
				nanb = std::stod(nan.attribute("b").value());
			}
			catch (const std::exception& e)
			{
				// Do nothing.  Use default NaN color.
			}
		}
		else if (ext == "json")
		{
			std::ifstream ifs(file);
			json inj;

			try
			{
				ifs >> inj;
				ifs.close();  // close to release lock on file
			}
			catch (const std::exception& e)
			{
				std::cout << "\nError:  cannot load JSON colormap file \"" << file << "\"." << std::endl;
				std::cout << e.what() << std::endl;
				return ERR_JSON;
			}

			//std::cout << "Input JSON =\n" << std::setw(4) << inj << "\n" << std::endl;

			unsigned int i = 0;
			while (inj[i]["Name"] != name)
			{
				i++;
				if (i >= inj.size())
				{
					std::cout << "\nError:  cannot find colormap name \"" << name << "\"." << std::endl;
					return ERR_JSON;
				}
			}

			//std::string name0 = inj[i]["Name"];
			//std::cout << "name0 = " << name0 << std::endl;

			std::vector<float> xrgb = inj[i]["RGBPoints"];
			//std::cout << "xrgb = " << xrgb << std::endl;
			for (unsigned int j = 0; j < xrgb.size();)
			{
				x.push_back(xrgb[j++]);
				r.push_back(xrgb[j++]);
				g.push_back(xrgb[j++]);
				b.push_back(xrgb[j++]);
			}

			space = inj[i]["ColorSpace"];
			try
			{
				std::vector<float> nan = inj[i]["NanColor"];
				nanr = nan[0];
				nang = nan[1];
				nanb = nan[2];
			}
			catch (const std::exception& e)
			{
				// Do nothing.  Use default NaN color.
			}
		}
		else
		{
			std::cout << "\nError:  unknown colormap file extension \"." << ext << "\"." << std::endl;
			return ERR_FILETYPE;
		}

		std::cout << "NaN RGB = " << nanr << nang << nanb << "\n";
		std::cout << "Colormap space = " << space << "\n";
		std::cout << "Number of interpolation points = " << x.size() << "\n";
		std::cout << std::endl;

		// Sort by x.
		auto idxc = sortidx(x);
		std::vector<double> tmp;
		tmp = x;
		for (int i = 0; i < x.size(); i++)
			x[i] = tmp[idxc[i]];
		tmp = r;
		for (int i = 0; i < r.size(); i++)
			r[i] = tmp[idxc[i]];
		tmp = g;
		for (int i = 0; i < g.size(); i++)
			g[i] = tmp[idxc[i]];
		tmp = b;
		for (int i = 0; i < b.size(); i++)
			b[i] = tmp[idxc[i]];

		// Map x to range [0, 1], since some XMLs use [-1, 1] instead.
		auto cxmin = x[0];
		auto cxmax = x[x.size() - 1];
		for (int i = 0; i < x.size(); i++)
			x[i] = (x[i] - cxmin) / (cxmax - cxmin);

		//std::cout << "Color map:\n";
		//for (int i = 0; i < x.size(); i++)
		//{
		//	std::cout << "x = " << x[i] << ", r = " << r[i]
		//	        << ", g = " << g[i] << ", b = " << b[i] << "\n";
		//}

		return 0;
	}

	std::vector<uint8_t> map(double y) const
	{
		// "const" after the function name here means that this
		// method does not modify the class.

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
		if (paraView)
		{
			auto iter = std::lower_bound(x.begin(), x.end(), y);  // binary search
			int i = std::distance(x.begin(), iter);
			//std::cout << "i = " << i << "\n";
			//std::cout << "y = " << y << "\n";
			//std::cout << "iter = " << *iter << "\n";
			if (i > 0)
			{
				double f = (y - x[i-1]) / (x[i] - x[i-1]);
				if (space == "Diverging")
				{
					// TODO:  this may not be 100% correct, it's hard to
					// tell.  The adjustHue function described in the
					// reference may be required.  But there's only one
					// colormap in my XML files that uses this space
					// anyway.

					std::vector<float> rgb0(3);
					rgb0[0] = r[i-1];
					rgb0[1] = g[i-1];
					rgb0[2] = b[i-1];

					std::vector<float> rgb1(3);
					rgb1[0] = r[i];
					rgb1[1] = g[i];
					rgb1[2] = b[i];

					auto msh0 = rgb2msh(rgb0);
					auto msh1 = rgb2msh(rgb1);

					// Interpolating between an angle close to
					// +pi and an angle close to -pi should not
					// give a result near zero.
					if (abs(msh1[2] - msh0[2]) > pi)
					{
						if (msh0[2] < msh1[2])
							msh0[2] += 2 * pi;
						else
							msh1[2] += 2 * pi;
					}

					std::vector<float> msh(3);
					msh[0] = (1.0 - f) * msh0[0] + f * msh1[0];
					msh[1] = (1.0 - f) * msh0[1] + f * msh1[1];
					msh[2] = (1.0 - f) * msh0[2] + f * msh1[2];

					auto rgbf = msh2rgb(msh);

					rgb[0] = twofivesix * rgbf[0];
					rgb[1] = twofivesix * rgbf[1];
					rgb[2] = twofivesix * rgbf[2];

					////auto xyz0 = rgb2xyz(rgb0);
					////auto xyz  = msh2xyz(msh);
					////auto xyz1 = rgb2xyz(rgb1);
					//std::cout << "rgb0 = " << rgb0 << std::endl;
					//std::cout << "rgb  = " << (int) rgb[0] << " " << (int) rgb[1] << " " << (int) rgb[2] << std::endl;
					//std::cout << "rgbf = " << rgbf << std::endl;
					//std::cout << "rgb1 = " << rgb1 << std::endl;
					//std::cout << "msh0 = " << msh0 << std::endl;
					//std::cout << "msh  = " << msh  << std::endl;
					//std::cout << "msh1 = " << msh1 << std::endl;
					////std::cout << "xyz0 = " << xyz0 << std::endl;
					////std::cout << "xyz  = " << xyz  << std::endl;
					////std::cout << "xyz1 = " << xyz1 << std::endl;
					//std::cout << std::endl;
				}
				else if (space == "Lab")
				{
					std::vector<float> rgb0(3);
					rgb0[0] = r[i-1];
					rgb0[1] = g[i-1];
					rgb0[2] = b[i-1];

					std::vector<float> rgb1(3);
					rgb1[0] = r[i];
					rgb1[1] = g[i];
					rgb1[2] = b[i];

					auto lab0 = rgb2lab(rgb0);
					auto lab1 = rgb2lab(rgb1);

					std::vector<float> lab(3);
					lab[0] = (1.0 - f) * lab0[0] + f * lab1[0];
					lab[1] = (1.0 - f) * lab0[1] + f * lab1[1];
					lab[2] = (1.0 - f) * lab0[2] + f * lab1[2];

					auto rgbf = lab2rgb(lab);

					rgb[0] = twofivesix * rgbf[0];
					rgb[1] = twofivesix * rgbf[1];
					rgb[2] = twofivesix * rgbf[2];

					//auto xyz0 = rgb2xyz(rgb0);
					//auto xyz  = lab2xyz(lab);
					//auto xyz1 = rgb2xyz(rgb1);
					//std::cout << "rgb0 = " << rgb0 << std::endl;
					//std::cout << "rgb  = " << (int) rgb[0] << " " << (int) rgb[1] << " " << (int) rgb[2] << std::endl;
					//std::cout << "rgbf = " << rgbf << std::endl;
					//std::cout << "rgb1 = " << rgb1 << std::endl;
					//std::cout << "lab0 = " << lab0 << std::endl;
					//std::cout << "lab  = " << lab  << std::endl;
					//std::cout << "lab1 = " << lab1 << std::endl;
					//std::cout << "xyz0 = " << xyz0 << std::endl;
					//std::cout << "xyz  = " << xyz  << std::endl;
					//std::cout << "xyz1 = " << xyz1 << std::endl;
					//std::cout << std::endl;
				}
				else if (space == "HSV")
				{
					std::vector<float> rgb0(3);
					rgb0[0] = r[i-1];
					rgb0[1] = g[i-1];
					rgb0[2] = b[i-1];

					std::vector<float> rgb1(3);
					rgb1[0] = r[i];
					rgb1[1] = g[i];
					rgb1[2] = b[i];

					auto hsv0 = rgb2hsv(rgb0);
					auto hsv1 = rgb2hsv(rgb1);

					std::vector<float> hsv(3);
					hsv[0] = (1.0 - f) * hsv0[0] + f * hsv1[0];
					hsv[1] = (1.0 - f) * hsv0[1] + f * hsv1[1];
					hsv[2] = (1.0 - f) * hsv0[2] + f * hsv1[2];

					auto rgbf = hsv2rgb(hsv);

					rgb[0] = twofivesix * rgbf[0];
					rgb[1] = twofivesix * rgbf[1];
					rgb[2] = twofivesix * rgbf[2];
				}
				else  // RGB
				{
					// TODO:  add a warning for unrecognized color spaces.
					rgb[0] = twofivesix * ((1.0 - f) * r[i-1] + f * r[i]);
					rgb[1] = twofivesix * ((1.0 - f) * g[i-1] + f * g[i]);
					rgb[2] = twofivesix * ((1.0 - f) * b[i-1] + f * b[i]);
				}
			}
			else
			{
				rgb[0] = twofivesix * r[i];
				rgb[1] = twofivesix * g[i];
				rgb[2] = twofivesix * b[i];
			}
		}
		else if (0 <= imap && imap < 3)
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
};

#endif  // INCLUDE_COLORMAP_H_

