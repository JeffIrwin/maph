# maph
*This isn't normal, but on maph it is*

Make a heatmap of everywhere you've worked out from Strava.

## Compile

Use CMake, or just run `build.sh` on Linux or `build.cmd` on Windows.

## Download your bulk GPX data from Strava

For instructions, see https://support.strava.com/hc/en-us/articles/216918437-Exporting-your-Data-and-Bulk-Export#Bulk

Be careful!  This contains sensitive personal data.

## Configure your JSON input for maph

Choose your colormap file and colormap name.  Otherwise, the default colormap is a full-spectrum red-to-red HSV rainbow.

Set a globbed path to your GPX files.  If you don't have data, some stripped-down samples are included in [data/export\_0/](data/export_0/).

By default, maph will snap to the bounds of your data.  If you have worked out in multiple cities, this won't be very helpful.  To zoom in, set lattitudes (y) and longitudes (x) with `Min x`, ..., `Max y`.  Google maps can be helpful for finding the coordinates that you want.

Here is an example:

    {
    	"Image size x"      : 3840,
    	"Image size y"      : 2160,
    
    	"Colormap file"     : "submodules/colormaps/ColorMaps5.0.0.xml",
    	"Colormap name"     : "Black-Body Radiation",
    
    	"GPX files"         : "data/export_0/activities/*.gpx",
    	"Min x"             : -80.95,
    	"Max x"             : -80.857,
    	"Min y"             : -14.414,
    	"Max y"             : -14.30,
    
    	"Verbosity"         : 0,
    	"File name prefix"  : "example-2a"
    }

## Run

    ./target/maph data/input_1.json

Maph will write the results to a PNG file.  A couple examples are below (we did the maph):

![Viridis colormap](https://raw.githubusercontent.com/JeffIrwin/maph/master/data/expected-output/example-1a.png)

![Black-body radiation colormap](https://raw.githubusercontent.com/JeffIrwin/maph/master/data/expected-output/example-2a.png)

![Blue to red rainbow colormap](https://raw.githubusercontent.com/JeffIrwin/maph/master/data/expected-output/example-3a.png)

## Platforms

This tool has been tested on:
- MSVC 19.0.24215.1 (Visual Studio 14 2015)
- MSVC 19.15.26726.0 (Visual Studio 15 2017)
- MSVC 19.24.28314.0 (Visual Studio 16 2019)
- Ubuntu 18.04.3 LTS (Bionic Beaver) with gcc (Ubuntu 7.4.0-1ubuntu1\~18.04.1) 7.4.0

Other platforms and compilers may work.

