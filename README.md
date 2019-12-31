# maph
*This isn't normal, but on maph it is*

## Compile

Use CMake, or just run `build.sh` on Linux.

## Download your bulk GPX data from Strava

For instructions, see https://support.strava.com/hc/en-us/articles/216918437-Exporting-your-Data-and-Bulk-Export#Bulk

Be careful!  This contains sensitive personal data.

## Configure your JSON input for maph

Choose your colormap file and colormap name.  Otherwise, the default colormap is a full-spectrum red-to-red HSV rainbow.

Set a globbed path to your GPX files.  If you don't have data, some stripped-down samples are included in data/export\_0/.

By default, maph will snap to the bounds of your data.  If you have worked out in multiple cities, this won't be very helpful.  To zoom in, set lattitudes and longitudes with `Min x`, ..., `Max y`.  Google maps can be helpful for finding the coordinates that you want.

## Run

    ./target/maph data/input_1.json

Maph will write the results to a PNG file.  A couple examples are below:

![Viridis colormap](https://raw.githubusercontent.com/JeffIrwin/maph/master/data/expected-output/example-1a.png)

![Black-body radiation colormap](https://raw.githubusercontent.com/JeffIrwin/maph/master/data/expected-output/example-2a.png)

![Blue to red rainbow colormap](https://raw.githubusercontent.com/JeffIrwin/maph/master/data/expected-output/example-3a.png)

