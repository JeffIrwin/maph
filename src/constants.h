
#ifndef INCLUDE_CONSTANTS_H_
#define INCLUDE_CONSTANTS_H_

#include <colormapper_constants.h>

#if defined(_WIN32) || defined(_WIN64)
	const std::string slash = "\\";
#else
	const std::string slash = "/";
#endif

// Error return codes
const int ERR_CMD_ARG   = -1;
const int ERR_XML_OPEN  = -3;
const int ERR_PNG       = -5;

#endif  // INCLUDE_CONSTANTS_H_

