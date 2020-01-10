
#include <maph.h>

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

