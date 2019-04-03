#include <iostream>
#include <coin/CoinUtilsConfig.h>
#include <coin/CbcConfig.h>

int main(int argc, char** argv) {
  std::cout << "coinutils version: " << COINUTILS_VERSION << std::endl;
  std::cout << "cbc version: " << CBC_VERSION << std::endl;
	return 0;
}

