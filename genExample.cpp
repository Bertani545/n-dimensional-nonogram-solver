#include<iostream>
#include"nonogram.hpp"


int main() {
	Nonogram nonogram;

	vector<int> dimSizes = {5,5,5,5,5};
	nonogram.buildRandom(5, dimSizes, 0.6);
	nonogram.printSolved();
	nonogram.printHints();
}