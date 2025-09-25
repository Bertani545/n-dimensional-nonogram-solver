#include<iostream>
#include"nonogram.hpp"



int main(int argc, char const *argv[])
{
	if (argc < 2) return 1;
	Nonogram nonogram;
	string path;
	
	cout << "From grid:" << endl;
	path = argv[1];
	nonogram.buildFromSquares(path);
	nonogram.printOriginal();
	nonogram.printSolved();


	cout << "From lists:" << endl;
	path = argv[2];
	if (nonogram.buildFromLists(path)) {
		for (int i = 0; i < nonogram.getDimensions(); i++) {
			vector<vector<int>> temp1 = nonogram.getHintsToSolveDimension(i);
			for (vector<int> v : temp1) {
				for (int e : v) cout << e << " ";
					cout << endl;
			}cout << endl;
		}
	} else {
		throw std::invalid_argument("Hints do not match the puzzle");
	}

	nonogram.printSolved();


	vector<int> dims = {5, 5};
	vector<int> data = {1,0,1,0,1,1,1,1,1,1,0,1,1,1,0,0,1,0,1,0,0,1,1,1,0};
	nonogram.buildFromData(2, dims, data);

	cout << "From Data:" << endl;
	nonogram.printOriginal();
	nonogram.printSolved();

	return 0;

}