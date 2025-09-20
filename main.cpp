#include<iostream>
#include <fstream>
#include<vector>
#include<string>
#include <sstream>

using namespace std;



class Nonogram
{
private:
	int numberDimensions;
	vector<int> dimensionsSize;

	// Hint to solve dimension dim, over all possible other dimensions
	/*
		Example, in a 2 x 2

		* *
		´ ´
		hints[0][1] gives {1} because there is 1 * on column 1 
		hints[1][0] gives {2} because there are 2 * on row 0

	*/
	vector<vector<vector<int>>> hints; // Dimension -> Row -> hint -> number
	vector<int> data;
	vector<int> strides;
	vector<vector<int>> hintStrides;

	void buildStrides() {
		this->strides = vector<int>(this->numberDimensions);
		this->strides[this->numberDimensions - 1] = 1;
		for (int i = this->numberDimensions - 2; i >= 0; i--) {
			this->strides[i] = this->strides[i + 1] * this->dimensionsSize[i + 1];
		}
	}


	vector<int> computeLine(vector<int>& line) {
		int n = line.size();

		int last_number = line[0];
		int segment_size = (last_number == 1) ? 1 : 0;
		int total_numbers = (last_number == 1) ? 1 : 0;
		vector<int> hints;
		for (int i = 1; i < n; i++) {
			if (line[i]) { // 1
				if (line[i] == last_number) {
					segment_size++;
				} else {
					// New block
					segment_size = 1;
					total_numbers++;
				}
			} else { // 0
				if (line[i] != last_number) { // 1 and then 0
					hints.push_back(segment_size);
				}
			}
			last_number = line[i];
		}
		if (last_number == 1) hints.push_back(segment_size); // Last correction
		if (hints.size() == 0) hints.push_back(0);

		return hints;
			
	}
	

	void buildLists() {

		long long totalMult = 1;
		for (int s : this->dimensionsSize) totalMult *= s;

		this->hints = vector<vector<vector<int>>>(this->numberDimensions);
		for (int d = 0; d < this->numberDimensions; d++) {

			int length = this->dimensionsSize[d];   // how long each line is
			int count  = totalMult / length;        // how many lines along d
			int step   = this->strides[d];          // stride along d
			this->hints[d] = vector<vector<int>>(count);

			for (int h = 0; h < count; h++) {
				vector<int> line;
				line.reserve(length);

				// compute starting offset in flat data
		        // "h" encodes the fixed coordinates of all other dims
		        int offset = 0;
		        int tmp = h;
		        for (int k = 0; k < this->numberDimensions; k++) {
		            if (k == d) continue;
		            int coord = tmp % this->dimensionsSize[k];
		            tmp /= this->dimensionsSize[k];
		            offset += coord * this->strides[k];
		        }

		        // collect values along dimension d
		        for (int i = 0; i < length; i++) {
		            line.push_back(this->data[offset + i * step]);
		        }


		        this->hints[d][h] = computeLine(line);

			}
			
		}

		/*
		for(int i = 0; i < dimensionsSize[0]; i++) {
			for (int j = 0; j < dimensionsSize[1]; j++) {
				cout << getValue({i,j}) << " ";
			}cout << endl;
		}
		*/

		// hintStrides[d][k] = stride contribution of dimension k when building hint index for dim d
		this->hintStrides = vector<vector<int>>(this->numberDimensions, vector<int>(this->numberDimensions, 0));
		for (int d = 0; d < numberDimensions; d++) {
		    int factor = 1;
		    for (int k = 0; k < numberDimensions; k++) {
		        if (k == d) continue;
		        hintStrides[d][k] = factor;
		        factor *= dimensionsSize[k];
		    }
		}
	}

	bool solve() { return false;}

	int getHintIndex(const vector<int>& idx, int dim) {
	    int h = 0;
	    for (int i = 0; i < numberDimensions; i++) {
	        if (i == dim) continue;
	        h += idx[i] * hintStrides[dim][i];
	    }
	    return h;
	}


	
public:
	Nonogram(){};
	~Nonogram(){};

	vector<int> getHintLine(const vector<int>& idx, int dim) {
		if (idx.size() != this->numberDimensions - 1) {
			throw std::invalid_argument("idx size must be numberDimensions - 1");
		}

		// Add trash to dim position
	    int ri = 0; // index into reducedIdx
	    vector<int> new_idx(this->numberDimensions);
	    for (int j = 0; j < this->numberDimensions; j++) {
	        if (j == dim) continue;
	        new_idx[j] = idx[ri];
	        ri++;
	    }

	    int h = getHintIndex(new_idx, dim);
	    return hints[dim][h];
	}

	int getValue(vector<int> idx) {
		if (idx.size() != this->numberDimensions) return -1;

		int dataPos = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			dataPos += idx[i] * this->strides[i];
		}
		if (dataPos >= this->data.size()) return -1;
		return this->data[dataPos];
	}

	bool clean() {
		numberDimensions = 0;
		dimensionsSize = vector<int>();
		hints = vector<vector<vector<int>>>();
		data = vector<int>();
		return false;
	}

	

	bool buildFromSquares(string pathToPuzzle) {
		ifstream file(pathToPuzzle);
		if (!file) return false;
		if (!(file >> this->numberDimensions)) return clean();
		
		// Parse dimension size
		this->dimensionsSize = vector<int>(this->numberDimensions);

		cout << this->numberDimensions << endl;
		int total = 1;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (!(file >> this->dimensionsSize[i])) return clean();
			total *= dimensionsSize[i];
		}


		// For idx calculations
		this->buildStrides();

		// Parse data
		this->data = vector<int>(total);
		for (int i = 0; i < total; i++)
			if (!(file >> this->data[i])) return clean();

		this->buildLists();

		return true;
	}

	bool buildFromLists(string pathToLists) {
		ifstream file(pathToLists);
		if (!file) return false;

		if (!(file >> this->numberDimensions)) return clean();
		
		// Parse dimension size
		this->dimensionsSize = vector<int>(this->numberDimensions);
		for (int i = 0; i < this->numberDimensions; i++) 
			if (!(file >> this->dimensionsSize[i])) return clean();

		this->buildStrides();
		
		// Parse hints per dimension
		this->hints = vector<vector<vector<int>>>(this->numberDimensions);
		string line;
		for (int d = 0; d < this->numberDimensions; d++) {
			int size = this->dimensionsSize[d];
			this->hints[d] = vector<vector<int>>(size);
			for (int i = 0; i < size; i++) {
				if (!getline(file, line)) return clean();
				istringstream iss(line);
				vector<int> currHints;
				int val;
				while(iss >> val) currHints.push_back(val);
				this->hints[d][i] = currHints;
			}

			
		}

		return true;
	}

	// Getters
	int getDimensions() { return numberDimensions;}

	// Copy
	vector<int> getData() {return data;}
	int getDimensionSize(int i) { 
		if (i >= this->dimensionsSize.size()) return -1;
		return dimensionsSize[i];
	}
	vector<vector<int>> getDimensionsHint(int i) {
		if (i >= this->dimensionsSize.size()) return vector<vector<int>>();
		return hints[i];
	}
	vector<int> getHintsDimensionRow(int dim, int r) {
		if (dim >= this->dimensionsSize.size()) return vector<int>();
		if (r >= this->hints[dim].size()) return  vector<int>();
		return hints[dim][r];
	}

	
};




int main(int argc, char const *argv[])
{
	if (argc < 2) return 1;
	string path = argv[1];
	Nonogram nonogram;
	nonogram.buildFromSquares(path);

	vector<vector<int>> temp1 = nonogram.getDimensionsHint(1);
	for (vector<int> v : temp1) {
		for (int e : v) cout << e << " ";
			cout << endl;
	}
/*
	vector<int> temp = nonogram.getHintLine({3}, 0);
	for (int e : temp) cout << e << " ";
		cout << endl;
	return 0;
*/
}