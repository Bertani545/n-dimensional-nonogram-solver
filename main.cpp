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
		hints[0] gives {2, 0}
		hints[1] gives {1 ,1}

	*/
	vector<vector<vector<int>>> hints;
	vector<int> data;
	vector<int> unsolvedData;
	vector<int> strides;
	vector<vector<int>> hintStrides;

	// For dimensions: x, y, z, ...
	void buildStrides() {
		this->strides = vector<int>(this->numberDimensions);
		this->strides[0] = 1;
		for (int i = 1; i < this->numberDimensions; i++) {
			this->strides[i] = this->strides[i - 1] * this->dimensionsSize[i - 1];
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

		
	}

	void buildHintStrides() {
		// hintStrides[d][k] = stride contribution of dimension k when building hint index for dim d
		this->hintStrides = vector<vector<int>>(this->numberDimensions, vector<int>(this->numberDimensions, 0));
		for (int d = 0; d < numberDimensions; d++) {
		    int factor = 1;
		    for (int k = 0; k < numberDimensions; k++) {
		        if (k == d) continue;
		        this->hintStrides[d][k] = factor;
		        factor *= dimensionsSize[k];
		    }
		}
	}


	

	int getHintIndex(int dim, const vector<int>& idx) {
	    int h = 0;
	    for (int i = 0; i < numberDimensions; i++) {
	        if (i == dim) continue;
	        h += idx[i] * hintStrides[dim][i];
	    }
	    return h;
	}

	void writeLine(int dim, vector<int>& idx, vector<int>& newValues) {
		int length = this->dimensionsSize[dim];
		
		// Get offset
		int offset = 0;
		int ri = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (i == dim) continue;
			if (idx[ri] >= this->dimensionsSize[i] || idx[ri] < 0) throw std::invalid_argument("idx content not valid for getLine");
			offset += idx[ri] * this->strides[i];
			ri++;
		}

		for (int i = 0; i < length; i++)
			this->unsolvedData[offset + i * this->strides[dim]] = newValues[i];

	}


	// Solves as much as possible from the line
	// returns false if nothing was placed or multiple solutions exists
	bool line_solver(int dim, vector<int>& idx) {
		vector<int> line = this->getLine(dim, idx, false);
		vector<int> lineHints = this->getHintLine(dim, idx);

		int n = line.size();
		int k = lineHints.size();

		// Leftmost-start
		vector<int> left(k);
		int pos = 0;
		for (int i = 0; i < k; i++) {
			while (pos < n && line[pos] == 0) pos++; // Skip empty cells
			left[i] = pos;
			pos += lineHints[i];
			if (pos < n) pos++;
		}

		// Rightmost
		vector<int> right(k);
		pos = n;
		for (int i = k-1; i >= 0; i--) {
			pos -= lineHints[i];
			while(pos > 0 && line[pos + lineHints[i] - 1] == 0) pos--;
			right[i] = pos;
			if (pos > 0) pos--;
		}


		// Overlap
		vector<int> newLine = line;
		for (int i = 0; i < k; i++) {

			int startMax = right[i];
			int endMin = left[i] + lineHints[i] - 1;


			while (endMin >= startMax) {
				newLine[endMin] = 1;
				endMin--;
			}
		}

		// No forced empties for the moment for i don't know how to do it myslef xd

		for (int pos = 0; pos < n; pos++) {
			bool canBeFilled = false;
			for (int i = 0; i < k; i++) {
				int L = left[i];
				int R = right[i];
				int startMin = L;
				int startMax = R;
				if (pos >= startMin && pos <= startMax + lineHints[i] - 1) {
					canBeFilled = true;
					break;
				}

			}
			if (!canBeFilled && newLine[pos] == -1) newLine[pos] = 0;
		}

		

		writeLine(dim, idx, newLine);
		return true;
	}

	float getHeuristic(int dim,  vector<int>& idx) {
		vector<int> line = this->getLine(dim, idx, false);
		vector<int> lineHints = this->getHintLine(dim, idx);

		float he = 0;
		// Many hints
		for (int e : lineHints) {
			he += ((float) e )/ ((float) this->dimensionsSize[dim]);
		}

		// Many solved squares
		for (int e : line) {
			if (e != 1) continue;
			he += ((float) e )/ ((float) this->dimensionsSize[dim]);
		}
		return he;
	}

	// Will try to solve the nonogram from the lists
	// if impossible or it does not equal the original data, returns false
	// else, it saves the data (if it didn't exists) and returns true
	bool solve() { 
		this->unsolvedData = vector<int>(this->data.size(), -1);

		// Create all the possible lines
		// Solve the obvious ones
		// insert the others in a priority queue
		for (int dim = 0; dim < this->numberDimensions; dim++) {

		}

		// Iterate until empty


		return true;
	}

	// Checks if the lists fit in the dimensions
	bool validate() {
		for (int i = 0; i < this->numberDimensions; i++) {
			for (vector<int> hint : this->hints[i]) {
				int total = hint.size() - 1;
				for (int e : hint) total += e;

				if (total > this->dimensionsSize[i]) return false;
			}
		}
		return true;
	}

	// Cuts the solved edges and re makes the lists
	void cut_edges() {


	}

	void line_density() {

	}

public:
	Nonogram(){};
	~Nonogram(){};

	// In which direction dim would you like to solve and which of all possible lines
	vector<int> getHintLine(int dim, const vector<int>& idx) {
		if (idx.size() != this->numberDimensions - 1) {
			throw std::invalid_argument("idx size must be numberDimensions - 1");
		}

		if (dim >= this->numberDimensions) {
			throw std::invalid_argument("No valid dimension");
		}

		// Add trash to dim position
	    int ri = 0; // index into reducedIdx
	    vector<int> new_idx(this->numberDimensions);
	    for (int j = 0; j < this->numberDimensions; j++) {
	        if (j == dim) continue;
	        if (idx[ri] < 0 || idx[ri] >= this->dimensionsSize[ri]) throw std::invalid_argument("idx content not valid for getHintLine");
	        new_idx[j] = idx[ri];
	        ri++;
	    }

	    int h = getHintIndex(dim, new_idx);
	    return hints[dim][h];
	}

	int getValue(vector<int> idx) {
		if (idx.size() != this->numberDimensions) return -1;

		int dataPos = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (idx[i] < 0 || idx[i] >= this->dimensionsSize[i]) throw std::invalid_argument("idx content not valid for getValue");
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
		clean();
		ifstream file(pathToPuzzle);
		if (!file) return false;
		if (!(file >> this->numberDimensions)) return clean();
		
		// Parse dimension size
		this->dimensionsSize = vector<int>(this->numberDimensions);

		
		int total = 1;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (!(file >> this->dimensionsSize[i])) return clean();
			total *= dimensionsSize[i];
		}
		cout << this->numberDimensions << endl;
		for (int dim : this->dimensionsSize) cout << dim << " ";
		cout << endl << endl;

		// For idx calculations
		this->buildStrides();

		// Parse data
		this->data = vector<int>(total);
		for (int i = 0; i < total; i++)
			if (!(file >> this->data[i])) return clean();

		this->buildLists();
		this->buildHintStrides();

		return solve();
	}

	bool buildFromLists(string pathToLists) {
		clean();
		ifstream file(pathToLists);
		if (!file) return false;

		string line;
		getline(file, line);
		istringstream firstLine(line);

		if (!(firstLine >> this->numberDimensions)) return clean();
		
		// Parse dimension size
		long long total_hints = 1;
		this->dimensionsSize = vector<int>(this->numberDimensions);
		for (int i = 0; i < this->numberDimensions; i++) {
			if (!(firstLine >> this->dimensionsSize[i])) return clean();
			total_hints *= this->dimensionsSize[i];
		}
		cout << this->numberDimensions << endl;
		for (int dim : this->dimensionsSize) cout << dim << " ";
		cout <<  "Total hints: "<< total_hints << endl;
		cout << endl;

		this->buildStrides();
		
		// Parse hints per dimension
		this->hints = vector<vector<vector<int>>>(this->numberDimensions);
		for (int d = 0; d < this->numberDimensions; d++) {
			int count = total_hints / this->dimensionsSize[d];
			this->hints[d] = vector<vector<int>>(count);
			for (int i = 0; i < count; i++) {
				if (!getline(file, line)) return clean();
				istringstream iss(line);
				vector<int> currHints;
				int val;
				while(iss >> val) currHints.push_back(val);
				this->hints[d][i] = currHints;
			}

			
		}

		this->buildHintStrides();

		//if (!this->validate()) return clean();

		return this->validate() ? solve() : clean();		
		//return solve();
	}

	// Getters
	int getDimensions() { return numberDimensions;}

	// Copy
	vector<int> getData() {return data;}
	int getDimensionSize(int i) { 
		if (i >= this->dimensionsSize.size()) return -1;
		return dimensionsSize[i];
	}
	vector<vector<int>> getHintsToSolveDimension(int i) {
		if (i >= this->dimensionsSize.size()) return vector<vector<int>>();
		return hints[i];
	}
	vector<int> getHintsDimensionRow(int dim, int r) {
		if (dim >= this->dimensionsSize.size()) return vector<int>();
		if (r >= this->hints[dim].size()) return  vector<int>();
		return hints[dim][r];
	}

	vector<int> getLine(int dim, vector<int>& idx, bool solved = true) {
		int length = this->dimensionsSize[dim];
		vector<int> line;
		
		// Get offset
		int offset = 0;
		int ri = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (i == dim) continue;
			if (idx[ri] >= this->dimensionsSize[i] || idx[ri] < 0) throw std::invalid_argument("idx content not valid for getLine");
			offset += idx[ri] * this->strides[i];
			ri++;
		}
		if (solved)
			for (int i = 0; i < length; i++) 
				line.push_back(this->data[offset + i * this->strides[dim]]);				
		else
			for (int i = 0; i < length; i++)
				line.push_back(this->unsolvedData[offset + i * this->strides[dim]]);

		return line;
	}

	vector<int> solveLine(int dim, vector<int>& idx) {
		if (line_solver(dim, idx))
		{
			return this->getLine(dim, idx, false);
		}
		else {
			cout << "No square was placed" << endl;
		}

		return {};
	}

};




int main(int argc, char const *argv[])
{
	if (argc < 2) return 1;
	Nonogram nonogram;
	string path;
	
	cout << "From grid:" << endl;
	path = argv[1];
	nonogram.buildFromSquares(path);
	for (int i = 0; i < nonogram.getDimensions(); i++) {
		vector<vector<int>> temp1 = nonogram.getHintsToSolveDimension(i);
		for (vector<int> v : temp1) {
			for (int e : v) cout << e << " ";
				cout << endl;
		}cout << endl;
	}

	vector<int> idx = {0};
	vector<int> line = nonogram.getLine(0, idx);
	vector<int> hintLine =  nonogram.getHintLine(0, idx);
	for (int e : line)
		cout << e << " ";
	cout << endl;
	for (int e : hintLine)
		cout << e << " ";
	cout << endl;

	vector<int> solvedLine = nonogram.solveLine(0, idx);
	cout << "Filled squares :" << endl;
	for (int e : solvedLine)
		cout << e << " ";
	cout << endl;

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

	

/*
	cout << endl;
	vector<int> temp = nonogram.getHintLine({1}, 0);
	for (int e : temp) cout << e << " ";
		cout << endl;
		*/
	return 0;

}