#ifndef NONOGRAM
#define NONOGRAM

#include<iostream>
#include <fstream>
#include<vector>
#include<string>
#include <sstream>
#include <queue>
#include<set>
#include<algorithm>
#include<random>
#include <chrono> 


using namespace std;

struct  LineDef {
	int dim;
	vector<int> idx;
};

enum class Color {
	Undefined = -1,
	White = 0,
	Black = 1,
};


class Nonogram
{
private:
	int numberDimensions;
	vector<int> dimensionsSize;

	// Hint to solve dimension dim, over all possible other dimensions
	/*
		Example, in a 2 x 2

		* *
		' '
		hints[0] gives {2, 0}
		hints[1] gives {1 ,1}

	*/
	vector<vector<vector<int>>> hints;
	vector<Color> data;
	vector<Color> unsolvedData;
	vector<int> strides;
	vector<vector<int>> hintStrides;
	int total_squares = 0;

	// For dimensions: x, y, z, ...
	void buildStrides() {
		this->strides = vector<int>(this->numberDimensions);
		this->strides[0] = 1;
		for (int i = 1; i < this->numberDimensions; i++) {
			this->strides[i] = this->strides[i - 1] * this->dimensionsSize[i - 1];
		}
	}


	vector<int> computeLine(vector<Color>& line) {
		int n = line.size();

		Color last_number = line[0];
		int segment_size = (last_number == Color::Black) ? 1 : 0;
		int total_numbers = (last_number == Color::Black) ? 1 : 0;
		vector<int> hints;
		for (int i = 1; i < n; i++) {
			if (line[i] == Color::Black) { // 1
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
		if (last_number == Color::Black) hints.push_back(segment_size); // Last correction
		if (hints.size() == 0) hints.push_back(0);

		return hints;
			
	}
	
	


	void buildLists() {
		long long totalMult = 1;
		for (int s : this->dimensionsSize) totalMult *= s;

		this->hints = vector<vector<vector<int>>>(this->numberDimensions);
		for (int d = 0; d < this->numberDimensions; d++) {

			int length = this->dimensionsSize[d];   // how long each line is
			int count  = totalMult / length;		// how many lines along d
			int step   = this->strides[d];		  // stride along d
			this->hints[d] = vector<vector<int>>(count);

			for (int h = 0; h < count; h++) {
				vector<Color> line;
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


	

	int getHintIndex(struct LineDef& def) {
		int h = 0;
		for (int i = 0; i < numberDimensions; i++) {
			if (i == def.dim) continue;
			h += def.idx[i] * hintStrides[def.dim][i];
		}
		return h;
	}

	void writeLine(struct LineDef& def, vector<Color>& newValues) {
		int length = this->dimensionsSize[def.dim];
		
		// Get offset
		int offset = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (i == def.dim) continue;
			if (def.idx[i] >= this->dimensionsSize[i] || def.idx[i] < 0) throw std::invalid_argument("idx content not valid for getLine");
			offset += def.idx[i] * this->strides[i];
		}

		for (int i = 0; i < length; i++)
			this->unsolvedData[offset + i * this->strides[def.dim]] = newValues[i];

	}

	// From "An Efficient Approach to Solving Nonograms" by Chen Wu et al.
	int Fix(int i, int j, vector<Color>& line, vector<int>& lineHints, vector<vector<int>>& dpTable) {
		if (dpTable[i][j] >= 0) return dpTable[i][j];

		// Base cases
		if (j == 0) {
			// If there are no more hints, the remaining cells must all be white (0)
			for (int k = 1; k <= i; ++k) {
				if (line[k] == Color::Black) {
					dpTable[i][j] = 0;
					return false;
				}
			}
			dpTable[i][j] = 1;
			return true;
		}
		if (i == 0) {
			// We have hints left but no cells. This is only valid if we have no hints (handled above).
			dpTable[i][j] = 0;
			return false;
		}

		//dpTable[i][j] = Fix0(i,j, line, lineHints, dpTable) || Fix1(i,j, line, lineHints, dpTable);
		bool res = false;

		// Option 1: cell i is white (0). Allowed if it's not explicitly black.
		if (line[i] != Color::Black) res = Fix(i - 1, j, line, lineHints, dpTable);

		// Option 2: place j-th block so it ends at i (block length h_j = lineHints[j-1])
		int h_j = lineHints[j-1];
		// Enough celss to fit
		if (!res && i >= h_j) {
			bool isBlock = true;
			for (int k = i - h_j + 1; k <= i; ++k) {
				if (line[k] == Color::White) { // Should be 1 or unknown, but if it's 0, it can't be a block.
					isBlock = false;
					break;
				}
			}

			if (isBlock) {
				int leftPos = i - h_j; // cell right before the block (0 if block starts at 1)
				if (leftPos == 0) { // The block starts at the beginning of the line
					res = res || Fix(0, j - 1, line, lineHints, dpTable);
				} else if (line[i - h_j] != Color::Black) { // separator cell (leftPos) must NOT be black (cannot be 1).
					res = res || Fix(i - h_j - 1, j - 1, line, lineHints, dpTable);
				}
			}
		}

		dpTable[i][j] = res;
		return res;
	}

	bool Coloring(int n, int k, vector<Color>& line, vector<int>& lineHints, vector<int>& modifiedEntries ) {
		bool changed = false;

		for (int i = 1; i <= n; ++i) {
			// If the cell is already determined, skip it
			if (line[i] != Color::Undefined) continue;

			// Try coloring the cell as '1' (black)
			line[i] = Color::Black;
			vector<vector<int>> dpTable1(n + 1, vector<int>(k + 1, -1));
			bool possible1 = Fix(n, k, line, lineHints, dpTable1);

			// Try coloring the cell as '0' (white)
			line[i] = Color::White;
			vector<vector<int>> dpTable0(n + 1, vector<int>(k + 1, -1));
			bool possible0 = Fix(n, k, line, lineHints, dpTable0);;

			// Apply the deduction rules
			if (possible1 && !possible0) {
				// It must be a '1' because '0' is not a valid solution
				line[i] = Color::Black;
				modifiedEntries.push_back(i-1);
				changed = true;
			} else if (!possible1 && possible0) {
				// It must be a '0' because '1' is not a valid solution
				line[i] = Color::White;
				modifiedEntries.push_back(i-1);
				changed = true;
			} else {
				line[i] = Color::Undefined;
			}
		}

		return changed;
	}

	bool lineSolver(struct LineDef& def, vector<int>& modifiedEntries) {
		vector<Color> line = getLine(def, false);
		line.insert(line.begin(), {Color::White}); 
		vector<int> lineHints = getHintLine(def);
/*
		cout << "Line: "<< endl;
		for (int e : line) cout << e  << " ";
			cout << endl;
		cout << "Hints:" << endl;
		for (int e : lineHints) cout << e  << " ";
			cout << endl;
*/
		int k = lineHints.size();
		int n = line.size() - 1;

		// DP method
		vector<vector<int>> dpTable(n +1 , vector<int>(k+1, -1)); 
		//if (!Fix(n - 1,k - 1, line, lineHints, dpTable)) return false;
		if (!Fix(n, k, line, lineHints, dpTable)) return false;

		bool changedInPass = true;
		bool totalChanged = false;
		while (changedInPass) {
			changedInPass = Coloring(n, k, line, lineHints, modifiedEntries);
			if (changedInPass) {
				totalChanged = true;
			}
		}

		if (totalChanged) {
			line.erase(line.begin());
			writeLine(def, line);
			return true;
		}
		//cout << "Nothig aout of this one" << endl;
		return false;
	}

	bool isTrivial (struct LineDef& def) {
		vector<int> lineHints = this->getHintLine(def);
		if (lineHints[0] == 0) return true;
		int total = lineHints.size() - 1;
		for (int hint : lineHints)
			total += hint;
		if (total == this->dimensionsSize[def.dim]) return true;
		return false;
	}

	// Assumed trivial fromt he beggining
	void solveTrivial(struct LineDef& def) {
		vector<int> lineHints = this->getHintLine(def);

		vector<Color> newLine;
		if (lineHints.size() == 1) {
			if (lineHints[0] == 0) 
				newLine = vector<Color>(this->dimensionsSize[def.dim], Color::White);
			else
				newLine = vector<Color>(this->dimensionsSize[def.dim], Color::Black);
		}
		else {
			newLine = vector<Color>(this->dimensionsSize[def.dim], Color::White);
			int i = 0;
			for (int hint : lineHints) {
				for (hint; hint > 0; hint--, i++) {
					newLine[i] = Color::Black;
				}
				i++;
			}
		}
		writeLine(def, newLine);
	}

	float getHeuristic(struct LineDef& def) {
		vector<Color> line = this->getLine(def, false);
		vector<int> lineHints = this->getHintLine(def);

		float he = 0;
		// Many hints
		for (int e : lineHints) {
			he += ((float) e )/ ((float) this->dimensionsSize[def.dim]);
		}

		// Many solved squares
		for (Color e : line) {
			if (e == Color::Undefined) continue;
			he += ((float) e )/ ((float) this->dimensionsSize[def.dim]);
		}
		return he;
	}

	vector<vector<int>> generateLineIndexTuples(int dim) {
		long long totalMult = 1;
		for (int s : this->dimensionsSize) totalMult *= s;
		
		int length = dimensionsSize[dim];	   // length of each line
		int count  = totalMult / length;		// how many lines along this dimension

		vector<vector<int>> tuples;
		tuples.reserve(count);

		for (int h = 0; h < count; h++) {
			vector<int> coords(dimensionsSize.size());
			int tmp = h;
			for (int k = 0; k < dimensionsSize.size(); k++) {
				if (k == dim) continue;
				coords[k] = tmp % dimensionsSize[k];
				tmp /= dimensionsSize[k];
			}
			// dim is "free" (0..length-1) when actually scanning line
			tuples.push_back(coords);
		}

		return tuples;
	}

	//using Item = std::pair<float, std::pair<int, std::vector<int>>>;
	using Item = std::pair<float, struct LineDef>;
	// We know a.second.second and b.second.second are equal in size but must differ somewhere
	struct Compare {
		bool operator()(const Item& a, const Item& b) const {
			if (a.first != b.first)
				return a.first < b.first; // bigger priority value = higher priority
			if (a.second.dim != b.second.dim)
				return a.second.dim < b.second.dim;
			int i = 0;
			while (a.second.idx[i] == b.second.idx[i]) i++;
			return a.second.idx[i] < b.second.idx[i];
		}
	};

	// Will try to solve the nonogram from the lists
	// if impossible or it does not equal the original data, returns false
	// else, it saves the data (if it didn't exists) and returns true
	bool solve() {
		int totalMult = 1;
		for (int s : this->dimensionsSize) totalMult *= s;

		this->unsolvedData = vector<Color>(totalMult, Color::Undefined);

		//priority_queue<Item, std::vector<Item>, Compare> nextAnalyze;
		set<Item, Compare> nextAnalyze; // To allow updates

		// Create all the possible lines
		// Solve the obvious ones
		// insert the others in a priority queue
		for (int dim = 0; dim < this->numberDimensions; dim++) {
			// All possible permutation for this dimension
			auto tuples = generateLineIndexTuples(dim);

			for (auto& coords : tuples) {
				LineDef def = {dim, coords};
				if (isTrivial(def)) {
					this->solveTrivial(def);
				} else {
					float heuristic = this->getHeuristic(def);
					Item toSolve = {heuristic, {dim, coords}};
					nextAnalyze.insert(toSolve);
				}
			}
		}
		


		vector<int> modifiedEntries;
		// Iterate until empty

		
		while (!nextAnalyze.empty()) {
			modifiedEntries.clear();
			const auto& top = *nextAnalyze.begin(); // Maybe
			auto [priority, inner] = top;
			struct LineDef old_line = inner;
			nextAnalyze.erase(nextAnalyze.begin());


 			
			bool solvedAnything = lineSolver(old_line, modifiedEntries);
			if (!solvedAnything) continue;
			//this->printSolved();

			// Build the next ones to try
			
			for (int idx : modifiedEntries) {
				vector<int> coords = old_line.idx;
				coords[old_line.dim] = idx;

				for (int d2 = 0; d2 < numberDimensions; d2++) {
					if (d2 == old_line.dim) continue;
					struct LineDef newLine = {d2, coords};
					float heuristic = this->getHeuristic(newLine);
					Item toSolve = {heuristic, {d2, coords}};

					// Chek for update
					auto it = std::find_if(nextAnalyze.begin(), nextAnalyze.end(), [&toSolve](const Item& x){
						return x.second.dim == toSolve.second.dim &&
								x.second.idx == toSolve.second.idx;
					});

					if (it != nextAnalyze.end()) {
						nextAnalyze.erase(it);
						nextAnalyze.insert(toSolve);
					} else {
						nextAnalyze.insert(toSolve);
					}

					
				}
			}
			
		}
		
		if (isSolved()) {
			this->data = vector<Color>(this->unsolvedData); // Consistent
			return true;
		}
		return false;
	}


	// Checks for remaining -1
	bool isSolved() {
		// Well, could be worse
		bool isSolved = true;
		for (Color e : this->unsolvedData) {
			if (e == Color::Undefined) {
				isSolved = false;
				break;
			}
			total_squares += (e == Color::Black); // We use this traverse to also check for this value
		}
		return isSolved;
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

	// Receives a modified line, constructs the hints with the mod
	// Brings back the ambiguity and checks if still valid
	//Only for building a random one
	bool lineChecker(struct LineDef& def) {
		// Modify the hint
		vector<Color> line = getLine(def, true);
		line.insert(line.begin(), {Color::White}); 
		int hintIdx = getHintIndex(def);
		this->hints[def.dim][hintIdx] = computeLine(line);
		vector<int> lineHints = getHintLine(def);


		line[def.dim+1] = Color::Undefined; // Insert the ambiguity
		int k = lineHints.size();
		int n = line.size() - 1;

		// DP method
		vector<vector<int>> dpTable(n +1 , vector<int>(k+1, -1)); 
		//if (!Fix(n - 1,k - 1, line, lineHints, dpTable)) return false;
		if (!Fix(n, k, line, lineHints, dpTable)) return false;

		bool changedInPass = true;
		bool totalChanged = false;
		vector<int> dummy;
		while (changedInPass) {
			changedInPass = Coloring(n, k, line, lineHints, dummy);
			if (changedInPass) {
				totalChanged = true;
			}
		}

		return totalChanged;
	}


			

	bool checkCellConsistency(int currIdx) {

		 // Decode currIdx into coordinates
	    vector<int> coords(this->numberDimensions);
	    for (int d = this->numberDimensions - 1; d >= 0; --d) {
	        coords[d] = currIdx / this->strides[d];
	        currIdx %= this->strides[d];
	    }

	    // For each dimension, check the line containing this cell
		for (int dim = 0; dim < this->numberDimensions; ++dim) {
	        vector<int> idx = coords; // cords[dim] is variable
	        struct LineDef def = {dim, idx};
	        if (!lineChecker(def)) {
	            return false; // line unsolvable
	        }
	    }

	    return true;
	}

public:
	Nonogram(){};
	~Nonogram(){};

	// In which direction dim would you like to solve and which of all possible lines
	vector<int> getHintLine(struct LineDef& def) {
		if (def.idx.size() != this->numberDimensions) {
			throw std::invalid_argument("idx size must be numberDimensions");
		}

		if (def.dim >= this->numberDimensions) {
			throw std::invalid_argument("No valid dimension");
		}


		int h = getHintIndex(def);
		return hints[def.dim][h];
	}

	int getCellId(vector<int>& idx) {
		int dataPos = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (idx[i] < 0 || idx[i] >= this->dimensionsSize[i]) throw std::invalid_argument("idx content not valid for getCellId");
			dataPos += idx[i] * this->strides[i];
		}
		if (dataPos >= this->data.size()) return -1;
		return dataPos;
	}

	Color getValue(vector<int>& idx) {
		if (idx.size() != this->numberDimensions) return Color::Undefined;
		int dataPos = getCellId(idx);
		return this->data[dataPos];
	}

	bool clean() {
		numberDimensions = 0;
		total_squares = 0;
		dimensionsSize.clear();
		hints.clear();
		data.clear();
		unsolvedData.clear();
		strides.clear();
		hintStrides.clear();

		return false;
	}

	
	bool buildFromData(int nDim, vector<int>& dimSizes, vector<int>& new_data) {
		clean();
		this->numberDimensions = nDim;
		if (dimSizes.size() != nDim) return clean();

		this->dimensionsSize = vector<int>(dimSizes);
		int total = 1;
		for (int d : this->dimensionsSize) total *= d;

		if (new_data.size() != total) return clean();
		this->buildStrides();
		this->data = vector<Color>(new_data.size());
		for (int i = 0; i < data.size(); ++i) { 
			if (new_data[i] < -1 || new_data[i] > 1) return clean();
			this->data[i] = static_cast<Color>(new_data[i]);
		}
		
		this->buildLists();
		this->buildHintStrides();
		return solve();
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
		/*
		cout << this->numberDimensions << endl;
		for (int dim : this->dimensionsSize) cout << dim << " ";
		cout << endl << endl;
		*/

		// For idx calculations
		this->buildStrides();

		// Parse data
		this->data = vector<Color>(total);
		int temp;
		for (int i = 0; i < total; i++){
			if (!(file >> temp)) return clean();
			if (temp < -1 || temp > 1) return clean();
			this->data[i] = static_cast<Color>(temp);
		}

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
		/*
		cout << this->numberDimensions << endl;
		for (int dim : this->dimensionsSize) cout << dim << " ";
		cout <<  "Total hints: "<< total_hints << endl;
		cout << endl;
		*/

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


	

	void buildRandom(int nDim, vector<int>& dimSizes, float appproxCoverPercent) {
		if (appproxCoverPercent > 1.0) return;

		clean();
		this->numberDimensions = nDim;
		if (dimSizes.size() != nDim) {clean(); return;}

		this->dimensionsSize = vector<int>(dimSizes);
		int total = 1;
		for (int d : this->dimensionsSize) total *= d;

		this->buildStrides();
		this->buildHintStrides();
		
		do {
			this->data = vector<Color>(total, Color::White); // Start at White
			double currCover = 0.0;
			double totalPlaced = 0;
			double floatTotal = (float) total;
			int currIdx = -1;

			auto const seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_int_distribution<int> distribution(0, total - 1);
			this->buildLists(); // All 0's

			while(currCover - appproxCoverPercent < 0) {
				// Place a square
				currIdx = distribution(generator);
				this->data[currIdx] = Color::Black;

				if (checkCellConsistency(currIdx)) {
					totalPlaced += 1.0;
					currCover = totalPlaced / floatTotal;
				} else {
					this->data[currIdx] = Color::White;
				}

			}
			this->buildHintStrides();
			this->buildLists();
		}while(!solve());
	}


	// Getters
	int getDimensions() { return numberDimensions;}

	// Copy
	vector<Color> getData() {return data;}
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

	vector<Color> getLine(struct LineDef& def, bool solved = true) {
		int length = this->dimensionsSize[def.dim];
		vector<Color> line;
		
		// Get offset
		int offset = 0;
		for (int i = 0; i < this->numberDimensions; i++) {
			if (i == def.dim) continue;
			if (def.idx[i] >= this->dimensionsSize[i] || def.idx[i] < 0) throw std::invalid_argument("idx content not valid for getLine");
			offset += def.idx[i] * this->strides[i];
		}
		if (solved)
			for (int i = 0; i < length; i++) 
				line.push_back(this->data[offset + i * this->strides[def.dim]]);				
		else
			for (int i = 0; i < length; i++)
				line.push_back(this->unsolvedData[offset + i * this->strides[def.dim]]);

		return line;
	}

	

	void printOriginal() {
		cout << endl << "Original:" << endl; 
		if (this->numberDimensions == 2) {
			for (int i = 0; i < this->dimensionsSize[1]; i++) {
				for (int j = 0; j < this->dimensionsSize[0]; j++) {
					int value = static_cast<int>(this->data[i * this->strides[1] + j]);
					cout << value << " ";
				} cout << endl;
			}cout << endl;
		}
		else {
			for (Color e : this->data) cout << static_cast<int>(e) << " ";
			cout << endl;
		}
	}

	void printSolved() {
		cout << endl << "This is what I got:" << endl; 
		if (this->numberDimensions == 2) {
			for (int i = 0; i < this->dimensionsSize[1]; i++) {
				for (int j = 0; j < this->dimensionsSize[0]; j++) {
					int value = static_cast<int>(this->unsolvedData[i * this->strides[1] + j]);
					if (value >= 0) cout << " " << value << " ";
					else cout << value << " ";
				} cout << endl;
			}cout << endl;
		}
		else {
			for (Color e : this->unsolvedData) cout << static_cast<int>(e) << " ";
			cout << endl;
		}
	}

	void printHints() {
		cout << "Hints: " << endl;
		for (vector<vector<int>> list : this->hints) {
			for (vector<int> lineHints : list) {
				for (int hint : lineHints) cout << hint << " ";
				cout << endl;
			}
			cout << endl;
		}
	}

};


#endif