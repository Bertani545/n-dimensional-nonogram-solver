#include<iostream>
#include <fstream>
#include<vector>


using namespace std;

struct LISTS {
	vector<vector<int>> vertical;
	vector<vector<int>> horizontal;
};


void solve(int n, int m, struct LISTS lists) {

}


struct LISTS* build_lists(vector<vector<int>>& grid){
	struct LISTS* lists = new struct LISTS;
	int n = grid.size();
	int m = grid[0].size();
	lists->horizontal = vector<vector<int>>(n);
	lists->vertical = vector<vector<int>>(m);

	
	int total_numbers;
	int numbers[max(n, m) / 2 + 1];

	// Build horizonta lists
	for (int i = 0; i < n; i++) {
		
		int last_number = grid[i][0];
		int segment_size = (last_number == 1) ? 1 : 0;
		total_numbers = (last_number == 1) ? 1 : 0;
		for (int j = 1; j < m; j++) {
			if (grid[i][j]) { // 1
				if (grid[i][j] == last_number) {
					segment_size++;
				} else {
					// New block
					segment_size = 1;
					total_numbers++;
				}
			} else { // 0
				if (grid[i][j] != last_number) { // 1 and then 0
					numbers[total_numbers - 1] = segment_size;
				}
			}
			last_number = grid[i][j];
		}
		if (last_number == 1) numbers[total_numbers - 1] = segment_size; // Last correction

		if (total_numbers) {
			for (int k = 0; k < total_numbers; k++) lists->horizontal[i].push_back(numbers[k]);
		} else {
			// No squares
			lists->horizontal[i].push_back(0);
		}
	}

	// Vertical lists
	for (int j = 0; j < m; j++) {
		
		int last_number = grid[0][j];
		int segment_size = (last_number == 1) ? 1 : 0;
		total_numbers = (last_number == 1) ? 1 : 0;
		for (int i = 1; i < n; i++) {
			if (grid[i][j]) { // 1
				if (grid[i][j] == last_number) {
					segment_size++;
				} else {
					// New block
					segment_size = 1;
					total_numbers++;
				}
			} else { // 0
				if (grid[i][j] != last_number) { // 1 and then 0
					numbers[total_numbers - 1] = segment_size;
				}
			}
			last_number = grid[i][j];
		}
		if (last_number == 1) numbers[total_numbers - 1] = segment_size; // Last correction

		if (total_numbers) {
			for (int k = 0; k < total_numbers; k++) lists->vertical[j].push_back(numbers[k]);
		} else {
			// No squares
			lists->vertical[j].push_back(0);
		}
	}

	// Check them
	for (vector<int> v : lists->vertical){
		for ( int e : v)
			cout << e << " ";
		cout << endl;
	}cout << endl;
	for (vector<int> v : lists->horizontal) {
		for ( int e : v)
			cout << e << " ";
		cout << endl;
	}
		


	return lists;
}


struct LISTS* read_puzzle(ifstream& file) {
	int n = 0, m = 0;
	if (!(file >> n >> m)) return NULL;
	if (!n || !m) return NULL;
	printf("%d %d\n", n,m);
	vector<vector<int>> grid(n, vector<int>(m));

	for (int i = 0; i < n; i ++) {
		for (int j = 0; j < m; j++) {
			if (! (file >> grid[i][j])) return NULL;
		}
	}
	file.close();

	for (int i = 0; i < n; i ++) {
		for (int j = 0; j < m; j++) {
			printf("%d ", grid[i][j]);
		}printf("\n");
	}printf("\n");
	
	return build_lists(grid);
}

void read_lists(ifstream& file) {
	int n = 0, m = 0;

	file.close();
}




int main(int argc, char const *argv[])
{
	ifstream nonogram("Examples/nonogram.txt");
	if (!nonogram.is_open()) { // check if opened successfully
        std::cerr << "Could not open the file!\n";
        return 1;
    }

	struct LISTS* lists = read_puzzle(nonogram);

	delete lists;
	return 0;
}