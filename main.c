#include<stdio.h>
#include <stdlib.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

struct LISTS {
	int** vertical;
	int** horizontal;
	int n;
	int m;
};


void free_lists(struct LISTS *l)
{
	for (int i = 0; i < l->n; i++) {
	    free(l->horizontal[i]);
	}
	free(l->horizontal);
	for (int i = 0; i < l->m; i++) {
	    free(l->vertical[i]);
	}
	free(l->vertical);
	
	l->vertical = NULL;
	l->horizontal = NULL;
}

void solve(int n, int m, struct LISTS lists) {

}


struct LISTS* build_lists(int n, int m, int grid[n][m]){
	int** vert = (int**) malloc(m * sizeof(int*));
	int** horz = (int**) malloc(n * sizeof(int*));

	// Build horizonta lists
	int total_numbers;
	int numbers[MAX(n,m) / 2 + 1];
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
			horz[i] = (int*) malloc(total_numbers * sizeof(int));
			for (int k = 0; k < total_numbers; k++) horz[i][k] = numbers[k];
		} else {
			// No squares
			horz[i] = (int*) malloc(sizeof(int));
			horz[i][0] = 0;
		}
		

		for (int k = 0; k < total_numbers; k++) printf("%d ", numbers[k]);
		printf("\n");
	}


	struct LISTS* lists = (struct LISTS*) malloc(sizeof(struct LISTS));
	lists->vertical = vert;
	lists->horizontal = horz;
	lists->n = n;
	lists->m = m;
	return lists;
}


struct LISTS* read_puzzle(FILE* file) {
	int n = 0, m = 0;
	if (fscanf(file, "%d %d", &n, &m) != 2) return NULL;
	if (!n || !m) return NULL;
	printf("%d %d\n", n,m);
	int grid[n][m];

	for (int i = 0; i < n; i ++) {
		for (int j = 0; j < m; j++) {
			if (!fscanf(file, "%d", &grid[i][j])) return NULL;
		}
	}
	fclose(file);

	for (int i = 0; i < n; i ++) {
		for (int j = 0; j < m; j++) {
			printf("%d ", grid[i][j]);
		}printf("\n");
	}printf("\n");
	
	return build_lists(n, m, grid);
}

void read_lists(FILE* file) {
	int n = 0, m = 0;

	fclose(file);
}




int main(int argc, char const *argv[])
{
	FILE* nonogram = fopen("Examples/nonogram.txt", "r");
	if (nonogram == NULL) {
		perror("Error opening puzzle");
		return 1;
	}

	struct LISTS* lists = read_puzzle(nonogram);
	//free_lists(lists);

	return 0;
}