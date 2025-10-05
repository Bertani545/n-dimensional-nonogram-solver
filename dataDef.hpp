#ifndef DATADEF
#define DATADEF


#include<vector>
#include<tuple>

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

enum class HintState {
	unsolved = 0,
	solved = 1
};


// ---------------- Item definition ----------------
// Def, heuristic, solvability
using Item = tuple<struct LineDef, float, bool>;

// Comparator for set (higher priority first)
struct Compare {
    bool operator()(const Item& a, const Item& b) const {
    	if (get<1>(a) != get<1>(b))
    		return get<1>(a) < get<1>(b);

    	struct LineDef lineA = get<0>(a);
    	struct LineDef lineB = get<0>(b);

    	if (lineA.dim != lineB.dim)
    		return lineA.dim < lineB.dim;
        return lineA.idx < lineB.idx;   // tie-breaker
    }
};

// Key that ignores priority (only d2 + coords matter)
struct Key {
    struct LineDef line;

    Key() = default;
    Key(const struct LineDef& l) : line(l) {}

    bool operator==(const Key& other) const {
        if (line.dim != other.line.dim) return false;
        bool eqVec;
        for (int i = 0; i < line.idx.size(); i++) {
            if (i == line.dim) continue;
            if (line.idx[i] != other.line.idx[i])
                return false;
        }
        return true;
    }
};

// Hash for Key (for unordered_map)
struct KeyHash {
    std::size_t operator()(const Key& k) const {
        std::size_t h1 = std::hash<int>()(k.line.dim);
        std::size_t h2 = 0;
        for (int i = 0; i < k.line.idx.size(); i++) {
            if (i == k.line.dim) continue;
            int c = k.line.idx[i];
             h2 ^= std::hash<int>()(c) + 0x9e3779b9 + (h2 << 6) + (h2 >> 2);
        }
        return h1 ^ h2;
    }
};


#endif