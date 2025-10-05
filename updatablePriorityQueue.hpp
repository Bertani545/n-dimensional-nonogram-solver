#ifndef UDQUEUE
#define UDQUEUE

#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>

#include"dataDef.hpp"

// ---------------- Priority Queue Wrapper ----------------
class UpdatablePriorityQueue {
    using SetType = std::set<Item, Compare>;
    SetType pq;
    std::unordered_map<Key, SetType::iterator, KeyHash> lookup;

public:
    bool empty() const { return pq.empty(); }
    size_t size() const { return pq.size(); }

    // Insert or update an element
    void insertOrUpdate(float priority, bool solvability, struct LineDef line) {
        Key key{line};

        if (auto it = lookup.find(key); it != lookup.end()) {
            pq.erase(it->second);
            lookup.erase(it);
        }

        auto [iter, _] = pq.insert({line, priority, solvability});
        lookup[key] = iter;
    }

    // Get top element (does not remove it)
    const Item& top() const {
        return *pq.begin();
    }

    // Pop top element
    void pop() {
        if (pq.empty()) return;
        auto top = *pq.begin();
        Key key{get<0>(top)};
        lookup.erase(key);
        pq.erase(pq.begin());
    }

    // Check if element exists
    bool contains(struct LineDef line) const {
        Key key{line};
        return lookup.find(key) != lookup.end();
    }
};


#endif