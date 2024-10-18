#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using ld = long double;

struct Place {
    ld y, x;
    Place() {}
    Place(ld y_, ld x_) : y(y_), x(x_) {}
    Place operator+(const std::vector<ld>& v) const { return Place(y + v[0], x + v[1]); }
    std::string output() { return "(" + std::to_string(y) + " " + std::to_string(x) + ")"; }

    ld calc_distance(const Place& other) const {
        return std::sqrt((y - other.y) * (y - other.y) + (x - other.x) * (x - other.x));
    }
};

class KDTree {
   public:
    KDTree(const std::vector<Place>& places) { root = buildTree(places, 0); }

    Place nearestNeighbor(const Place& target) { return nearest(root, target, 0).place; }

   private:
    struct KDNode {
        Place place;
        KDNode* left;
        KDNode* right;

        KDNode(const Place& p) : place(p), left(nullptr), right(nullptr) {}
    };
    KDNode* root;

    KDNode* buildTree(const std::vector<Place>& places, int depth) {
        if (places.empty()) return nullptr;

        auto placesCopy = places;
        int axis = depth % 2;
        auto comparator = [axis](const Place& a, const Place& b) {
            return (axis == 0 ? a.y < b.y : a.x < b.x);
        };
        std::sort(placesCopy.begin(), placesCopy.end(), comparator);
        size_t median = placesCopy.size() / 2;

        KDNode* node = new KDNode(placesCopy[median]);
        std::vector<Place> leftPlaces(placesCopy.begin(), placesCopy.begin() + median);
        std::vector<Place> rightPlaces(placesCopy.begin() + median + 1, placesCopy.end());

        node->left = buildTree(leftPlaces, depth + 1);
        node->right = buildTree(rightPlaces, depth + 1);
        return node;
    }

    KDNode nearest(KDNode* node, const Place& target, int depth) {
        if (!node) {
            return {Place()};
        }

        int axis = depth % 2;
        KDNode* nextBranch = nullptr;
        KDNode* otherBranch = nullptr;

        if ((axis == 0 && target.y < node->place.y) || (axis == 1 && target.x < node->place.x)) {
            nextBranch = node->left;
            otherBranch = node->right;
        } else {
            nextBranch = node->right;
            otherBranch = node->left;
        }

        KDNode best = nearest(nextBranch, target, depth + 1);
        if (best.place.calc_distance(target) == 0 ||
            node->place.calc_distance(target) < best.place.calc_distance(target)) {
            best = *node;
        }

        if (otherBranch &&
            std::abs((axis == 0 ? target.y - node->place.y : target.x - node->place.x)) <
                best.place.calc_distance(target)) {
            KDNode possibleBetter = nearest(otherBranch, target, depth + 1);
            if (possibleBetter.place.calc_distance(target) < best.place.calc_distance(target)) {
                best = possibleBetter;
            }
        }
        return best;
    }
};

int main() {
    std::vector<Place> places = {
        {2.1, 3.2}, {5.4, 2.1}, {9.2, 8.3}, {4.1, 7.3}, {8.4, 1.2},
    };

    KDTree tree(places);

    Place target = {2.0, 3.0};
    Place nearest = tree.nearestNeighbor(target);

    std::cout << "Nearest place to " << target.output() << " is " << nearest.output() << "\n";

    return 0;
}
