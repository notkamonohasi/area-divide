#include <bits/stdc++.h>
#define ll long long
#define ld long double
#define rep(i, n) for (ll i = 0; i < (ll)(n); i++)
#define frep(i, k, n) for (ll i = (ll)k; i < (ll)(n); i++)
#define brep(i, n) for (ll i = (ll)(n - 1); i >= 0; i--)
#define irep(it, mp) for (auto it = mp.begin(); it != mp.end(); it++)
#define pprint(x) cout << fixed << setprecision(12) << x << endl
#define BS(st, key) (st.find(key) != st.end())
#define ALL(v) v.begin(), v.end()
using namespace std;

void print() { cout << endl; }
template <class Head, class... Tail>
void print(Head&& head, Tail&&... tail) {
    if (sizeof...(tail) == 0)
        cout << head << endl;
    else {
        cout << head << " ";
        print(std::forward<Tail>(tail)...);
    }
}
template <typename T>
void vprint(vector<T>& v) {
    rep(i, v.size()) cout << v[i] << (i != v.size() - 1 ? " " : "");
    cout << endl;
}

using Graph = vector<vector<ll>>;

struct Edge {
    ll from, to;
    ld cost;

    Edge() {}
    Edge(ll from_, ll to_, ld cost_) : from(from_), to(to_), cost(cost_) {}
    bool operator>(const Edge& another) const { return cost > another.cost; }
    bool operator<(const Edge& another) const { return cost < another.cost; }
};
using EdgeGraph = vector<vector<Edge>>;

struct Place {
    ld y, x;
    Place() {}
    Place(ld y_, ld x_) : y(y_), x(x_) {}

    Place operator+(const vector<ld>& v) const { return Place(y + v[0], x + v[1]); }
    string output() { return "(" + to_string(y) + " " + to_string(x) + ")"; }

    ld calc_distance(const Place& other) const {
        return std::sqrt((y - other.y) * (y - other.y) + (x - other.x) * (x - other.x));
    }
};

namespace Algorithm {
struct Dijkstra {
    EdgeGraph g;
    ll start;
    vector<ld> dist;
    vector<ll> prev;

    const ld INF = INT_MAX;

    Dijkstra(const EdgeGraph& g_, ll start_) : g(g_), start(start_) { init(); }

    void init() {
        using P = pair<ld, ll>;
        ll N = g.size();
        dist.resize(N, INF);
        prev.resize(N, -1);            // 初期化
        vector<bool> visit(N, false);  // 確定
        priority_queue<P, vector<P>, greater<P>> pq;
        dist[start] = 0.0;
        pq.push(P(dist[start], start));

        while (!pq.empty()) {
            P p = pq.top();
            pq.pop();
            ll v = p.second;
            if (visit[v]) continue;  // これがないと確定済み地点を何回も見ることになる
            visit[v] = true;
            for (auto& e : g[v]) {
                ld next_cost = dist[v] + e.cost;
                if (!visit[e.to] && dist[e.to] > next_cost) {
                    dist[e.to] = next_cost;
                    prev[e.to] = v;  // 頂点 v を通って e.to にたどり着いた
                    pq.push(P(next_cost, e.to));
                }
            }
        }
    }

    // goalからstartへの経路を復元する
    vector<ll> get_path(ll goal) {
        vector<ll> ret;
        ll pos = goal;
        assert(dist[goal] != INF);
        while (pos >= 0) {
            ret.push_back(pos);
            pos = prev[pos];
        }
        reverse(ALL(ret));
        return ret;
    }

    // startからtへの距離を取得する
    ld get_dist(ll t) { return dist[t]; }
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
}  // namespace Algorithm

//
// ==========================================================================================================================================

// param
constexpr ld EPS = 1e-18;

namespace Utils {
ld calc_chebyshev_dist(Place& p1, Place& p2) { return max(abs(p1.y - p2.y), abs(p1.x - p2.x)); }
ld calc_euclid_dist(Place& p1, Place& p2) {
    return sqrt(pow(p1.y - p2.y, 2) + pow(p1.x - p2.x, 2));
}
};  // namespace Utils

struct NumberGenerator {
    ll pos = 0;

    NumberGenerator() {}

    ll gen() {
        pos++;
        return pos - 1;
    }
};

struct Obstacle {
    ll number;
    Place place;

    Obstacle() {}
    Obstacle(ll number_, Place place_) : number(number_), place(place_) {}
};

struct Area {
    ll number;
    Place place;                                // 正方形の中心
    ld length;                                  // 正方形の一辺の「半分」の長さ
    vector<Obstacle> obstacle_vec;              // 内部にある障害物
    unordered_set<ll> neighbor_area_number_st;  // 隣接するareaの番号

    Area() {}
    Area(ll number_, Place place_, ld length_) : number(number_), place(place_), length(length_) {
        assert(length > EPS * 10);  // 誤差で壊れるため
    }

    void add_neighbor(ll area_number) {
        ll before = neighbor_area_number_st.size();
        neighbor_area_number_st.insert(area_number);
        assert(neighbor_area_number_st.size() == before + 1);
    }

    void del_neighbor(ll area_number) {
        auto it = neighbor_area_number_st.find(area_number);
        assert(it != neighbor_area_number_st.end());
        neighbor_area_number_st.erase(it);
    }

    void add_obstacle(Obstacle& obstacle) {
        assert(is_inside(obstacle.place));
        obstacle_vec.push_back(obstacle);
    }

    // area同士が隣接しているか
    bool is_neighbor(Area& another) {
        ld y = place.y;
        ld x = place.x;
        ld l = length;
        ld dy = another.place.y;
        ld dx = another.place.x;
        ld dl = another.length;

        /**
         * NOTE: EPS小さく見積もることで、対角線の関係を除外する
         * NOTE: EPS大きく見積もることで、隣接させる
         */

        // x方向に接している
        if (min(y + l, dy + dl) > max(y - l, dy - dl) + EPS) {
            if (abs(x - dx) <= l + dl + EPS) return true;
        }

        // y方向に接している
        if (min(x + l, dx + dl) > max(x - l, dx - dl) + EPS) {
            if (abs(y - dy) <= l + dl + EPS) return true;
        }

        return false;
    }

    // pがこの中にあるか（辺上がどちらになるかは運）
    bool is_inside(Place& p) { return Utils::calc_chebyshev_dist(place, p) < length; }

    ll get_n_obstacles() const { return obstacle_vec.size(); }

    bool is_obstacle_empty() const { return obstacle_vec.empty(); }
};

struct State {
    unordered_map<ll, Area> area_mp;
    NumberGenerator area_number_generator;

    State() {}
    State(vector<Obstacle> obstacle_vec, ld xy) {
        // NOTE: 正方形に限定
        Area area(area_number_generator.gen(), Place(xy / 2.0, xy / 2.0), xy / 2.0);
        for (Obstacle& obstacle : obstacle_vec) area.add_obstacle(obstacle);

        area_mp[area.number] = area;
    }

    // 対象のareaを4分割する
    // 追加したAreaを返す
    vector<Area> divide_area(ll area_number) {
        assert(BS(area_mp, area_number));

        Area target_area = area_mp[area_number];

        // 削除
        area_mp.erase(area_number);
        irep(it, target_area.neighbor_area_number_st) {
            area_mp[*it].neighbor_area_number_st.erase(target_area.number);
        }

        vector<Area> next_area_vec;
        const vector<vector<ld>> NEIGHBOR_VEC = {{1, 1}, {1, -1}, {-1, -1}, {-1, 1}};
        for (auto v : NEIGHBOR_VEC) {
            ld l = target_area.length / 2.0;
            Place place = target_area.place + vector<ld>{l * v[0], l * v[1]};

            Area area(area_number_generator.gen(),
                      target_area.place + vector<ld>{l * v[0], l * v[1]}, l);

            // 障害物
            for (Obstacle& obstacle : target_area.obstacle_vec) {
                if (area.is_inside(obstacle.place)) area.add_obstacle(obstacle);
            }

            // area
            irep(it, target_area.neighbor_area_number_st) {
                assert(BS(area_mp, *it));
                if (area.is_neighbor(area_mp[*it])) {
                    area.add_neighbor(*it);
                    area_mp[*it].add_neighbor(area.number);
                }
            }

            next_area_vec.push_back(area);
        }

        // 今回追加するarea同士の隣接関係
        rep(i, 4) {
            frep(j, i + 1, 4) {
                if (next_area_vec[i].is_neighbor(next_area_vec[j])) {
                    assert(abs(i - j) != 2);  // debug
                    next_area_vec[i].add_neighbor(next_area_vec[j].number);
                    next_area_vec[j].add_neighbor(next_area_vec[i].number);
                }
            }
        }

        // 処理が終わったので記録
        for (Area& area : next_area_vec) {
            assert(!BS(area_mp, area.number));
            area_mp[area.number] = area;
        }

        return next_area_vec;
    }

    // placeを含むarea_numberを探す
    ll find_area_number(Place& place) {
        ll ret = -1;
        irep(it, area_mp) {
            if (it->second.is_inside(place)) {
                assert(ret < 0);
                ret = it->second.number;
            }
        }
        assert(ret >= 0);
        return ret;
    }

    ll get_n_areas() { return area_mp.size(); }

    void print_state() {
        print(get_n_areas());
        irep(it, area_mp) {
            Area area = it->second;
            Place place = it->second.place;
            ld l = it->second.length;
            print(place.y - l, place.x - l, place.y + l, place.x + l);
        }
    }
};

struct AreaDivisionSolver {
    State state;

    AreaDivisionSolver(vector<Obstacle> obstacle_vec, ld xy) { state = State(obstacle_vec, xy); }

    void solve(ll max_n_areas) {
        // 障害物を含む点のうち、最も面積の大きなareaを分割
        auto compare = [](const Area& a, const Area& b) {
            return (!a.is_obstacle_empty()) * pow(a.length, 2) <
                   (!b.is_obstacle_empty()) * pow(b.length, 2);
        };
        priority_queue<Area, vector<Area>, decltype(compare)> pq(compare);
        pq.push(state.area_mp.begin()->second);

        while (state.get_n_areas() < max_n_areas) {
            assert(!pq.empty());

            Area target_area = pq.top();
            pq.pop();

            vector<Area> next_area_vec = state.divide_area(target_area.number);
            for (Area& area : next_area_vec) pq.push(area);
        }
    }
};

struct RouteSolver {
    vector<Obstacle> obstacle_vec;
    State state;
    Place start_place, goal_place;

    RouteSolver(vector<Obstacle>& obstacle_vec_, State& state_, Place start_place_,
                Place goal_place_) {
        obstacle_vec = obstacle_vec_;
        state = state_;
        start_place = start_place_;
        goal_place = goal_place_;
    }

    void solve(ll max_n_areas) {
        while (state.get_n_areas() < max_n_areas) solve_turn();

        vector<ll> ans_path_vec = solve_shortest_path(state.find_area_number(start_place),
                                                      state.find_area_number(goal_place));

        // 出力
        state.print_state();
        print_path(ans_path_vec);
    }

    void solve_turn() {
        vector<ll> path_vec = solve_shortest_path(state.find_area_number(start_place),
                                                  state.find_area_number(goal_place));
        for (ll area_number : path_vec) state.divide_area(area_number);
    }

    vector<ll> solve_shortest_path(ll start_area_number, ll goal_area_number) {
        unordered_map<ll, ll> area_number_mp;  // map -> vector
        vector<ll> inv_area_number_vec;        // vector -> map
        ll pos = 0;
        irep(it, state.area_mp) {
            ll area_number = it->second.number;
            area_number_mp[area_number] = pos;
            inv_area_number_vec.push_back(area_number);
            pos++;
        }
        assert(area_number_mp.find(start_area_number) != area_number_mp.end());
        assert(area_number_mp.find(goal_area_number) != area_number_mp.end());
        ll start_area_idx = area_number_mp[start_area_number];
        ll goal_area_idx = area_number_mp[goal_area_number];

        vector<Place> obstacle_place_vec;
        for (Obstacle& obstacle : obstacle_vec) obstacle_place_vec.push_back(obstacle.place);
        Algorithm::KDTree kdtree(obstacle_place_vec);

        auto calc_obstacle_cost = [&](Place& place) -> ld {
            Place nearest_obstacle_place = kdtree.nearestNeighbor(place);
            ld dist = Utils::calc_euclid_dist(place, nearest_obstacle_place);
            if (dist < 0.20)
                return 1e3;
            else
                return 0.0;
        };

        vector<vector<Edge>> edge_graph(area_number_mp.size());
        irep(it, state.area_mp) {
            if (!it->second.obstacle_vec.empty()) continue;
            ll from_area_number = it->first;
            ll from_area_idx = area_number_mp[from_area_number];
            Place from_place = it->second.place;
            for (ll to_area_number : it->second.neighbor_area_number_st) {
                assert(state.area_mp.find(to_area_number) != state.area_mp.end());
                if (!state.area_mp[to_area_number].obstacle_vec.empty()) continue;
                ll to_area_idx = area_number_mp[to_area_number];
                Place to_place = state.area_mp[to_area_number].place;

                // cost
                ld dist_cost = Utils::calc_euclid_dist(from_place, to_place);
                ld obstacle_cost = calc_obstacle_cost(to_place);
                ld cost = dist_cost + obstacle_cost;

                edge_graph[from_area_idx].push_back(Edge(from_area_idx, to_area_idx, cost));
            }
        }

        Algorithm::Dijkstra dijkstra(edge_graph, start_area_idx);
        vector<ll> idx_path_vec = dijkstra.get_path(goal_area_idx);

        // 元の番号に戻す
        vector<ll> path_vec;
        for (ll x : idx_path_vec) path_vec.push_back(inv_area_number_vec[x]);

        return path_vec;
    }

    void print_path(vector<ll> path_vec) {
        print(path_vec.size() + 2);
        print(start_place.y, start_place.x);
        for (ll area_number : path_vec) {
            Place place = state.area_mp[area_number].place;
            print(place.y, place.x);
        }
        print(goal_place.y, goal_place.x);
    }
};

//
// =========================================================================================================
vector<Obstacle> OBSTACLE_VEC;
ld XY;
ll FIRST_MAX_N_AREAS, SECOND_MAX_N_AREAS;
Place START_PLACE, GOAL_PLACE;
void input() {
    ll n_obstacles;
    cin >> XY >> n_obstacles >> FIRST_MAX_N_AREAS >> SECOND_MAX_N_AREAS;
    cin >> START_PLACE.y >> START_PLACE.x >> GOAL_PLACE.y >> GOAL_PLACE.x;

    NumberGenerator obstacle_number_generator;

    rep(_, n_obstacles) {
        ld y, x;
        cin >> y >> x;
        OBSTACLE_VEC.push_back(Obstacle(obstacle_number_generator.gen(), Place(y, x)));
    }
}

int main() {
    input();
    AreaDivisionSolver area_division_solver(OBSTACLE_VEC, XY);
    area_division_solver.solve(FIRST_MAX_N_AREAS);
    // area_division_solver.state.print_state();

    RouteSolver route_solver(OBSTACLE_VEC, area_division_solver.state, START_PLACE, GOAL_PLACE);
    route_solver.solve(SECOND_MAX_N_AREAS);
}
