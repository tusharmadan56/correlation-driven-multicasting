#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

const double INF = 1e18;
const double EPS = 1e-9;

/* graph[u] = {v, p} */
vector<vector<pair<int,double>>> readGraph(int n, int m) {
    vector<vector<pair<int,double>>> graph(n);

    for (int i = 0; i < m; i++) {
        int u, v;
        double p;
        cin >> u >> v >> p;

        // store inverted probability
        graph[u].push_back({v, 1.0 / p});
    }
    return graph;
}

double calcCost(
    int u,
    int dest,
    const vector<double>& cost,
    const vector<vector<pair<int,double>>>& graph
) {
    if (u == dest) return 0.0;

    vector<pair<double,double>> next;

    for (auto &e : graph[u]) {
        int v = e.first;
        double p = e.second;

        if (cost[v] < cost[u]) {
            next.push_back({cost[v], p});
        }
    }

    if (next.empty()) return cost[u];

    sort(next.begin(), next.end());

    double fail = 1.0;
    double sum = 0.0;

    for (auto &x : next) {
        sum += fail * x.second * x.first;
        fail *= (1.0 - x.second);
    }

    double ok = 1.0 - fail;
    if (ok < EPS) return cost[u];

    return (1.0 + sum) / ok;
}

vector<double> solve(
    int n,
    int dest,
    const vector<vector<pair<int,double>>>& graph
) {
    vector<double> cost(n, INF);
    cost[dest] = 0.0;

    for (int it = 0; it < n - 1; it++) {
        bool changed = false;

        for (int u = 0; u < n; u++) {
            double val = calcCost(u, dest, cost, graph);

            if (val + EPS < cost[u]) {
                cost[u] = val;
                changed = true;
            }
        }

        if (!changed) break;
    }
    return cost;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    cin >> n >> m;

    auto graph = readGraph(n, m);

    int src, dest;
    cin >> src >> dest;

    auto cost = solve(n, dest, graph);

    cout << fixed << setprecision(6);

    if (cost[src] >= INF / 2) {
        cout << "Source cannot reach destination\n";
    } else {
        cout << cost[src] << "\n";
    }

    return 0;
}