#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

const double INF = 1e18;


int n, m;
vector<vector<int>> adj;
vector<vector<double>> cost;
vector<vector<double>> prob;

// check edge
bool hasEdge(int u, int v) {
    for (int x : adj[u])
        if (x == v) return true;
    return false;
}

//  probability of edge
double probEdge(int u, int v) {
    return prob[u][v];
}


double solve(int u, vector<int> dests) {

    if(dests.size()==0) return 0;

    // BASE CASE
    if (dests.size() == 1) {
        int d = dests[0];

        // if direct edge exists, return its cost
        if (hasEdge(u, d)) {
            return cost[u][d];
        }
        // else: fall through and solve normally
    }

    // -------- CHECK: is any destination directly reachable? 
    bool destReachable = false;
    for (int d : dests) {
        if (hasEdge(u, d)) {
            destReachable = true;
            break;
        }
    }


    // cout<<u<<endl;

    // for(auto x:dests){
    //     cout<<x<<" ";
    // }

    // cout<<endl;

    // cout<<destReachable<<endl;

    // ================= CASE 1 =================
    // at least one destination reachable from u


    if (destReachable) {

        double A = 0.0;
        double B = 0.0;

        int k = adj[u].size();



        for (int mask = 0; mask < (1 << k); mask++) {
            double pcur = 1.0;
            vector<int> remaining = dests;

            for (int i = 0; i < k; i++) {
                int v = adj[u][i];
                double p = probEdge(u, v);

                if (mask & (1 << i)) {
                    pcur *= p;
                    auto it = find(remaining.begin(), remaining.end(), v);
                    if (it != remaining.end())
                        remaining.erase(it);
                } else {
                    pcur *= (1.0 - p);
                }
            }

            if (remaining.size() == dests.size()) {
                A += pcur;
            } else {
                B += pcur * (1 + solve(u, remaining));
            }
        }

        return (A + B) / (1.0 - A);
    }

    // ================= CASE 2 =================
    // no destination reachable from u
    double A = 0.0;
    double B = 0.0;
    int k = adj[u].size();

    for (int mask = 0; mask < (1 << k); mask++) {
        double pcur = 1.0;
        double best = INF;

        for (int i = 0; i < k; i++) {
            int v = adj[u][i];
            double p = probEdge(u, v);

            //cout<<v<<" "<<p<<endl;

            if (mask & (1 << i)) {
                pcur *= p;
                best = min(best, solve(v, dests));
            } else {
                pcur *= (1.0 - p);
            }
        }



        if (best == INF) {
            A += pcur;
        } else {
            B += pcur * (1 + best);
        }

        //cout<<A<<" "<<B<<endl;
    }

    return (A + B) / (1.0 - A);
}

// ================= MAIN =================
int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    #ifndef ONLINE_JUDGE
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
        freopen("Error.txt", "w", stderr);
    #endif

    cin >> n >> m;

    adj.resize(n);
    cost.assign(n, vector<double>(n, 0.0));
    prob.assign(n, vector<double>(n, 0.0));

    // input edges
    // u v cost probability
    for (int i = 0; i < m; i++) {
        int u, v;
        double c;
        cin >> u >> v >> c;
        adj[u].push_back(v);
        cost[u][v] = c;
        prob[u][v] = 1.0/c;
    }

    int source;
    cin >> source;

    int k;
    cin >> k;
    vector<int> dests(k);
    for (int i = 0; i < k; i++)
        cin >> dests[i];

    double ans = solve(source, dests);
    cout << fixed << setprecision(6) << ans << "\n";
}