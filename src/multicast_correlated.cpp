#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

const double INF = 1e18;


int n, m;
vector<vector<int>> adj;
vector<vector<double>> cost;
vector<double> psim; 


bool hasEdge(int u, int v) {
    for (int x : adj[u]) if (x == v) return true;
    return false;
}

// probability from inverse cost 
double probEdge(int u, int v) {
    return 1.0 / cost[u][v];
}


double solve(int u, vector<int> dests) {

    // ---- Base case 
    if(dests.size()==0) return 0;

    if (dests.size() == 1) {
        int d = dests[0];
        if (hasEdge(u, d)) {
            return cost[u][d];   // direct edge → done
        }
        // else fall through and solve normally
    }

    int k = adj[u].size();

  
    // SIMILARITY CASE 
    
    if (psim[u] > 0.0 && dests.size() == 1 && k >= 2) {
        //int d = dests[0];

        // For now, we use the first two outgoing edges 
        int v1 = adj[u][0];
        int v2 = adj[u][1];

        double p1 = probEdge(u, v1);
        double p2 = probEdge(u, v2);
        double s  = psim[u];

        double x1 = solve(v1, dests); 
        double x2 = solve(v2, dests); 

        // y = (1 - (p1+p2-s))*(1+y)
        //   + (p1-s)*(1+x1)
        //   + (p2-s)*(1+x2)
        //   + s*(1+min(x1,x2))

        //cout<<x1<<" "<<x2<<endl;

        double A = 1.0 - (p1 + p2 - s);
        double B = (p1 - s) * (1 + x1)
                 + (p2 - s) * (1 + x2)
                 + s * (1 + min(x1, x2));

       
        A = min(A, 1.0 - 1e-12); 
        return (A + B) / (1.0 - A);
    }

    
    // ORIGINAL ALGORITHM 
   

    // Check if any destination directly reachable
    bool destReachable = false;
    for (int d : dests) {
        if (hasEdge(u, d)) {
            destReachable = true;
            break;
        }
    }

    // ---------- CASE 1: destination reachable ----------
    if (destReachable) {
        double A = 0.0, B = 0.0;

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
                A += pcur; // no destination reached
            } else {
                B += pcur * (1 + solve(u, remaining));
            }
        }

        A = min(A, 1.0 - 1e-12);
        return (A + B) / (1.0 - A);
    }

    // ---------- CASE 2: no destination reachable ----------
    double A = 0.0, B = 0.0;

    for (int mask = 0; mask < (1 << k); mask++) {
        double pcur = 1.0;
        double best = INF;

        for (int i = 0; i < k; i++) {
            int v = adj[u][i];
            double p = probEdge(u, v);

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
    }

    A = min(A, 1.0 - 1e-12);
    return (A + B) / (1.0 - A);
}


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    #ifndef ONLINE_JUDGE
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
        freopen("Error.txt", "w", stderr);
    #endif

    cin >> n >> m;

    adj.resize(n);
    cost.assign(n, vector<double>(n, 0.0));
    psim.assign(n, 0.0);

    // Read edges: u v cost
    for (int i = 0; i < m; i++) {
        int u, v;
        double c;
        cin >> u >> v >> c;
        adj[u].push_back(v);
        cost[u][v] = c;
    }

    // Ask psim for nodes with outdegree > 1
    for (int u = 0; u < n; u++) {
        if ((int)adj[u].size() > 1) {
            //cout<<"Enter psim for node "<<u<<endl;
            cin >> psim[u]; 
        } else {
            psim[u] = 0.0;
        }
    }

    int source;
    cin >> source;

    int k;
    cin >> k;
    vector<int> dests(k);
    for (int i = 0; i < k; i++) cin >> dests[i];

    double ans = solve(source, dests);
    cout << fixed << setprecision(6) << ans << "\n";
    return 0;
}