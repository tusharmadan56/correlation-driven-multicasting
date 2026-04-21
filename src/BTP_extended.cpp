#include <bits/stdc++.h>
using namespace std;

#define INF 1e18

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
uniform_real_distribution<double> dist(0.0, 1.0);

struct Edge {
    int to;
    double p;
};

vector<Edge> adj[10];

double dp_L2[10][8];
double dp_L1[10][8];

int bitOf_L3(int n) { return (n >= 7 && n <= 9) ? n - 7 : -1; }
int bitOf_L2(int n) { return (n >= 4 && n <= 6) ? n - 4 : -1; }

void printVal(double x) {
    if (x >= INF / 2) cout << "INF";
    else cout << fixed << setprecision(4) << x;
}

string maskStr(int mask, int base) {
    string s = "{";
    bool first = true;
    for (int b = 0; b < 3; b++) {
        if (mask & (1 << b)) {
            if (!first) s += ",";
            s += to_string(b + base);
            first = false;
        }
    }
    return s + "}";
}

void computeDP(int u, double dp[][8], int (*bitMapper)(int)) {
    vector<Edge> edges;
    for (auto &e : adj[u])
        if (bitMapper(e.to) >= 0)
            edges.push_back(e);

    int k = edges.size();

    for (int mask = 0; mask < 8; mask++) {
        if (mask == 0) { dp[u][mask] = 0; continue; }

        bool ok = true;
        for (int b = 0; b < 3; b++) {
            if (!(mask & (1 << b))) continue;
            bool found = false;
            for (auto &e : edges)
                if (bitMapper(e.to) == b) { found = true; break; }
            if (!found) { ok = false; break; }
        }
        if (!ok) { dp[u][mask] = INF; continue; }

        double A = 0, B = 0;
        for (int m = 0; m < (1 << k); m++) {
            double pc = 1.0;
            int nm = mask;
            for (int i = 0; i < k; i++) {
                int b = bitMapper(edges[i].to);
                if (m & (1 << i)) {
                    pc *= edges[i].p;
                    if (b >= 0 && (nm & (1 << b))) nm ^= (1 << b);
                } else {
                    pc *= (1.0 - edges[i].p);
                }
            }
            if (nm == mask) A += pc;
            else B += pc * (1 + dp[u][nm]);
        }

        dp[u][mask] = (A >= 1.0 - 1e-12) ? INF : (A + B) / (1.0 - A);
    }
}

struct Assignment {
    map<int, int> node_dest;
    double cost;
};

Assignment optimalAssignment(int l2_mask) {
    vector<int> nodes;
    for (int b = 0; b < 3; b++)
        if (l2_mask & (1 << b))
            nodes.push_back(b + 4);

    Assignment best;
    best.cost = INF;

    auto canReach = [](int nd, int dest) {
        for (auto &e : adj[nd])
            if (e.to == dest) return true;
        return false;
    };

    function<void(int, map<int, int> &)> solve =
        [&](int d, map<int, int> &a) {
            if (d == 3) {
                double c = 0;
                for (auto &[n, m] : a) {
                    c += dp_L2[n][m];
                    if (c >= INF / 2) return;
                }
                if (c < best.cost) {
                    best.cost = c;
                    best.node_dest = a;
                }
                return;
            }
            int dest = d + 7;
            for (int n : nodes) {
                if (!canReach(n, dest)) continue;
                a[n] |= (1 << d);
                solve(d + 1, a);
                a[n] ^= (1 << d);
                if (a[n] == 0) a.erase(n);
            }
        };

    map<int, int> a;
    solve(0, a);
    return best;
}

struct Strategy {
    int fwd;
    int l2_tgt;
    Assignment asgn;
    double cost;
    double c0, c1, c2;
};

Strategy bestStrategy(int f) {
    double p0f = 0;
    for (auto &e : adj[0])
        if (e.to == f) { p0f = e.p; break; }
    double c0 = (p0f > 0) ? 1.0 / p0f : INF;

    vector<int> nbrs;
    for (auto &e : adj[f])
        if (bitOf_L2(e.to) >= 0)
            nbrs.push_back(e.to);

    int k = nbrs.size();
    Strategy best;
    best.cost = INF;
    best.fwd = f;

    for (int sub = 1; sub < (1 << k); sub++) {
        int lm = 0;
        for (int i = 0; i < k; i++)
            if (sub & (1 << i))
                lm |= (1 << bitOf_L2(nbrs[i]));

        int reach = 0;
        for (int b = 0; b < 3; b++)
            if (lm & (1 << b))
                for (auto &e : adj[b + 4])
                    if (bitOf_L3(e.to) >= 0)
                        reach |= (1 << bitOf_L3(e.to));
        if (reach != 7) continue;

        double c1 = dp_L1[f][lm];
        Assignment asgn = optimalAssignment(lm);
        double tot = c0 + c1 + asgn.cost;

        if (tot < best.cost)
            best = {f, lm, asgn, tot, c0, c1, asgn.cost};
    }

    return best;
}

int simMulticast(int u, int mask, int (*bmap)(int)) {
    int steps = 0;
    while (mask) {
        steps++;
        for (auto &e : adj[u]) {
            int b = bmap(e.to);
            if (b >= 0 && (mask & (1 << b)) && dist(rng) < e.p)
                mask ^= (1 << b);
        }
        if (steps > 100000) break;
    }
    return steps;
}

int simConservative(const Strategy &s) {
    int tot = 0;

    double p = 0;
    for (auto &e : adj[0])
        if (e.to == s.fwd) { p = e.p; break; }
    while (true) {
        tot++;
        if (dist(rng) < p) break;
        if (tot > 100000) return tot;
    }

    tot += simMulticast(s.fwd, s.l2_tgt, bitOf_L2);

    for (auto &[n, m] : s.asgn.node_dest)
        tot += simMulticast(n, m, bitOf_L3);

    return tot;
}

int main() {
    cout << fixed << setprecision(4);

    cout << "============================================================\n";
    cout << " BTP Extended: Correlation-Driven Multicasting (3 Levels)\n";
    cout << " Graph: 0 -> {1,2,3} -> {4,5,6} -> {7,8,9}\n";
    cout << " Destinations: {7, 8, 9}\n";
    cout << "============================================================\n\n";

    adj[0] = {{1, 1.0 / 3}, {2, 1.0 / 6}, {3, 1.0 / 2}};

    adj[1] = {{4, 1.0 / 3}, {5, 1.0 / 2}};
    adj[2] = {{4, 1.0 / 4}, {5, 1.0 / 3}, {6, 1.0 / 3}};
    adj[3] = {{5, 1.0 / 3}, {6, 1.0 / 4}};

    adj[4] = {{7, 1.0 / 3}, {8, 1.0 / 2}};
    adj[5] = {{8, 1.0 / 3}, {9, 1.0 / 4}};
    adj[6] = {{9, 1.0 / 3}};

    cout << "--- GRAPH EDGES ---\n";
    for (int u = 0; u <= 6; u++)
        for (auto &e : adj[u])
            cout << "  " << u << " -> " << e.to
                 << "  (p = " << e.p << ")\n";
    cout << "\n";

    for (int u = 7; u <= 9; u++)
        for (int mask = 0; mask < 8; mask++) {
            if (mask == 0)
                dp_L2[u][mask] = 0;
            else {
                int b = bitOf_L3(u);
                dp_L2[u][mask] = (mask == (1 << b)) ? 0 : INF;
            }
        }

    cout << "--- LEVEL 2 DP: nodes {4,5,6} -> destinations {7,8,9} ---\n";
    for (int u = 4; u <= 6; u++) {
        computeDP(u, dp_L2, bitOf_L3);
        cout << "Node " << u << ":\n";
        for (int mask = 1; mask < 8; mask++) {
            cout << "  dp_L2[" << u << "][" << maskStr(mask, 7) << "] = ";
            printVal(dp_L2[u][mask]);
            cout << "\n";
        }
    }
    cout << "\n";

    //cout << "--- LEVEL 1 DP: nodes {1,2,3} -> level-2 {4,5,6} ---\n";
    for (int u = 1; u <= 3; u++) {
        computeDP(u, dp_L1, bitOf_L2);
        cout << "Node " << u << ":\n";
        for (int mask = 1; mask < 8; mask++) {
            cout << "  dp_L1[" << u << "][" << maskStr(mask, 4) << "] = ";
            printVal(dp_L1[u][mask]);
            cout << "\n";
        }
    }
    cout << "\n";

    //cout << "--- FORWARDING SET ANALYSIS ---\n";

    vector<int> valid_sets, minimal_sets;

    for (int mask = 1; mask < 8; mask++) {
        int l2_reach = 0;
        for (int b = 0; b < 3; b++) {
            if (!(mask & (1 << b))) continue;
            int node = b + 1;
            for (auto &e : adj[node])
                if (bitOf_L2(e.to) >= 0)
                    l2_reach |= (1 << bitOf_L2(e.to));
        }

        int l3_reach = 0;
        for (int b = 0; b < 3; b++)
            if (l2_reach & (1 << b))
                for (auto &e : adj[b + 4])
                    if (bitOf_L3(e.to) >= 0)
                        l3_reach |= (1 << bitOf_L3(e.to));

        // cout << "  Set " << maskStr(mask, 1) << ": L2 reach = "
        //      << maskStr(l2_reach, 4) << ", L3 reach = "
        //      << maskStr(l3_reach, 7);

        if (l3_reach == 7) {
            //cout << "  -> VALID\n";
            valid_sets.push_back(mask);
        } else {
           // cout << "  -> INVALID\n";
        }
    }

    for (int s : valid_sets) {
        bool is_super = false;
        for (int t : valid_sets)
            if (t != s && (t & s) == t) { is_super = true; break; }
        if (!is_super) minimal_sets.push_back(s);
    }

    cout << "\nMinimal forwarding sets:\n";
    for (int m : minimal_sets)
        cout << "  " << maskStr(m, 1) << "\n";
    cout << "\n";

    cout << "--- COST COMPUTATION FOR MINIMAL SETS ---\n\n";

    vector<Strategy> strategies;

    for (int m : minimal_sets) {
        if (__builtin_popcount(m) == 1) {
            int b = __builtin_ctz(m);
            int f = b + 1;

            Strategy s = bestStrategy(f);
            strategies.push_back(s);

            cout << "Forwarding Set {" << f << "}:\n";
            cout << "  Cost 0->" << f << " = 1/p = ";
            printVal(s.c0);
            cout << "\n";

            cout << "  Level-2 target: " << maskStr(s.l2_tgt, 4)
                 << ", cost (L1 multicast) = ";
            printVal(s.c1);
            cout << "\n";

            cout << "  Destination assignment (L2->L3):\n";
            for (auto &[node, dmask] : s.asgn.node_dest) {
                cout << "    Node " << node << " -> "
                     << maskStr(dmask, 7) << ", cost = ";
                printVal(dp_L2[node][dmask]);
                cout << "\n";
            }
            cout << "  Assignment cost = ";
            printVal(s.c2);
            cout << "\n";

            cout << "  >>> TOTAL COST = ";
            printVal(s.cost);
            cout << "\n\n";
        }
    }

    Strategy *best = nullptr;
    for (auto &s : strategies)
        if (!best || s.cost < best->cost) best = &s;

    if (best) {
        cout << "=== BEST STRATEGY: Forwarding Set {"
             << best->fwd << "} ===\n";
        cout << "  Total analytical cost = ";
        printVal(best->cost);
        cout << "\n\n";
    }

    cout << "--- SIMULATION VERIFICATION ---\n";
    int runs = 100000;

    for (auto &s : strategies) {
        double totCons = 0;
        for (int i = 0; i < runs; i++) {
            totCons += simConservative(s);
        }
        double avgCons = totCons / runs;

        cout << "Forwarding Set {" << s.fwd << "}:\n";
        cout << "  Analytical (DP)     = ";
        printVal(s.cost);
        cout << "\n";
        cout << "  Sim Cost            = " << avgCons
             << "  (" << runs << " runs)\n";
        cout << "\n";
    }

    return 0;
}
