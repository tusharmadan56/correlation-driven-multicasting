#include <bits/stdc++.h>
using namespace std;

#define INF 1e18

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
uniform_real_distribution<double> dist(0.0, 1.0);

struct Edge {
    int to;
    double p;
};

vector<Edge> adj[7];
double dp[7][8];

// 4->0, 5->1, 6->2
int bitOf(int node) {
    if (node == 4) return 0;
    if (node == 5) return 1;
    if (node == 6) return 2;
    return -1;
}

// ---------------- LEVEL 2 DP ----------------
void computeLevel2(int u) {

    for (int mask = 0; mask < 8; mask++) {

        if (mask == 0) {
            dp[u][mask] = 0;
            continue;
        }

        int k = adj[u].size();
        double A = 0.0, B = 0.0;

        for (int m = 0; m < (1 << k); m++) {

            double pcur = 1.0;
            int newMask = mask;

            for (int i = 0; i < k; i++) {
                int v = adj[u][i].to;
                double p = adj[u][i].p;

                if (m & (1 << i)) {
                    pcur *= p;
                    int b = bitOf(v);
                    if (mask & (1 << b))
                        newMask ^= (1 << b);
                } else {
                    pcur *= (1.0 - p);
                }
            }

            if (newMask == mask) {
                A += pcur;
            } else {
                B += pcur * (1 + dp[u][newMask]);
            }
        }

        if (A == 1.0) dp[u][mask] = INF;
        else dp[u][mask] = (A + B) / (1.0 - A);
    }
}

// ---------------- PRINT ----------------
void printVal(double x) {
    if (x >= INF/2) cout << "INF";
    else cout << fixed << setprecision(6) << x;
}

// ---------------- SOLVE {2} ----------------
double solve_2() {

    double p2 = adj[0][1].p; // 0->2

    double P_none = (1 - p2);
    double P_2 = p2;

    double E = 0;

    double c = dp[2][7]; // {4,5,6}

    if (c >= INF) E += P_2 * INF;
    else E += P_2 * (1 + c);

    if (P_none == 1.0) return INF;
    if (E >= INF/2) return INF;

    return (E + P_none) / (1.0 - P_none);
}

// ---------------- SOLVE {1,3} ----------------
double solve_13() {

    double p1 = adj[0][0].p;
    double p3 = adj[0][2].p;
    
    //cout<<p1<<" "<<p3<<endl;

    double P_none = (1-p1)*(1-p3);
    double P_1 = p1*(1-p3);
    double P_3 = (1-p1)*p3;
    double P_both = p1*p3;

    double E = 0;

    // 1 only → needs {4,5} from 1 and {6} from 3
    double c1 = dp[1][3] + dp[3][4] + 1/p3;
    
    //cout<<c1<<" "<<dp[1][3]<<" "<<dp[3][4]<<endl;

    if (c1 >= INF) E += P_1 * INF;
    else E += P_1 * (1 + c1);

    // 3 only → needs {5,6} from 3 and {4} from 1
    double c3 = dp[3][6] + dp[1][1] + 1/p1;
    
   // cout<<c3<<" "<<dp[3][6]<<" "<<dp[1][1]<<endl;

    if (c3 >= INF) E += P_3 * INF;
    else E += P_3 * (1 + c3);

    // both → choose better
    double best = min(c1, c3);

    if (best >= INF) E += P_both * INF;
    else E += P_both * (1 + best);

    if (P_none == 1.0) return INF;
    if (E >= INF/2) return INF;

    return (E + P_none) / (1.0 - P_none);
}

// ---------------- GENERIC SIMULATION ----------------
int simulate_node(int u, int mask) {

    int steps = 0;

    while (mask) {

        steps++;

        for (auto &e : adj[u]) {

            if (dist(rng) < e.p) {
                int b = bitOf(e.to);
                if (mask & (1 << b))
                    mask ^= (1 << b);
            }
        }

        if (steps > 10000) break;
    }

    return steps;
}

// ---------------- SIMULATE {2} ----------------
int simulate_once_2() {

    int steps = 0;
    double p2 = adj[0][1].p;

    while (true) {
        steps++;

        if (dist(rng) < p2) {
            return steps + simulate_node(2, 7);
        }
    }
}


int simulate_from_0_to_4() {

    int steps = 0;

    while (true) {

        double p=0;
        steps++;
        p=adj[0][0].p;

       // cout<<p<<endl;

        while(true){
            
            steps++;

            if(dist(rng) < p){
                break;
            }
            
        }

       // cout<<steps<<endl;

        p=adj[1][0].p;
        while(true){
            
            steps++;

            if(dist(rng) < p){
                break;
            }
            
        }

       // cout<<p<<endl;

        //cout<<steps<<endl;
        break;
    }

    return steps ;
}

int simulate_from_0_to_6() {

    int steps = 0;

    while (true) {

        double p=0;

        p=adj[0][2].p;
       // steps++;

        while(true){
            
             steps++;
            if(dist(rng) < p){
                break;
            }
           
        }

        p=adj[3][1].p;
        while(true){
            
            steps++;
            if(dist(rng) < p){
                break;
            }
            
        }
        break;
    }

    return steps ;
}



// ---------------- SIMULATE {1,3} ----------------
int simulate_once_13() {

    int steps = 0;

    double p1 = adj[0][0].p; // 0->1
    double p3 = adj[0][2].p; // 0->3

    double p_only1 = p1 * (1 - p3);
    double p_only3 = (1 - p1) * p3;
    double p_both  = p1 * p3;
    double p_none  = (1 - p1) * (1 - p3);

    double c1 = p_only1;
    double c2 = c1 + p_only3;
    double c3 = c2 + p_both;

    while (true) {

        steps++;

        double r = dist(rng);

        // ONLY 1
        if (r < c1) {

            int cost =
                simulate_node(1, 3) +   // {4,5}
                simulate_from_0_to_6();         // {6} via node 3

            return steps + cost;
        }

        // ONLY 3
        else if (r < c2) {

            int cost =
                simulate_node(3, 6) +   // {5,6}
                simulate_from_0_to_4();         // {4} via node 1

            return steps + cost;
        }

        // BOTH
        else if (r < c3) {

            int cost1 =
                simulate_node(1, 3) +
                simulate_from_0_to_6();

            int cost2 =
                simulate_node(3, 6) +
                simulate_from_0_to_4();

            return steps + min(cost1, cost2);
        }

        
    }
}

void run_simulation() {

    int runs = 10000;
    double total2 = 0, total13 = 0;

    for (int i = 0; i < runs; i++) {
        total2 += simulate_once_2();
        total13 += simulate_once_13();
    }

    cout << "\nSimulation {2} = " << total2 / runs << "\n";
    cout << "Simulation {1,3} = " << total13 / runs << "\n";
}


int main() {

    cout << fixed << setprecision(6);

    // GRAPH
    adj[0] = {
        {1, 1.0/3.0},
        {2, 1.0/6.0},
        {3, 1.0/2.0}
    };

    adj[1] = { {4,1.0/3.0}, {5,1.0/2.0} };

    
    adj[2] = {
        {4,1.0/4.0},
        {5,1.0/3.0},
        {6,1.0/3.0}   
    };

    adj[3] = { {5,1.0/3.0}, {6,1.0/4.0} };

    // BASE
    for (int u = 4; u <= 6; u++) {
        for (int mask = 0; mask < 8; mask++) {
            if (mask == 0) dp[u][mask] = 0;
            else {
                int b = bitOf(u);
                if (mask == (1 << b)) dp[u][mask] = 0;
                else dp[u][mask] = INF;
            }
        }
    }

    // DP
    computeLevel2(1);
    computeLevel2(2);
    computeLevel2(3);

    cout << "\n--- RESULTS ---\n";

    double ans2 = solve_2();
    double ans13 = solve_13();

    cout << "{2} -> ";
    printVal(ans2);
    cout << "\n";

    cout << "{1,3} -> ";
    printVal(ans13);
    cout << "\n";

    cout << "Best = ";
    if (ans2 < ans13) cout << "{2}\n";
    else cout << "{1,3}\n";

    run_simulation();

    return 0;
}