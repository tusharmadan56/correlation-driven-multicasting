#include <bits/stdc++.h>
using namespace std;

#define INF 1e18

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
uniform_real_distribution<double> dist(0.0, 1.0);

// =====================================================================
//  BTP Extended: Correlation-Driven Multicasting (3-Level Graph)
//
//  Graph Topology:
//
//               0 (source)
//             / | \
//           1   2   3          (Level 1 - forwarding candidates)
//          /|  /|\  |\
//         4  5  6              (Level 2 - intermediate relays)
//        /|  |\  \
//       7  8  9                (Level 3 - destinations)
//
//  Edges:
//    Level 0->1: 0-1(1/3), 0-2(1/6), 0-3(1/2)
//    Level 1->2: 1-4(1/3), 1-5(1/2), 2-4(1/4), 2-5(1/3), 2-6(1/3),
//                3-5(1/3), 3-6(1/4)
//    Level 2->3: 4-7(1/3), 4-8(1/2), 5-8(1/3), 5-9(1/4), 6-9(1/3)
//
//  Heuristic: Enumerate all forwarding sets at Level 1,
//             prune supersets, keep minimal sets that reach all
//             destinations {7,8,9}.
// =====================================================================

struct Edge {
    int to;
    double p;
    double ett; // Expected Time per Transmission = 1/p (geometric mean trials)
};

vector<Edge> adj[10]; // nodes 0-9

// DP tables
// dp_L2[u][mask]: expected tx from level-2 node u to deliver destinations
//                 in mask (bitmask over {7,8,9}: 7->bit0, 8->bit1, 9->bit2)
// dp_L1[u][mask]: expected tx from level-1 node u to deliver level-2 nodes
//                 in mask (bitmask over {4,5,6}: 4->bit0, 5->bit1, 6->bit2)
double dp_L2[10][8];
double dp_L1[10][8];

// ==================== BIT MAPPINGS ====================
int bitOf_L3(int n) { return (n >= 7 && n <= 9) ? n - 7 : -1; }
int bitOf_L2(int n) { return (n >= 4 && n <= 6) ? n - 4 : -1; }

// ==================== UTILITIES ====================
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

// ==================== GENERIC DP ====================
// Compute expected transmissions from node u to deliver all targets in mask.
// bitMapper maps neighbor node IDs to bit positions (0,1,2).
// Works for both Level-2->Level-3 and Level-1->Level-2.
void computeDP(int u, double dp[][8], int (*bitMapper)(int)) {
    // Collect edges relevant to this level
    vector<Edge> edges;
    for (auto &e : adj[u])
        if (bitMapper(e.to) >= 0)
            edges.push_back(e);

    int k = edges.size();

    for (int mask = 0; mask < 8; mask++) {
        if (mask == 0) { dp[u][mask] = 0; continue; }

        // Check every bit in mask is reachable via some edge
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

// ==================== ASSIGNMENT OPTIMIZATION ====================
// Given a set of active level-2 nodes (l2_mask over {4,5,6}), find
// the minimum-cost assignment of destinations {7,8,9} to those nodes.
struct Assignment {
    map<int, int> node_dest; // level-2 node -> dest mask (bits over {7,8,9})
    double cost;
};

Assignment optimalAssignment(int l2_mask) {
    vector<int> nodes;
    for (int b = 0; b < 3; b++)
        if (l2_mask & (1 << b))
            nodes.push_back(b + 4);

    Assignment best;
    best.cost = INF;

    // Data-driven reachability
    auto canReach = [](int nd, int dest) {
        for (auto &e : adj[nd])
            if (e.to == dest) return true;
        return false;
    };

    // Enumerate assignments: each of 3 dests assigned to one available node
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

// ==================== STRATEGY ====================
struct Strategy {
    int fwd;        // level-1 forwarder node
    int l2_tgt;     // mask of level-2 nodes to target
    Assignment asgn; // destination assignment
    double cost;     // total analytical cost
    double c0, c1, c2; // cost breakdown: 0->f, f->L2, L2->L3
};

// Find the best strategy for a single level-1 forwarder f
Strategy bestStrategy(int f) {
    double p0f = 0;
    for (auto &e : adj[0])
        if (e.to == f) { p0f = e.p; break; }
    double c0 = (p0f > 0) ? 1.0 / p0f : INF;

    // Level-2 neighbors of f
    vector<int> nbrs;
    for (auto &e : adj[f])
        if (bitOf_L2(e.to) >= 0)
            nbrs.push_back(e.to);

    int k = nbrs.size();
    Strategy best;
    best.cost = INF;
    best.fwd = f;

    // Try every subset of f's level-2 neighbors
    for (int sub = 1; sub < (1 << k); sub++) {
        int lm = 0;
        for (int i = 0; i < k; i++)
            if (sub & (1 << i))
                lm |= (1 << bitOf_L2(nbrs[i]));

        // Check if these level-2 nodes can collectively reach all {7,8,9}
        int reach = 0;
        for (int b = 0; b < 3; b++)
            if (lm & (1 << b))
                for (auto &e : adj[b + 4])
                    if (bitOf_L3(e.to) >= 0)
                        reach |= (1 << bitOf_L3(e.to));
        if (reach != 7) continue; // need all 3 destinations

        double c1 = dp_L1[f][lm];
        Assignment asgn = optimalAssignment(lm);
        double tot = c0 + c1 + asgn.cost;

        if (tot < best.cost)
            best = {f, lm, asgn, tot, c0, c1, asgn.cost};
    }

    return best;
}

// ==================== SIMULATION ====================

// Simulate multicast from node u until all targets in mask are delivered
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

// Conservative simulation: each level-2 node independently delivers
// its assigned destinations. Matches analytical model exactly.
int simConservative(const Strategy &s) {
    int tot = 0;

    // Step 1: 0 -> forwarder (geometric)
    double p = 0;
    for (auto &e : adj[0])
        if (e.to == s.fwd) { p = e.p; break; }
    while (true) {
        tot++;
        if (dist(rng) < p) break;
        if (tot > 100000) return tot;
    }

    // Step 2: forwarder -> level-2 targets
    tot += simMulticast(s.fwd, s.l2_tgt, bitOf_L2);

    // Step 3: level-2 nodes -> destinations (independent)
    for (auto &[n, m] : s.asgn.node_dest)
        tot += simMulticast(n, m, bitOf_L3);

    return tot;
}

// Opportunistic simulation: level-2 nodes transmit sequentially.
// Each node delivers its EXCLUSIVE remaining destinations, but we
// track ALL destinations that happen to receive (broadcast bonus).
// This exploits wireless broadcast advantage across levels.
int simOpportunistic(const Strategy &s) {
    int tot = 0;

    // Step 1: 0 -> forwarder
    double p = 0;
    for (auto &e : adj[0])
        if (e.to == s.fwd) { p = e.p; break; }
    while (true) {
        tot++;
        if (dist(rng) < p) break;
        if (tot > 100000) return tot;
    }

    // Step 2: forwarder -> level-2 targets
    tot += simMulticast(s.fwd, s.l2_tgt, bitOf_L2);

    // Step 3: level-2 nodes deliver sequentially with cross-tracking
    int rem = 7; // remaining destinations {7,8,9} = 0b111
    vector<int> l2_nodes;
    for (int b = 0; b < 3; b++)
        if (s.l2_tgt & (1 << b))
            l2_nodes.push_back(b + 4);

    // Phase A: each node delivers exclusive remaining dests, tracking bonuses
    for (int node : l2_nodes) {
        if (rem == 0) break;

        // Find exclusive remaining destinations for this node
        int exclusive = 0;
        for (auto &e : adj[node]) {
            int b3 = bitOf_L3(e.to);
            if (b3 < 0 || !(rem & (1 << b3))) continue;
            bool only_me = true;
            for (int o : l2_nodes) {
                if (o == node) continue;
                for (auto &e2 : adj[o])
                    if (e2.to == e.to) { only_me = false; break; }
                if (!only_me) break;
            }
            if (only_me) exclusive |= (1 << b3);
        }

        if (exclusive == 0) continue;

        // Transmit until exclusive dests reached, track ALL received
        int reached = 0;
        while ((exclusive & ~reached) != 0) {
            tot++;
            for (auto &e : adj[node]) {
                int b3 = bitOf_L3(e.to);
                if (b3 >= 0 && (rem & (1 << b3)) && !(reached & (1 << b3)))
                    if (dist(rng) < e.p)
                        reached |= (1 << b3);
            }
            if (tot > 200000) return tot;
        }
        rem &= ~reached;
    }

    // Phase B: handle any remaining shared destinations
    for (int node : l2_nodes) {
        if (rem == 0) break;
        int target = 0;
        for (auto &e : adj[node]) {
            int b3 = bitOf_L3(e.to);
            if (b3 >= 0 && (rem & (1 << b3)))
                target |= (1 << b3);
        }
        if (target == 0) continue;
        tot += simMulticast(node, target, bitOf_L3);
        rem &= ~target;
    }

    return tot;
}

// =====================================================================
//  PACKET TRANSMISSION STRATEGY
//
//  For a given Strategy (forwarder f, L2 target mask, dest assignment),
//  we derive the complete step-by-step packet transmission plan:
//
//  Phase 0: Source (0) -> Forwarder (f)
//    - Keep retransmitting until f receives the packet.
//    - Expected slots: ETT(0->f) = 1/p(0,f)
//
//  Phase 1: Forwarder (f) -> Relay nodes (L2 targets)
//    - f multicasts to all assigned L2 relays simultaneously each slot.
//    - Each relay independently receives with its edge probability.
//    - Keep multicasting until ALL assigned relays have the packet.
//    - ETT-priority: relays with higher ETT (lower p) are the bottleneck.
//    - Expected slots: dp_L1[f][l2_mask]
//
//  Phase 2: Relay nodes (L2) -> Destinations (L3) [PARALLEL]
//    - Each relay independently multicasts to its assigned destinations.
//    - Relays operate in PARALLEL (time-division not needed; different links).
//    - Within each relay: multicast until all its dests have the packet.
//    - ETT-priority: destinations with higher ETT are served first if shared.
//    - Expected slots (bottleneck): max over all relays of dp_L2[relay][dest_mask]
//    - Total Phase 2 cost (conservative model): sum of dp_L2 (sequential relays)
//
//  Multirate Scheduling Rule:
//    At each phase, sort edges by ETT descending (worst link first).
//    The packet is broadcast on ALL edges simultaneously each slot —
//    this is the wireless broadcast model. Highest-ETT edge determines
//    expected completion; other edges finish earlier (bonus deliveries).
// =====================================================================

struct PacketSchedule {
    // Phase 0
    int   ph0_src, ph0_dst;
    double ph0_ett, ph0_cost;

    // Phase 1: forwarder -> relays
    int   ph1_src;
    vector<pair<int,double>> ph1_links; // (relay_node, ETT), sorted by ETT desc
    double ph1_cost;                    // = dp_L1[f][l2_mask]

    // Phase 2: relay -> dest (one entry per relay)
    struct RelayPhase {
        int relay;
        vector<pair<int,double>> links; // (dest_node, ETT), sorted ETT desc
        double cost;                    // = dp_L2[relay][dest_mask]
    };
    vector<RelayPhase> ph2;
    double ph2_parallel_cost;   // max of relay costs (true parallel)
    double ph2_serial_cost;     // sum of relay costs (conservative)

    double total_parallel;      // ph0+ph1+ph2_parallel
    double total_serial;        // ph0+ph1+ph2_serial (= analytical DP cost)
};

PacketSchedule buildSchedule(const Strategy &s) {
    PacketSchedule sch;

    // --- Phase 0 ---
    sch.ph0_src = 0; sch.ph0_dst = s.fwd;
    for (auto &e : adj[0])
        if (e.to == s.fwd) { sch.ph0_ett = e.ett; break; }
    sch.ph0_cost = s.c0;

    // --- Phase 1 ---
    sch.ph1_src = s.fwd;
    for (auto &e : adj[s.fwd]) {
        int b = bitOf_L2(e.to);
        if (b >= 0 && (s.l2_tgt & (1 << b)))
            sch.ph1_links.push_back({e.to, e.ett});
    }
    // Sort by ETT descending (bottleneck first)
    sort(sch.ph1_links.begin(), sch.ph1_links.end(),
         [](auto &a, auto &b){ return a.second > b.second; });
    sch.ph1_cost = s.c1;

    // --- Phase 2 ---
    double maxRelCost = 0;
    for (auto &[relay, dmask] : s.asgn.node_dest) {
        PacketSchedule::RelayPhase rp;
        rp.relay = relay;
        rp.cost  = dp_L2[relay][dmask];
        for (auto &e : adj[relay]) {
            int b3 = bitOf_L3(e.to);
            if (b3 >= 0 && (dmask & (1 << b3)))
                rp.links.push_back({e.to, e.ett});
        }
        sort(rp.links.begin(), rp.links.end(),
             [](auto &a, auto &b){ return a.second > b.second; });
        sch.ph2.push_back(rp);
        maxRelCost = max(maxRelCost, rp.cost);
    }
    sch.ph2_parallel_cost = maxRelCost;
    sch.ph2_serial_cost   = s.c2;

    sch.total_serial   = sch.ph0_cost + sch.ph1_cost + sch.ph2_serial_cost;
    sch.total_parallel = sch.ph0_cost + sch.ph1_cost + maxRelCost;

    return sch;
}

void printPacketStrategy(const Strategy &s) {
    PacketSchedule sch = buildSchedule(s);

    cout << "╔══════════════════════════════════════════════════════════════╗\n";
    cout << "║   PACKET TRANSMISSION STRATEGY — Forwarding Set {"
         << s.fwd << "}           ║\n";
    cout << "╚══════════════════════════════════════════════════════════════╝\n\n";

    // ── PHASE 0 ────────────────────────────────────────────────────────
    cout << "  PHASE 0 │ Source → Forwarder\n";
    cout << "  ─────────────────────────────────────────────────\n";
    cout << "    Action : Node 0 transmits a packet to Node " << sch.ph0_dst << "\n";
    cout << "    Link   : 0 → " << sch.ph0_dst
         << "  (p = " << fixed << setprecision(4) << (1.0/sch.ph0_ett)
         << ",  ETT = " << sch.ph0_ett << " slots)\n";
    cout << "    Rule   : Retry every slot until Node " << sch.ph0_dst << " ACKs.\n";
    cout << "    E[slots]: ";  printVal(sch.ph0_cost); cout << "\n\n";

    // ── PHASE 1 ────────────────────────────────────────────────────────
    cout << "  PHASE 1 │ Forwarder " << sch.ph1_src << " → Relay nodes\n";
    cout << "  ─────────────────────────────────────────────────\n";
    cout << "    Action : Node " << sch.ph1_src
         << " broadcasts to relay set " << maskStr(s.l2_tgt, 4) << "\n";
    cout << "    Rule   : Every slot, all undelivered relays attempt reception\n";
    cout << "             simultaneously (wireless broadcast model).\n";
    cout << "             Stop when ALL relays in set have the packet.\n";
    cout << "    Links (ETT-priority order — highest ETT = bottleneck first):\n";
    int rank = 1;
    for (auto &[nd, ett] : sch.ph1_links) {
        double p = 1.0 / ett;
        cout << "      #" << rank++ << "  Node " << sch.ph1_src << " → " << nd
             << "  p=" << fixed << setprecision(4) << p
             << "  ETT=" << ett << "  ["
             << (rank == 2 ? "BOTTLENECK" : "faster") << "]\n";
    }
    cout << "    E[slots]: ";  printVal(sch.ph1_cost); cout << "\n\n";

    // ── PHASE 2 ────────────────────────────────────────────────────────
    cout << "  PHASE 2 │ Relay nodes → Destinations (PARALLEL)\n";
    cout << "  ─────────────────────────────────────────────────\n";
    cout << "    Action : Each relay independently multicasts to its\n";
    cout << "             assigned destinations. All relays run in parallel.\n";
    cout << "    Rule   : Within each relay, broadcast every slot to all\n";
    cout << "             undelivered destinations until all received.\n\n";

    for (auto &rp : sch.ph2) {
        cout << "    Relay " << rp.relay << " → " << maskStr(s.asgn.node_dest.at(rp.relay), 7) << ":\n";
        cout << "      Links (ETT-priority order):\n";
        int r2 = 1;
        for (auto &[nd, ett] : rp.links) {
            double p = 1.0 / ett;
            cout << "        #" << r2++ << "  Node " << rp.relay << " → " << nd
                 << "  p=" << fixed << setprecision(4) << p
                 << "  ETT=" << ett << "\n";
        }
        cout << "      E[slots for this relay]: "; printVal(rp.cost); cout << "\n\n";
    }

    cout << "    Phase 2 E[slots] (parallel — bottleneck relay): ";
    printVal(sch.ph2_parallel_cost); cout << "\n";
    cout << "    Phase 2 E[slots] (serial conservative model)  : ";
    printVal(sch.ph2_serial_cost); cout << "\n\n";

    // ── TOTAL ──────────────────────────────────────────────────────────
    cout << "  TOTAL EXPECTED SLOTS\n";
    cout << "  ─────────────────────────────────────────────────\n";
    cout << "    [ parallel Phase 2 ]  Ph0 + Ph1 + Ph2(parallel) = ";
    printVal(sch.ph0_cost); cout << " + ";
    printVal(sch.ph1_cost); cout << " + ";
    printVal(sch.ph2_parallel_cost); cout << " = ";
    printVal(sch.total_parallel); cout << "\n";

    cout << "    [ serial/DP model  ]  Ph0 + Ph1 + Ph2(serial)   = ";
    printVal(sch.ph0_cost); cout << " + ";
    printVal(sch.ph1_cost); cout << " + ";
    printVal(sch.ph2_serial_cost); cout << " = ";
    printVal(sch.total_serial); cout << "\n\n";

    // ── MULTIRATE SCHEDULING ADVICE ────────────────────────────────────
    cout << "  MULTIRATE SCHEDULING ADVICE\n";
    cout << "  ─────────────────────────────────────────────────\n";

    // Find bottleneck edge overall
    double worst_ett = 0; int worst_from = -1, worst_to = -1;
    auto checkEdge = [&](int u, int v) {
        for (auto &e : adj[u])
            if (e.to == v && e.ett > worst_ett)
                { worst_ett = e.ett; worst_from = u; worst_to = v; }
    };
    checkEdge(0, s.fwd);
    for (auto &[nd, _] : sch.ph1_links) checkEdge(s.fwd, nd);
    for (auto &rp : sch.ph2)
        for (auto &[nd, __] : rp.links) checkEdge(rp.relay, nd);

    cout << "    Global bottleneck edge : " << worst_from << " → " << worst_to
         << "  (ETT = " << fixed << setprecision(4) << worst_ett << " slots)\n";
    cout << "    → Allocate highest TX power / best modulation to this link.\n";
    cout << "    → In multirate CSMA: assign lowest rate (most robust MCS)\n";
    cout << "      to bottleneck; faster links can use higher MCS.\n";

    // Phase 2 parallel speedup
    double speedup = sch.ph2_serial_cost - sch.ph2_parallel_cost;
    if (speedup > 1e-9) {
        cout << "    → Phase 2 parallelism saves " << fixed << setprecision(4)
             << speedup << " expected slots vs sequential relay operation.\n";
    }
    cout << "\n";
}

// ==================== MAIN ====================
int main() {
    cout << fixed << setprecision(4);

    cout << "============================================================\n";
    cout << " BTP Extended: Correlation-Driven Multicasting (3 Levels)\n";
    cout << " Graph: 0 -> {1,2,3} -> {4,5,6} -> {7,8,9}\n";
    cout << " Destinations: {7, 8, 9}\n";
    cout << "============================================================\n\n";

    // ===================== BUILD GRAPH =====================
    // Helper lambda to build an Edge with auto-computed ETT
    auto E = [](int to, double p) -> Edge {
        return {to, p, (p > 0) ? 1.0 / p : INF};
    };

    // Level 0 -> Level 1
    adj[0] = {E(1, 1.0/3), E(2, 1.0/6), E(3, 1.0/2)};

    // Level 1 -> Level 2
    adj[1] = {E(4, 1.0/3), E(5, 1.0/2)};
    adj[2] = {E(4, 1.0/4), E(5, 1.0/3), E(6, 1.0/3)};
    adj[3] = {E(5, 1.0/3), E(6, 1.0/4)};

    // Level 2 -> Level 3
    adj[4] = {E(7, 1.0/3), E(8, 1.0/2)};
    adj[5] = {E(8, 1.0/3), E(9, 1.0/4)};
    adj[6] = {E(9, 1.0/3)};

    // Print graph with ETT per edge
    cout << "--- GRAPH EDGES & ETT (Expected Transmissions per Edge) ---\n";
    cout << "+--------+--------+----------+----------+\n";
    cout << "| From   | To     |    p     |  ETT=1/p |\n";
    cout << "+--------+--------+----------+----------+\n";
    for (int u = 0; u <= 6; u++) {
        for (auto &e : adj[u]) {
            printf("| %3d    | %3d    | %8.4f | %8.4f |\n",
                   u, e.to, e.p, e.ett);
        }
    }
    cout << "+--------+--------+----------+----------+\n\n";

    // ===================== ETT SUMMARY TABLE =====================
    cout << "--- ETT (Expected Time per Transmission) per Edge ---\n";
    cout << "  ETT(u->v) = 1/p  [expected single-link trials before success]\n\n";

    auto printETTLevel = [&](const string &label, int from_lo, int from_hi) {
        cout << "  " << label << ":\n";
        for (int u = from_lo; u <= from_hi; u++) {
            for (auto &e : adj[u]) {
                cout << "    ETT(" << u << "->" << e.to << ") = 1/"
                     << fixed << setprecision(4) << e.p << " = ";
                printVal(e.ett);
                cout << " tx\n";
            }
        }
        cout << "\n";
    };

    printETTLevel("Level 0->1  (source -> forwarders)",    0, 0);
    printETTLevel("Level 1->2  (forwarders -> relays)",    1, 3);
    printETTLevel("Level 2->3  (relays -> destinations)",  4, 6);

    // ===================== BASE CASES =====================
    // Destinations 7,8,9: self-delivery is free
    for (int u = 7; u <= 9; u++)
        for (int mask = 0; mask < 8; mask++) {
            if (mask == 0)
                dp_L2[u][mask] = 0;
            else {
                int b = bitOf_L3(u);
                dp_L2[u][mask] = (mask == (1 << b)) ? 0 : INF;
            }
        }

    // ===================== LEVEL 2 DP =====================
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

    // ===================== LEVEL 1 DP =====================
    cout << "--- LEVEL 1 DP: nodes {1,2,3} -> level-2 {4,5,6} ---\n";
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

    // ===================== FORWARDING SET ANALYSIS =====================
    cout << "--- FORWARDING SET ANALYSIS (Heuristic: Superset Pruning) ---\n";

    vector<int> valid_sets, minimal_sets;

    for (int mask = 1; mask < 8; mask++) {
        // Level-2 nodes reachable from this forwarding set
        int l2_reach = 0;
        for (int b = 0; b < 3; b++) {
            if (!(mask & (1 << b))) continue;
            int node = b + 1;
            for (auto &e : adj[node])
                if (bitOf_L2(e.to) >= 0)
                    l2_reach |= (1 << bitOf_L2(e.to));
        }

        // Level-3 destinations reachable via those level-2 nodes
        int l3_reach = 0;
        for (int b = 0; b < 3; b++)
            if (l2_reach & (1 << b))
                for (auto &e : adj[b + 4])
                    if (bitOf_L3(e.to) >= 0)
                        l3_reach |= (1 << bitOf_L3(e.to));

        cout << "  Set " << maskStr(mask, 1) << ": L2 reach = "
             << maskStr(l2_reach, 4) << ", L3 reach = "
             << maskStr(l3_reach, 7);

        if (l3_reach == 7) {
            cout << "  -> VALID\n";
            valid_sets.push_back(mask);
        } else {
            cout << "  -> INVALID (can't reach all destinations)\n";
        }
    }

    // Remove supersets
    for (int s : valid_sets) {
        bool is_super = false;
        for (int t : valid_sets)
            if (t != s && (t & s) == t) { is_super = true; break; }
        if (!is_super) minimal_sets.push_back(s);
    }

    cout << "\nMinimal forwarding sets (after superset pruning):\n";
    for (int m : minimal_sets)
        cout << "  " << maskStr(m, 1) << "\n";
    cout << "\n";

    // ===================== COST COMPUTATION =====================
    cout << "--- COST COMPUTATION FOR MINIMAL SETS ---\n\n";

    vector<Strategy> strategies;

    for (int m : minimal_sets) {
        if (__builtin_popcount(m) == 1) {
            int b = __builtin_ctz(m);
            int f = b + 1;

            Strategy s = bestStrategy(f);
            strategies.push_back(s);

            cout << "Forwarding Set {" << f << "}:\n";

            // ETT for the source->forwarder edge
            double ett_0f = 0;
            for (auto &e : adj[0])
                if (e.to == f) { ett_0f = e.ett; break; }
            cout << "  Hop 0->" << f << ":  p = ";
            for (auto &e : adj[0]) if (e.to == f)
                cout << fixed << setprecision(4) << e.p;
            cout << ",  ETT = " << fixed << setprecision(4) << ett_0f
                 << " tx,  Cost(0->" << f << ") = ";
            printVal(s.c0);
            cout << "\n";

            cout << "  Hop " << f << "->L2 (multicast):  L2 target = "
                 << maskStr(s.l2_tgt, 4) << "\n";
            // Print ETT for each edge from f to its L2 targets
            for (auto &e : adj[f]) {
                if (bitOf_L2(e.to) >= 0 && (s.l2_tgt & (1 << bitOf_L2(e.to)))) {
                    cout << "    ETT(" << f << "->" << e.to << ") = "
                         << fixed << setprecision(4) << e.ett << " tx\n";
                }
            }
            cout << "  Cost(L1 multicast) = ";
            printVal(s.c1);
            cout << "\n";

            cout << "  Hop L2->L3 (destination delivery):\n";
            for (auto &[node, dmask] : s.asgn.node_dest) {
                // Print ETT for each edge from this L2 node to its assigned dests
                for (auto &e : adj[node]) {
                    int b3 = bitOf_L3(e.to);
                    if (b3 >= 0 && (dmask & (1 << b3))) {
                        cout << "    ETT(" << node << "->" << e.to << ") = "
                             << fixed << setprecision(4) << e.ett << " tx\n";
                    }
                }
                cout << "    Node " << node << " -> "
                     << maskStr(dmask, 7) << ", Cost(L2->L3) = ";
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

    // ===================== BEST STRATEGY =====================
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

    // ===================== PACKET TRANSMISSION STRATEGY =====================
    cout << "============================================================\n";
    cout << " PACKET TRANSMISSION STRATEGY (all minimal forwarding sets)\n";
    cout << "============================================================\n\n";
    for (auto &s : strategies)
        printPacketStrategy(s);

    if (best) {
        cout << "============================================================\n";
        cout << " RECOMMENDED STRATEGY: Forwarding Set {" << best->fwd << "}\n";
        cout << " Reason: Lowest total expected transmissions = ";
        printVal(best->cost);
        cout << "\n";
        cout << "============================================================\n\n";
    }

    // ===================== SIMULATION =====================
    cout << "--- SIMULATION VERIFICATION ---\n";
    int runs = 100000;

    for (auto &s : strategies) {
        double totCons = 0, totOppo = 0;
        for (int i = 0; i < runs; i++) {
            totCons += simConservative(s);
            totOppo += simOpportunistic(s);
        }
        double avgCons = totCons / runs;
        double avgOppo = totOppo / runs;

        cout << "Forwarding Set {" << s.fwd << "}:\n";
        cout << "  Analytical (DP)     = ";
        printVal(s.cost);
        cout << "\n";
        cout << "  Sim Conservative    = " << avgCons
             << "  (" << runs << " runs)\n";
        cout << "  Sim Opportunistic   = " << avgOppo
             << "  (" << runs << " runs)\n";
        cout << "\n";
    }

    // ===================== SUMMARY TABLE =====================
    cout << "--- SUMMARY (ETT-based costs: expected transmissions to deliver all destinations) ---\n";
    cout << "+----------+-----------+----------------+------------------+\n";
    cout << "| Fwd Set  | Analytical| Sim(Conserv.)  | Sim(Opportun.)   |\n";
    cout << "+----------+-----------+----------------+------------------+\n";

    for (auto &s : strategies) {
        double totC = 0, totO = 0;
        // Use same RNG seed for consistency? No, already ran above.
        // Re-run quickly for table (smaller run count)
        int tRuns = 50000;
        for (int i = 0; i < tRuns; i++) {
            totC += simConservative(s);
            totO += simOpportunistic(s);
        }
        printf("|   {%d}    | %9.4f | %14.4f | %16.4f |\n",
               s.fwd, s.cost, totC / tRuns, totO / tRuns);
    }
    cout << "+----------+-----------+----------------+------------------+\n";

    return 0;
}