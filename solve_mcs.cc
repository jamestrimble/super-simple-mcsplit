#include "graph.hh"
#include "solve_mcs.hh"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

std::atomic<bool> abort_due_to_timeout;

namespace {
    struct VarAssignment
    {
        Assignment assignment;
    };

    struct VarAssignments
    {
        vector<VarAssignment> var_assignments;

        auto push(VarAssignment a) -> void
        {
            var_assignments.push_back(a);
        }

        auto pop() -> void
        {
            var_assignments.pop_back();
        }

        auto get_var_assignments() const -> const vector<VarAssignment> &
        {
            return var_assignments;
        }

        auto get_num_vtx_assignments() -> unsigned
        {
            return var_assignments.size();
        }

        auto clear() -> void
        {
            var_assignments.clear();
        }
    };

    struct Bidomain
    {
        int l,        r;        // start indices of left and right sets
        int left_len, right_len;
        bool is_adjacent;
    };

    auto calculate_degrees(const Graph & g) -> vector<int>
    {
        vector<int> degree(g.n, 0);
        for (int v=0; v<g.n; v++) {
            for (int w=0; w<g.n; w++) {
                unsigned int mask = 0xFFFFu;
                if (g.adjmat[v][w] & mask) degree[v]++;
                if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
            }
        }
        return degree;
    }

    class MCS
    {
        const Graph & g0;
        const Graph & g1;
        const Params params;
        McsStats stats;
        vector<int> left;
        vector<int> right;
        vector<Assignment> incumbent;

        auto show(const vector<VarAssignment>& current, const vector<Bidomain> &domains) -> void
        {
            cout << "Nodes: " << stats.nodes << endl;
            cout << "Length of current assignment: " << current.size() << endl;
            cout << "Current assignment:";
            for (unsigned int i=0; i<current.size(); i++) {
                cout << "  (" << current[i].assignment.v << " -> " << current[i].assignment.w << ")";
            }
            cout << endl;
            for (unsigned int i=0; i<domains.size(); i++) {
                struct Bidomain bd = domains[i];
                cout << "Left  ";
                for (int j=0; j<bd.left_len; j++)
                    cout << left[bd.l + j] << " ";
                cout << endl;
                cout << "Right  ";
                for (int j=0; j<bd.right_len; j++)
                    cout << right[bd.r + j] << " ";
                cout << endl;
            }
            cout << "\n" << endl;
        }

        auto calc_bound(const vector<Bidomain>& domains) -> int
        {
            int bound = 0;
            for (const Bidomain &bd : domains)
                bound += std::min(bd.left_len, bd.right_len);
            return bound;
        }

        // precondition: len > 0
        auto find_min_value(const vector<int>& arr, int start_idx, int len) -> int
        {
            return *std::min_element(arr.begin() + start_idx, arr.begin() + start_idx + len);
        }

        auto select_bidomain(const vector<Bidomain>& domains, int current_matching_size) -> int
        {
            // Select the bidomain with the smallest max(leftsize, rightsize), breaking
            // ties on the smallest vertex index in the left set
            std::pair<int, int> best_score { std::numeric_limits<int>::max(), std::numeric_limits<int>::max() };
            int best = -1;
            for (unsigned int i=0; i<domains.size(); i++) {
                const Bidomain &bd = domains[i];
                if (params.connected && current_matching_size>0 && !bd.is_adjacent) continue;
                int len = params.heuristic == Heuristic::min_max ?
                        std::max(bd.left_len, bd.right_len) :
                        bd.left_len * bd.right_len;
                int tie_breaker = find_min_value(left, bd.l, bd.left_len);
                auto score = std::make_pair( len, tie_breaker );
                if (score < best_score) {
                    best_score = score;
                    best = i;
                }
            }
            return best;
        }

        // Returns length of left half of array
        auto partition(vector<int>& all_vv, int start, int len, const vector<unsigned int> & adjrow) -> int
        {
            auto it = std::partition(
                    all_vv.begin() + start,
                    all_vv.begin() + start + len,
                    [&](const int elem){ return 0 != adjrow[elem]; });
            return it - (all_vv.begin() + start);
        }

        // multiway is for directed and/or labelled graphs
        auto filter_domains(const vector<Bidomain> & d, int v, int w, bool multiway) -> vector<Bidomain>
        {
            vector<Bidomain> new_d;
            new_d.reserve(d.size());
            for (const Bidomain &old_bd : d) {
                int l = old_bd.l;
                int r = old_bd.r;
                // After these two partitions, left_len and right_len are the lengths of the
                // arrays of vertices with edges from v or w (int the directed case, edges
                // either from or to v or w)
                int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
                int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
                int left_len_noedge = old_bd.left_len - left_len;
                int right_len_noedge = old_bd.right_len - right_len;
                if (left_len_noedge && right_len_noedge)
                    new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
                if (multiway && left_len && right_len) {
                    auto& adjrow_v = g0.adjmat[v];
                    auto& adjrow_w = g1.adjmat[w];
                    auto l_begin = std::begin(left) + l;
                    auto r_begin = std::begin(right) + r;
                    std::sort(l_begin, l_begin+left_len, [&](int a, int b)
                            { return adjrow_v[a] < adjrow_v[b]; });
                    std::sort(r_begin, r_begin+right_len, [&](int a, int b)
                            { return adjrow_w[a] < adjrow_w[b]; });
                    int l_top = l + left_len;
                    int r_top = r + right_len;
                    while (l<l_top && r<r_top) {
                        unsigned int left_label = adjrow_v[left[l]];
                        unsigned int right_label = adjrow_w[right[r]];
                        if (left_label < right_label) {
                            l++;
                        } else if (left_label > right_label) {
                            r++;
                        } else {
                            int lmin = l;
                            int rmin = r;
                            do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
                            do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
                            new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                        }
                    }
                } else if (left_len && right_len) {
                    new_d.push_back({l, r, left_len, right_len, true});
                }
            }
            return new_d;
        }

        auto remove_vtx_from_left_domain(Bidomain& bd, int v) -> void
        {
            int i = 0;
            while(left[bd.l + i] != v) i++;
            std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
            bd.left_len--;
        }

        auto remove_bidomain(vector<Bidomain>& domains, int idx) -> void
        {
            domains[idx] = domains[domains.size()-1];
            domains.pop_back();
        }

        auto update_incumbent(VarAssignments & current) -> void
        {
            incumbent.clear();
            for (auto a : current.get_var_assignments())
                if (a.assignment.w != -1)
                    incumbent.push_back(a.assignment);

            stats.time_of_last_incumbent_update = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - params.start_time).count();

            if (!params.quiet) cout << "New incumbent " << incumbent.size() << " at time " <<
                    stats.time_of_last_incumbent_update << " ms" << endl;
        }

        enum class Search
        {
            Aborted,
            Done,
            Restart,
            ReachedGoal
        };

        auto is_bound_minimum_possible(unsigned bound, unsigned matching_size_goal) -> bool
        {
            return bound == incumbent.size() + 1 ||
                    (params.mcsplit_down && bound == matching_size_goal);
        }

        auto search(
                VarAssignments & current,
                vector<Bidomain> & domains,
                unsigned int matching_size_goal) -> Search
        {
            if (abort_due_to_timeout)
                return Search::Aborted;

            ++stats.nodes;

            if (params.verbose)
                show(current.get_var_assignments(), domains);

            if (current.get_num_vtx_assignments() > incumbent.size()) {
                update_incumbent(current);
                if (incumbent.size() == matching_size_goal)
                    return Search::ReachedGoal;
            }

            unsigned int bound = current.get_num_vtx_assignments() + calc_bound(domains);
            if (bound <= incumbent.size() || (params.mcsplit_down && bound < matching_size_goal))
                return Search::Done;

            int bd_idx = select_bidomain(domains, current.get_num_vtx_assignments());
            if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
                return Search::Done;
            Bidomain &bd = domains[bd_idx];

            int v = find_min_value(left, bd.l, bd.left_len);

            // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
            std::vector<int> possible_values(right.begin() + bd.r,  // the vertices in the colour class beginning at bd.r
                    right.begin() + bd.r + bd.right_len);
            std::sort(possible_values.begin(), possible_values.end());

            if (!is_bound_minimum_possible(bound, matching_size_goal) || bd.left_len > bd.right_len)
                possible_values.push_back(-1);

            remove_vtx_from_left_domain(domains[bd_idx], v);
            bd.right_len--;

            for (int w : possible_values) {
                Search search_result;
                if (w != -1) {
                    current.push({{v, w}});
                    // swap w to the end of its colour class
                    auto it = std::find(right.begin() + bd.r, right.end(), w);
                    *it = right[bd.r + bd.right_len];
                    right[bd.r + bd.right_len] = w;

                    auto new_domains = filter_domains(domains, v, w, params.directed || params.edge_labelled);
                    search_result = search(current, new_domains, matching_size_goal);
                    current.pop();
                } else {
                    bd.right_len++;
                    if (bd.left_len == 0)
                        remove_bidomain(domains, bd_idx);
                    search_result = search(current, domains, matching_size_goal);
                }
                switch (search_result)
                {
                case Search::Restart:     break;   // this should never happen
                case Search::Aborted:     return Search::Aborted;
                case Search::ReachedGoal: return Search::ReachedGoal;
                case Search::Done:        break;
                }
            }

            return Search::Done;
        }

        auto run_search(vector<Bidomain> domains, unsigned int matching_size_goal) -> void
        {
            VarAssignments current;
            search(current, domains, matching_size_goal);
        }

    public:
        MCS(Graph & g0, Graph & g1, Params params)
            : g0(g0), g1(g1), params(params)
        { }

        auto run() -> std::pair<vector<Assignment>, McsStats>
        {
            auto domains = vector<Bidomain> {};

            std::set<unsigned int> left_labels;
            std::set<unsigned int> right_labels;
            for (unsigned int label : g0.label) left_labels.insert(label);
            for (unsigned int label : g1.label) right_labels.insert(label);
            std::set<unsigned int> labels;  // labels that appear in both graphs
            std::set_intersection(std::begin(left_labels),
                                  std::end(left_labels),
                                  std::begin(right_labels),
                                  std::end(right_labels),
                                  std::inserter(labels, std::begin(labels)));

            // Create a bidomain for each label that appears in both graphs
            for (unsigned int label : labels) {
                int start_l = left.size();
                int start_r = right.size();

                for (int i=0; i<g0.n; i++)
                    if (g0.label[i]==label)
                        left.push_back(i);
                for (int i=0; i<g1.n; i++)
                    if (g1.label[i]==label)
                        right.push_back(i);

                int left_len = left.size() - start_l;
                int right_len = right.size() - start_r;
                domains.push_back({start_l, start_r, left_len, right_len, false});
            }


	    if (params.mcsplit_down) {
		for (unsigned int goal = std::min(g0.n, g1.n) ; goal > 0 ; --goal) {
		    if (incumbent.size() == goal) break;
		    run_search(domains, goal);
		    if (incumbent.size() == goal || abort_due_to_timeout) break;
		    if (!params.quiet) cout << "Upper bound: " << goal-1 << std::endl;
		}
	    } else {
                run_search(domains, std::min(g0.n, g1.n));
	    }

            return {incumbent, stats};
        }
    };

    auto sum(const vector<int> & vec) -> int
    {
        return std::accumulate(std::begin(vec), std::end(vec), 0);
    }
};

auto solve_mcs(Graph & g0, Graph & g1, Params params)
		-> std::pair<vector<Assignment>, McsStats>
{
    auto g0_deg = calculate_degrees(g0);
    auto g1_deg = calculate_degrees(g1);

    // As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) > g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) > g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
    });

    struct Graph g0_sorted = induced_subgraph(g0, vv0);
    struct Graph g1_sorted = induced_subgraph(g1, vv1);

    auto solution = MCS(g0_sorted, g1_sorted, params).run();

    // Convert to indices from original, unsorted graphs
    for (auto& vtx_pair : solution.first) {
        vtx_pair.v = vv0[vtx_pair.v];
        vtx_pair.w = vv1[vtx_pair.w];
    }

    return solution;
}
