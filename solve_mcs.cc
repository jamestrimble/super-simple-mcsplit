#include "graph.hh"
#include "solve_mcs.hh"

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

std::atomic<bool> abort_due_to_timeout;

namespace {
    struct Bidomain
    {
        // start indices of left and right sets
        int l;
        int r;

        int left_len;
        int right_len;
    };

    struct DomainStore
    {
        vector<int> left;
        vector<int> right;
        vector<Bidomain> domains;
    };

    class MCS
    {
        const Graph & g0;
        const Graph & g1;
        const Params params;
        McsStats stats;
        vector<Assignment> incumbent;

        auto print_state(const vector<Assignment>& current, const DomainStore &domain_store) -> void
        {
            cout << "Nodes: " << stats.nodes << endl;
            cout << "Length of current assignment: " << current.size() << endl;
            cout << "Current assignment:";
            for (unsigned int i=0; i<current.size(); i++) {
                cout << "  (" << current[i].v << " -> " << current[i].w << ")";
            }
            cout << endl;
            for (unsigned int i=0; i<domain_store.domains.size(); i++) {
                struct Bidomain bd = domain_store.domains[i];
                cout << "Left  ";
                for (int j=0; j<bd.left_len; j++)
                    cout << domain_store.left[bd.l + j] << " ";
                cout << endl;
                cout << "Right  ";
                for (int j=0; j<bd.right_len; j++)
                    cout << domain_store.right[bd.r + j] << " ";
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

        auto select_bidomain(const DomainStore & domain_store, int current_matching_size) -> int
        {
            // Select the bidomain with the smallest max(leftsize, rightsize), breaking
            // ties on the smallest vertex index in the left set
            std::pair<int, int> best_score { std::numeric_limits<int>::max(), std::numeric_limits<int>::max() };
            int best = -1;
            for (unsigned int i=0; i<domain_store.domains.size(); i++) {
                const Bidomain &bd = domain_store.domains[i];
                int len = params.heuristic == Heuristic::min_max ?
                        std::max(bd.left_len, bd.right_len) :
                        bd.left_len * bd.right_len;
                int tie_breaker = find_min_value(domain_store.left, bd.l, bd.left_len);
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

        auto refined_domains(const DomainStore & ds, int v, int w) -> DomainStore
        {
            DomainStore new_ds { ds.left, ds.right, {} };
            new_ds.domains.reserve(ds.domains.size());
            for (const Bidomain &old_bd : ds.domains) {
                int l = old_bd.l;
                int r = old_bd.r;
                // After these two partitions, left_len and right_len are the lengths of the
                // arrays of vertices with edges from v or w
                int left_len = partition(new_ds.left, l, old_bd.left_len, g0.adjmat[v]);
                int right_len = partition(new_ds.right, r, old_bd.right_len, g1.adjmat[w]);
                int left_len_noedge = old_bd.left_len - left_len;
                int right_len_noedge = old_bd.right_len - right_len;
                if (left_len_noedge && right_len_noedge)
                    new_ds.domains.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge});
                if (left_len && right_len)
                    new_ds.domains.push_back({l, r, left_len, right_len});
            }
            return new_ds;
        }

        auto remove_vtx_from_left_domain(vector<int> & left, Bidomain & bd, int v) -> void
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

        auto update_incumbent(vector<Assignment> & current) -> void
        {
            incumbent = current;

            stats.time_of_last_incumbent_update = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - params.start_time).count();

            if (!params.quiet) cout << "New incumbent " << incumbent.size() << " at time " <<
                    stats.time_of_last_incumbent_update << " ms" << endl;
        }

        enum class Search
        {
            Aborted,
            Done,
            ReachedGoal
        };

        auto search(
                vector<Assignment> & current,
                DomainStore & domain_store,
                unsigned int matching_size_goal) -> Search
        {
            if (abort_due_to_timeout)
                return Search::Aborted;

            ++stats.nodes;

            if (params.verbose)
                print_state(current, domain_store);

            if (current.size() > incumbent.size()) {
                update_incumbent(current);
                if (incumbent.size() == matching_size_goal)
                    return Search::ReachedGoal;
            }

            unsigned int bound = current.size() + calc_bound(domain_store.domains);
            if (bound <= incumbent.size() || (params.mcsplit_down && bound < matching_size_goal))
                return Search::Done;

            int bd_idx = select_bidomain(domain_store, current.size());
            Bidomain &bd = domain_store.domains[bd_idx];

            int v = find_min_value(domain_store.left, bd.l, bd.left_len);

            // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
            auto right_label_class_begin = domain_store.right.begin() + bd.r;
            auto right_label_class_end = right_label_class_begin + bd.right_len;
            int & right_label_class_last_elem = *std::prev(right_label_class_end);
            std::sort(right_label_class_begin, right_label_class_end);

            remove_vtx_from_left_domain(domain_store.left, bd, v);
            bd.right_len--;

            for (auto it=right_label_class_begin; it!=right_label_class_end; ++it) {
                int w = *it;
                current.push_back({v, w});
                std::swap(*it, right_label_class_last_elem); // swap w to the end of its label class

                auto new_domain_store = refined_domains(domain_store, v, w);
                Search search_result = search(current, new_domain_store, matching_size_goal);
                current.pop_back();
                std::swap(*it, right_label_class_last_elem); // swap w back to its correct place in the sorted label class
                switch (search_result)
                {
                case Search::Aborted:     return Search::Aborted;
                case Search::ReachedGoal: return Search::ReachedGoal;
                case Search::Done:        break;
                }
            }

            // Try assigning v to \bot (i.e. not using it in our mapping)
            bd.right_len++;
            if (bd.left_len == 0)
                remove_bidomain(domain_store.domains, bd_idx);
            return search(current, domain_store, matching_size_goal);
        }

    public:
        MCS(Graph & g0, Graph & g1, Params params)
            : g0(g0), g1(g1), params(params)
        { }

        auto run() -> std::pair<vector<Assignment>, McsStats>
        {
            DomainStore domain_store;

            for (int i=0; i<g0.n; i++)
                domain_store.left.push_back(i);
            for (int i=0; i<g1.n; i++)
                domain_store.right.push_back(i);

            domain_store.domains.push_back({0, 0, g0.n, g1.n});

            vector<Assignment> current;

	    if (params.mcsplit_down) {
		for (unsigned int goal = std::min(g0.n, g1.n) ; goal > 0 ; --goal) {
		    if (incumbent.size() == goal) break;
                    auto domain_store_copy = domain_store;
                    search(current, domain_store_copy, goal);
		    if (incumbent.size() == goal || abort_due_to_timeout) break;
		    if (!params.quiet) cout << "Upper bound: " << goal-1 << std::endl;
                    cout << stats.nodes << std::endl;
		}
	    } else {
                search(current, domain_store, std::min(g0.n, g1.n));
	    }

            return {incumbent, stats};
        }
    };
};

auto vertices_sorted_by_degree(Graph & g) -> vector<int>
{
    auto deg = calculate_degrees(g);
    vector<int> vv(g.n);
    std::iota(std::begin(vv), std::end(vv), 0);
    std::stable_sort(std::begin(vv), std::end(vv), [&](int a, int b) {
        return deg[a] > deg[b];
    });
    return vv;
}

auto solve_mcs(Graph & g0, Graph & g1, Params params)
		-> std::pair<vector<Assignment>, McsStats>
{
    auto vv0 = vertices_sorted_by_degree(g0);
    auto vv1 = vertices_sorted_by_degree(g1);

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
