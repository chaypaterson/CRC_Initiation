#include <iostream>
#include <map>
#include <vector>
#include <tuple>
#include <set>

/* Enumerate
 * Unit test program to enumerate all possible paths on a directed graph
 * defined by the transition rate matrices in params.h
 */

#include "params.h"
using namespace std;

// Definition of the derived graph:
// Two nodes q and q', each vectors of nodes on the underlying graph,
// are connected iff K(q->last -> q'->last) != 0

// Probability flows between nodes if one integer is different at rate
// rates[i][i'] or rates2[j][j'] or rates3[k][k']
// (which may be zero)
//
// In the graph G', each node is a sequence of nodes in G.

typedef tuple<int,int,int> State;

typedef vector<State> History;

int distance(State p, State q) {
    int d = 0;

    d += abs(get<0>(p)-get<0>(q));
    d += abs(get<1>(p)-get<1>(q));
    d += abs(get<2>(p)-get<2>(q));

    return d;
}

int connected(State p, State q) {
    int i1 = get<0>(p);
    int i2 = get<0>(q);
    int j1 = get<1>(p);
    int j2 = get<1>(q);
    int k1 = get<2>(p);
    int k2 = get<2>(q);

    if ((rates[i1][i2] > 0)&&(j1==j2)&&(k1==k2)) {
        return 1;
    } else if ((i1==i2)&&(rates2[j1][j2] > 0)&&(k1==k2)) {
        return 1;
    } else if ((i1==i2)&&(j1==j2)&&(rates3[k1][k2] > 0)) {
        return 1;
    } else {
        return 0;
    }
}

int main()
{
    // A graph is a set with a map from the set to itself which gives the
    // edges. We use the rate matrices + von Neumann neighbours to determine
    // the connectivity: if vN neighbours + rate != 0, then neighbours.

    set<State> genotypes;

    // populate genotypes with possible genotypes
    // We have a set of nodes: these can be indexed by tuples
    // {[0-4],[0-4],[0-1]}

    for (int i=0; i<5; i++) {
        for (int j=0; j<5; j++) {
            for (int k=0; k<2; k++) {
                State p = make_tuple(i,j,k);
                genotypes.emplace(p);
            }
        }
    }

    set<State> endnodes;
    for (int i=3; i<5; i++) {
        for (int j=3; j<5; j++) {
            State p = make_tuple(i,j,1);
            endnodes.emplace(p);
        }
    }
    cout << genotypes.size() << "," << endnodes.size() << endl;

    // iteratively build the lineages graph:

    set<History> lineages;

    // the root node is in lineages
    History root;
    root.push_back(make_tuple(0,0,0));
    lineages.emplace(root);

    int newnodes = 1;

    while (newnodes==1) {
        newnodes = 0;
        int oldsize = lineages.size();
        // Loop over each lineage l currently in lineages:
        for (auto l=lineages.begin(); l != lineages.end(); ++l) {
            // For each lineage l, iterate over nodes g in genotypes.
            for (auto g=genotypes.begin(); g != genotypes.end(); ++g) {
                State last = (*l).back();
                // If l.back() is connected to g, form a new history: 
                if (connected(last,(*g))==1) {
                    History h = (*l); 
                    h.push_back(*g);
                    // and add h to lineages
                    lineages.emplace(h);
                }
            }
        }

        // Do this until no new nodes are added
        if (lineages.size()==oldsize) {
            newnodes = 0;
        } else {
            newnodes = 1;
        }
    }

    for (auto p =lineages.begin(); p != lineages.end(); ++p) {
    	if (endnodes.count((*p).back())>0) {
	    for (auto q = (*p).begin(); q != (*p).end(); ++q) {
	        cout << get<0>(*q);
		cout << get<1>(*q);
		cout << get<2>(*q);
		if ((*q) != (*p).back()) cout << "-";
	    }
            cout << endl;
        }
    }
}
