/*
 * Header file containing functions used in other simulations.
 */

// A Genotype is 3 ints: each int labels a node on the underlying graph.
typedef tuple<int,int,int> Genotype;
typedef Genotype State;

// A History is a sequence of States
typedef vector<State> History;

// An Ensemble is a map from Genotypes to integers (populations)
// In the notation on my blackboard and notes, the Ensemble is the set of 
// occupation/population numbers for each genotype g, denoted \{n_g\}
typedef map<Genotype,uint64_t> Ensemble;

// A PathEnsemble is the equivalent for the path graph: a map from a path (a
// History) to an integer population.
typedef map<History,uint64_t> PathEnsemble;

// Initialisation functions:
void init_genotypes(set<State> &q)
{
    for (int i=0; i<5; i++) {
        for (int j=0; j<5; j++) {
            for (int k=0; k<2; k++) {
                State p = make_tuple(i,j,k);
                q.emplace(p);
            }
        }
    }   
}

void init_ends(set<State> &q)
{
    for (int i=3; i<5; i++) {
        for (int j=3; j<5; j++) {
            State p = make_tuple(i,j,1);
            q.emplace(p);
        }
    }
}

void init_ensemble(Ensemble &q, set<State> G)
{
    for (auto g = G.begin(); g != G.end(); ++g)
    {
        // for genotype in set of genotypes, set all genotype populations to zero
        q[*g] = 0; 
    } 
}

void reset_ensemble(Ensemble &q)
{
    for (int i=0; i<5; i++) {
        for (int j=0; j<5; j++) {
            for (int k=0; k<2; k++) {
                q[make_tuple(i,j,k)]=0;
            }
        }
    }
}

void reset_path_ensemble(PathEnsemble &popln, set<History> L)
{
    for (auto lineage = L.begin(); lineage != L.end(); ++lineage)
    {
        popln[*lineage] = 0;
    }
}

long double flow(State p, State q) {
    int i1 = get<0>(p);
    int i2 = get<0>(q);
    int j1 = get<1>(p);
    int j2 = get<1>(q);
    int k1 = get<2>(p);
    int k2 = get<2>(q);

    if ((rates[i1][i2] > 0)&&(j1==j2)&&(k1==k2)) {
        return rates[i1][i2];
    } else if ((i1==i2)&&(rates2[j1][j2] > 0)&&(k1==k2)) {
        return rates2[j1][j2];
    } else if ((i1==i2)&&(j1==j2)&&(rates3[k1][k2] > 0)) {
        return rates3[k1][k2];
    } else {
        return 0;
    }
}

void trace_paths(set<History> &q, set<State> genotypes)
{
    // unit test for lineage-tracking mode: populates the set &q with all
    // possible lineages (paths on the graph of genotypes G)
    History root;
    root.push_back(make_tuple(0,0,0));
    q.emplace(root);

    int newnodes = 1;
    while (newnodes==1) {
        size_t oldsize = q.size();

        for (auto l = q.begin(); l != q.end(); ++l) {
            // For each lineage l, iterate over nodes g in genotypes.
            for (auto g = genotypes.begin(); g!= genotypes.end(); ++g) {
                State last = (*l).back();
                // If l.back() is connected to g, form a new history:
                if (flow(last,(*g))>0) {
                    History h = (*l);
                    h.push_back(*g);
                    // and add h to lineages
                    q.emplace(h);
                }
            }
        }

        // Do this until no new nodes are added
        if (q.size()==oldsize) {
            newnodes = 0;
        } else {
            newnodes = 1;
        }
    }
}

map<State,set<State>> generate_neighbours(set<State> G) {
    // find and return a map from G to subsets of G: the map of neighbouring
    // sites
    map<State,set<State>> Neighbours;
    for (auto p = G.begin(); p != G.end(); ++p) {
        // calculate possible neighbours of p and add them to a set
        set<State> nbs;

        for (auto q = G.begin(); q != G.end(); ++q) {
            // if q is a neighbour of p, add it to this set
            if (flow(*p,*q)>0) nbs.emplace(*q);
        }
        Neighbours[*p] = nbs;
    }

    return Neighbours;
}
