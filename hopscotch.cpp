#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>
#include <set>
// Multithreading
#include <thread>
// Thread monitoring
#include <atomic>
#include <chrono>

/* HOPSCOTCH
 * This simulates a compound stochastic process.
 * Some beans can hop on a hopscotch court: nodes on a graph.
 * The rate of each hop is independent and tunable.
 *
 * We want to know the probability that if we have N such beans, at least one 
 * bean will have reached the end of the hopscotch court (and how it did so). 
 * To find this, we simulate many replicates of a game with N beans, and 
 * count how many beans reached the end of the court and had different life 
 * histories.
 * 
 * Each state countains a number of beans. Beans can hop between different
 * states with different rates.
 * */

// Example game: court has three dimensions, each dimension is a graph.
// One dimension is a simple one-step process, the others are two-step
// processes. Each two step process has three possible states in the first
// layer, and two possible states in the end layer. Each state has possible
// "next" states in the layer above, and associated rates for these.



// User-installed headers
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// function to return random integers from a RNG_t and a dist.
// NB: do NOT copy r

// Load parameters:
#include "params.h" // TODO load from XML file or something
using namespace std;

#include "initialisers.h"

struct simConstants {
    double ttmax;
    double dt;
    int runs;
    set<State> genotypes;
    map<State,set<State>> Neighbours;
    int mode;
};

void gillespie_beans (vector<Ensemble> *tResults, int seed, 
        atomic<int> *progress, simConstants Params) {
    // Pull in constant parameters
    double ttmax = Params.ttmax;
    int runs = Params.runs;
    set<State> G = Params.genotypes;
    map<State,set<State>> Neighbours = Params.Neighbours;

    // the main object being manipulated is a mapping of the state to the number
    // of beans in that state: this is the population.
    // declare core objects, initialize inside.

    vector<vector<vector<long double>>> 
        weights(5,vector<vector<long double>>(5,vector<long double>(2,0)));

    // precompute weights
    for (int i=0; i<5; i++) {
        for (int j=0; j<5; j++) {
            for (int k=0; k<2; k++) {
                State p = make_tuple(i,j,k);

                for (auto q = G.begin(); q != G.end(); ++q) {
                    weights[i][j][k] += flow(p,*q);
                }

                // crypt fission propensities
                if (i>2) {
                    if (k==0) {
                        weights[i][j][k] += rate_APClost;
                    } else {
                        weights[i][j][k] += rate_BOTHlost;
                    }
                } else {
                    if (k>0) {
                        weights[i][j][k] += rate_RASlost;
                    }
                }
            cout << "w(" << i << ", " << j << ", " << k << ") = " << weights[i][j][k] << endl;
            }
        }
    }

    // auxiliary arrays
    vector<vector<vector<long double>>> S(5,vector<vector<long double>>(5,vector<long double>(2,0)));

    vector<long double> R(5,0);

    // main simulation loop
    // Initialize RNG
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    // Perform runs simulations:
    for (int run=0; run < runs; run++) {
        // perform in a timed while loop instead? prevents stochastic hanging
 
		// declare dynamical variables

        double tt = 0;
        int itt = 0; // number of writes
        long double nTotal = Nbeans; // number of cells

        vector<vector<vector<long double>>> 
            n(5,vector<vector<long double>>(5,vector<long double>(2,0)));
        // initial condition:
        n[0][0][0]=Nbeans; 
        // Sporadic cancer: all N cells initially on {0,0,0}
        long double Gamma = 0;

        vector<vector<long double>> Sij(5,vector<long double>(5,0));
        vector<long double> Si(5,0);

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                for (int k=0; k<2; k++) {
                    S[i][j][k] = n[i][j][k]*weights[i][j][k];
                    Sij[i][j] += S[i][j][k];
                }
                Si[i] += Sij[i][j];
            }
            Gamma += Si[i];
        }

        // Pre-calculate axis hopping propensities
        // Here we exploit the product structure of the space: i.e. each gene
        // mutates independently of the others.

        vector<long double> RatesI(5,0);
        vector<long double> RatesJ(5,0);
        vector<long double> RatesK(2,0);
        
        for (int i=0; i<5; i++) {
            RatesI[i] = -rates[i][i];
        }

        for (int i=0; i<5; i++) {
            RatesJ[i] = -rates2[i][i];
        }

        for (int i=0; i<2; i++) {
            RatesK[i] = -rates3[i][i];
        }

        bool Continue = 1;

        int Z = runs;
        
    
        // Main loop:
        
        while (Continue) {
            if (tt > ttmax)
                Continue = 0;
            
            long double Gamma_prev = Gamma;

            long double u = gsl_ran_flat(r, 0., Gamma);

            // convert random number into event

            int ie=0;
            int je=0;
            int ke=0;

            // determine hopping site from u
            // prob(i,j,k) \propto S[i][j][k]
            // so how long does it take to iterate through S?
            
            long double remainder = u;

            while ((ie < 5) && (remainder > Si[ie])) {
                remainder -= Si[ie];
                ie++;
            }
            if (ie > 4)
                ie = 4; // correct for long double --> double precision
            
            while ((je < 5) && (remainder > Sij[ie][je])) {
                remainder -= Sij[ie][je];
                je++;
            }
            if (je > 4)
                je = 4; // correct for long double --> double precision
            
            while ((ke < 2) && (remainder > S[ie][je][ke])) {
                remainder -= S[ie][je][ke];
                ke++;
            }
            if (ke > 1)
                ke = 1; // correct for long double --> double precision
            
            // we now have the id of the next event:
            // ie,je,ke
            // we now have to find the target of the event.

            // target selection works 
            int it = ie;
            int jt = je;
            int kt = ke;

            // determine axis of hopping
            // this is the hopping propensity:

            // prop_hop = (ratesi[ie]+ratesj[je]+ratesk[ke])*n[ie][je][ke]
            // the total prospensity is S[ie][je][ke]

            // we have a uniform random variable called "remainder" in between 0 
            // and S[ie][je][ke]. this is formed from the random variable u, 
            // so we have called the RNG once.

            remainder /= n[ie][je][ke];

            // hopping
            if (remainder < RatesI[ie] + RatesJ[je] + RatesK[ke]) {
                if (remainder < RatesI[ie]) {
                    // "it" is the hopping axis.
                    // where is it hopping to?
                    // define the RV w by w = remainder - 0
                    // bracket w by rates[ie][it]
                    while ((++it < 5) && (remainder > rates[ie][it])) {
                        remainder -= rates[ie][it];
                    }
                    if (it > 4)
                        it = 4;
                } else {
                    remainder -= RatesI[ie];
                    if (remainder < RatesJ[je]) {
                        // "jt" is the hopping axis
                        // where is it hopping to?
                        // define the RV w by w = remainder - RatesI[ie]n[ie,je,ke]
                        while ((++jt < 5) && (remainder > rates2[je][jt])) {
                            remainder -= rates2[je][jt];
                        }
                        if (jt > 4)
                            jt = 4;
                    }
                    else {
                        // kt is the hopping axis
                        kt = 1; // k is special: the only possible site to hop to is kt=1
                    }
                }

                // make leap
                
                n[ie][je][ke]--; // leaves this node...
                n[it][jt][kt]++; // arrives this node

                // update propensities

                long double we = weights[ie][je][ke];
                S[ie][je][ke] -= we;
                Sij[ie][je] -= we;
                Si[ie] -= we;
                Gamma -= we;
                
                long double wt = weights[it][jt][kt];
                S[it][jt][kt] += wt;
                Sij[it][jt] += wt;
                Si[it] += wt;
                Gamma += wt;
            } else {
                // fission
                                
                n[ie][je][ke]++;
                // update propensities
                long double we = weights[ie][je][ke];
                S[ie][je][ke] += we;
                Sij[ie][je] += we;
                Si[ie] += we;
                Gamma += we;

                nTotal += 1.0;
            }

            // NEW: if the target is an end state (331, 341 etc), then 
            // skip the rest of the simulation.

            double Deltat;
            int IfEnd = 0;
            IfEnd += (n[3][3][1] > 0);
            IfEnd += (n[3][4][1] > 0);
            IfEnd += (n[4][3][1] > 0);
            IfEnd += (n[4][4][1] > 0);
            
            if (IfEnd>0) {
                Deltat = ttmax;
                tt = ttmax + 0.5;
                Continue = 0;
            } else {
                double Tau = 1./Gamma_prev;

                // get random time step
                Deltat = gsl_ran_exponential(r, Tau);
            }

            // each time tt crosses 1.0, write out data
            while (tt>=(double)itt) {
                // Write out nTotal also
                (*tResults)[itt][make_tuple(-1,-1,-1)] = Z;

                for (int i=0; i<5; i++) {
                    for (int j=0; j<5; j++) {
                        for (int k=0; k<2; k++) {
                            if (n[i][j][k]>0) {
                                (*tResults)[itt][make_tuple(i,j,k)]+=1.0;
                            }
                        }
                    }
                }

                // how far through this run are we?
                double frac = tt/ttmax;
                double percent = (100*((double)run+frac)/(double)runs);
                int ipct = (int)percent;

                progress->store(ipct);

                itt++;
            }

            tt += Deltat;
            
            if (tt >= ttmax + 1)
                Continue = 0;
		}
    }

    // Report that we are done:
    progress->store(100);
    gsl_rng_free(r);
}

void tauleap_beans (vector<Ensemble> *tResults, int seed, 
        atomic<int> *progress, simConstants Params) {
    // Pull in constant parameters passed in from outside: 
    double ttmax = Params.ttmax;
    double dt = Params.dt;
    int runs = Params.runs;
    set<State> G = Params.genotypes;
    map<State,set<State>> Neighbours = Params.Neighbours;

    // the main object being manipulated is a mapping of the state to the number
    // of beans in that state: this is the population.
    // declare core objects, initialize inside.
    Ensemble population;

    // main simulation loop
    // Initialize RNG
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    // Initializations complete.

    // Perform runs simulations:
    for (int run=0; run<runs; run++) {
        double tt = 0;
        State origin = make_tuple(0,0,0); // NEW: PLANE NOT CUBE

        // Initialize population vector:
        reset_ensemble(population);
        if (Params.mode == 0) {
            population[origin]=(uint64_t)Nbeans; // all N beans initially on {0,0,0}
        }
        if (Params.mode == 1) {
            population[origin]=(uint64_t)shiftNbeans;
        }

        // Main simulation loop:
        // count how many times we visit the last node:
        int lastvisit = 0;
        int TrackLastMut = 0; // 0: raw prevalence, 1: last mutation to cause cancer

        while (tt < ttmax) {
            // for each possible genotype,
            for (auto p = G.begin(); p != G.end(); ++p) {
                // perform MUTATION:
                // for each neighbouring genotype,
                set<State> H = Neighbours[*p];

                for (auto q = H.begin(); q != H.end(); ++q) {
                    // By iterating over neighbours only we save some time.
                    double prob = flow(*p,*q)*dt;
                    // Also, prob > 0 due to how we calculate neighbours
                    // but watch out for large time steps etc
                    if (prob > 0.01) prob = 1-exp(-prob);

                    // move some number from p to q
                    uint64_t dn = gsl_ran_binomial(r,prob,population[*p]);

                    if (dn > 0) {
                        population[*p] -= dn;
                        population[*q] += dn;
                        // new mode: if this is the end state, then where did this bean originate?
                        // add a counter for how many times we visit the last node:
                        lastvisit++;
                        
                        if (TrackLastMut) {
                            if ((*q)==make_tuple(3,3,1)) {
                                int it = (int)(tt+0.5);

                                while (it < ttmax) {
                                    // skip the rest of the simulation:
                                    (*tResults)[it][*p] += 1;
                                    it++;
                                }
                                tt += ttmax; //skip rest of simulation
                                // as soon as we produce a cancerous cell
                            } else if ((*q)==make_tuple(3,4,1)) {
                                int it = (int)(tt+0.5);

                                while (it < ttmax) {
                                    // skip the rest of the simulation:
                                    (*tResults)[it][*p] += 1;
                                    it++;
                                }
                                tt += ttmax; //skip rest of simulation
                                // as soon as we produce a cancerous cell
                            } else if ((*q)==make_tuple(4,3,1)) {
                                int it = (int)(tt+0.5);

                                while (it < ttmax) {
                                    // skip the rest of the simulation:
                                    (*tResults)[it][*p] += 1;
                                    it++;
                                }
                                tt += ttmax; //skip rest of simulation
                                // as soon as we produce a cancerous cell
                            } else if ((*q)==make_tuple(4,4,1)) {
                                int it = (int)(tt+0.5);

                                while (it < ttmax) {
                                    // skip the rest of the simulation:
                                    (*tResults)[it][*p] += 1;
                                    it++;
                                }
                                tt += ttmax; //skip rest of simulation
                                // as soon as we produce a cancerous cell
                            }
                        }
                        
                    }
                }

                // then perform FISSION:
                if (tt < ttmax) {
                    // if the genotype has lost APC (get<0>(*p)>2)
                    if ((get<0>(*p)>2)&&(rate_APClost>0)) {
                        // then add a number of new crypts here
                        // if KRAS is not mutated:
                        double prob = exp(-rate_APClost*dt);

                        if (get<2>(*p)==1) {
                            // if it is:
                            prob = exp(-rate_BOTHlost*dt);
                        }
                        uint64_t dn = gsl_ran_negative_binomial(r,prob,population[*p]);

                        if (dn > 0) {
                            population[*p] += dn;
                        }
                    }

                    if ((get<0>(*p)<=2)&&(get<2>(*p)>0)&&(rate_RASlost>0)) {
                        // then add a number of new crypts here
                        double prob = exp(-rate_RASlost*dt);

                        uint64_t dn = gsl_ran_negative_binomial(r,prob,population[*p]);

                        if (dn > 0) {
                            population[*p] += dn;
                        }
                    }
                }
            }

            if (fmod(tt,1.0)<dt) {
                // Record results from this run:
                int it = (int)(tt+0.5);
                double nTotal = 0;

                // old mode, P(visit)
                if (!TrackLastMut) {
                    for (auto p = G.begin(); p != G.end(); ++p) {
                        if (population[*p]>0) {
                            if (Params.mode == 0) {
                                (*tResults)[it][*p]++;
                            } 
                            if (Params.mode == 1) {
                                (*tResults)[it][*p] += population[*p];
                            }

                            nTotal += population[*p]; // calculate total pop.
                        }
                    }
                }
                //(*tResults)[it][make_tuple(-1,-1,-1)] += nTotal;

                double frac = tt/ttmax;
                double percent = (100*((double)run+frac)/(double)runs);

                int ipct = (int)percent;

                progress->store(ipct); // write progress to external monitor
            }

            // increment time
            tt += dt;
        }
    }

    // Report that we are done:
    progress->store(100);

    // Free RNG
    gsl_rng_free(r);
}

void tauleap_lineages(vector<PathEnsemble> *tResults, int seed, atomic<int> *progress, simConstants Params)
{
    // Lineage-tracking mode. 
    // First, as always, pull in parameters from outside:
    double ttmax = Params.ttmax;
    double dt = Params.dt;
    int runs = Params.runs;
    set<State> G = Params.genotypes;
    map<State,set<State>> Neighbours = Params.Neighbours;

    // Instead of hopping around nodes on the graph, cells hop between nodes
    // on a tree that represents all possible paths on this graph.

    PathEnsemble population;

    // Before we can initialise values, we need to construct the "path graph", 
    // consisting of all possible
    // lineages/paths on the underlying graph of genotypes G. This "path
    // graph" L can also be seen as a graph: nodes on L are a list (vector) of nodes
    // on G. We will then need to iterate over the path graph L to set initial
    // values for all the possible paths.
    
    // L is a set of lineages (histories):
    set<History> L;
    trace_paths(L,G); // builds L from the set of genotypes G
    // The root history is the starting point:
    History root;
    root.push_back(make_tuple(0,0,0));

    // The population is distributed over lineages (histories of nodes). so it
    // should be a map that accepts a lineage, and returns an integer.
    //
    // When cells divide, the population of the same lineage goes up:
    //      n[l] += dn;
    //
    // When cells hop to a neighbouring site, the population of a neighbouring
    // lineage goes up, and that of the parent lineage goes down:
    //      n[ln] += dn;
    //      n[l] -= dn;
    // ln is a leaf node of l: ln = l.push_back(target_genotype).
    //
    // These two events should be incorporated in the algorithm as follows:
    //      Iterate over all lineages l in L.
    //      For each lineage l:
    //          Check if the last node in l, l.back(), has an advantage.
    //          If it does, divide, and increase n[l].
    //
    //          Then, iterate over all other nodes on G.
    //          For each node g on G:
    //              Check if g neighbours l.back(): is there a flow from
    //              l.back() into g?
    //              If there is, define ln = l.push_back(g), and move dn cells
    //              from l to ln.
    //
    //      Increment time by dt.
    //
    // Main simulation loop
    // initialize RNG
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    // perform simulations:
    for (int run=0; run<runs; run++) {
        double tt = 0; // set time
        int itt = 0; // used for data writing

        // re-zero population vector:
        reset_path_ensemble(population,L);
        // place N beans on the root history:
        population[root] = Nbeans;

        while (tt < ttmax) {
            // for each possible lineage:
            for (auto l = L.begin(); l != L.end(); ++l) {
                // check if the most recently visited node has an advantage.
                // the most recently visited node is (*l)back():
                Genotype Last = (*l).back();

                // perform fission if it has an advantage:
                if ((get<0>(Last)>2)&&(rate_APClost>0)) {
                    // if KRAS is not mutated:
                    double prob = exp(-rate_APClost*dt);
                    if (get<2>(Last)==1) {
                        prob = exp(-rate_BOTHlost*dt);
                    }
                    uint64_t dn = gsl_ran_negative_binomial(r,prob,population[*l]);

                    if (dn > 0) population[*l] += dn;
                }

                if ((get<0>(Last)<=2)&&(get<2>(Last)>0)&&(rate_RASlost>0)) {
                    // then add a number of new crypts here
                    double prob = exp(-rate_RASlost*dt);

                    uint64_t dn = gsl_ran_negative_binomial(r,prob,population[*l]);

                    if (dn > 0) {
                        population[*l] += dn;
                    }
                }

                // then perform mutation:
                // iterate over neighbours only to save time
                set<State> H = Neighbours[Last];

                for (auto q = H.begin(); q != H.end(); ++q) {
                    double prob = flow(Last,*q)*dt;
                    if (prob > 0.01) prob = 1-exp(-prob); //more accurate

                    uint64_t dn = gsl_ran_binomial(r,prob,population[*l]);
                    if (dn > 0) {
                        History Target = (*l); //.push_back(*q) is a method
                        Target.push_back(*q); // Target is the next history
                        // Target must already be a member of L.
                        population[*l] -= dn;
                        population[Target] += dn;

                        // If the target node (*q) is a cancerous node, skip
                        // the rest of the simulation:
                        if ((*q)==make_tuple(3,3,1)) {
                            tt += ttmax; //skip rest of simulation
                        } else if ((*q)==make_tuple(3,4,1)) {
                            tt += ttmax; //skip rest of simulation
                        } else if ((*q)==make_tuple(4,3,1)) {
                            tt += ttmax; //skip rest of simulation
                        } else if ((*q)==make_tuple(4,4,1)) {
                            tt += ttmax; //skip rest of simulation
                        }
                    }
                }
            }

            while ((tt>=(double)(itt))&&((double)itt<ttmax)) {
                // Record results from this run:
                // Prob(pop on site l.back having lineage l is >0)
                for (auto l = L.begin(); l != L.end(); ++l) {
                    if (population[*l]>0) {
                        (*tResults)[itt][*l]++;
                    }
                }

                // report progress
                double frac = tt/ttmax;
                double percent = (100*((double)run+frac)/(double)runs);
                int ipct = (int) percent;

                progress->store(ipct); // write to external monitor

                itt++;
            }

            // increment time
            tt += dt;
        }
    }

    // To quit process, report that this thread is complete:
    progress->store(100);

    // and free RNG:
    gsl_rng_free(r);
}

int main (int argc, char** argv) {
    
    // Detect number of available cores
    unsigned int nCores = thread::hardware_concurrency();
    int nThreads = nCores-2; // AGGRESSIVE RARRR
    if (nThreads < 2) nThreads = nCores;

    // Get hostname id and PID
    int hostname = 0; 

    // Intialize constant parameters:
    double ttmax = 80.01; // simulation time in years.
    // ^^ NB: we add a small extra amount to make sure we capture the data
    // point at exactly 80 years.
    double dt = 0.01;
    int runs = 100000;     // runs per thread.
    int mode;
    string outputfile = "";

    // Get command-line arguments:
    bool Default=0; // Tau Leaping
    bool Gillespie=0;
    bool Lineages=0;
    // Other modes:
    bool PrintKey=0;
    bool Polite = 0;
     
    for (int i=1; i<argc; i++) {
        // Get mode
        if ((string)argv[i] == "--default") Default = 1;
        if ((string)argv[i] == "--gillespie") Gillespie = 1;
        if ((string)argv[i] == "--polite") Polite = 1;
        if ((string)argv[i] == "--withkey") PrintKey = 1;
        if ((string)argv[i] == "--lineages") Lineages = 1;

        // Errors
        if ((string)argv[i] == "-o") {
            if (i == argc-1) {
                cout << "Error: no output file" << endl;
                return 1;
            } else {
                // set outputfile
                outputfile = (string)argv[i+1];
            }
        }

        if ((string)argv[i] == "--host") {
            if (i == argc-1) {
                cout << "Error: empty hostname" << endl;
                return 1;
            } else {
                hostname = atoi(argv[i+1]);
            }
        }

		if ((string)argv[i] == "--cores") {
			if (i == argc-1) {
				cout << "Error: enter number of threads" << endl;
				return 1;
			} else {
                // set number of threads
				nThreads = atoi(argv[i+1]);
			}
		}
    }

    if (outputfile == "") {
        cout << "Error: no output file" << endl;
        return 1;
    }

    if ((Default+Gillespie+Lineages)==0) {
        cout << "Defaulting to tau leaping without transform" <<endl;
        Default = 1;
    }

    if ((Default+Gillespie+Lineages)>1) {
        cout << "Error: cannot select more than 1 mode" << endl;
        return 1;
    } else {
        if (Default) mode = 0;
        if (Gillespie) mode = 2;
        if (Lineages) mode = 3;
    }

    if (Polite) nThreads = nCores/2;

    // Flag check done

    cout << "Mode: " << mode << endl;
    cout << hostname << endl;

    // Initialize constant objects:
    // Set of states to iterate over
    set<State> G;
    init_genotypes(G);

    // Initialize output files
    ofstream output;
    output.open(outputfile);

    // Set of possible end states ([3|4][3|4]1):
    set<State> ends;
    init_ends(ends);

    // Store all possible neighbours:
    map<State,set<State>> Neighbours = generate_neighbours(G);

    // Bundle constant parameters and objects:
    simConstants P;
    P.ttmax = ttmax;
    P.dt = dt;
    P.runs = runs;
    P.genotypes = G;
    P.Neighbours = Neighbours;
    P.mode = mode;

    // Initialize dynamical objects:
    // Initialize results vector: vector of vector of frequencies
    // For storing simulation statistics
    int bins = (int)(ttmax+0.5+1);
    vector<vector<Ensemble>> results(nThreads);

    Ensemble tmp;
    reset_ensemble(tmp); // how fast is this compared to init_ensemble?
    //init_ensemble(tmp,G); // these two are equally fast
    // Add a dummy space for storing total population
    tmp[make_tuple(-1,-1,-1)]= 0.;

    for (int i=0; i<nThreads; i++) {
        for (int it=0; it<(bins+1); it++) {
            results[i].push_back(tmp);
        }
    }

    // The lineage-tracking mode needs a new kind of results vector, which
    // maps histories to probabilities.
    vector<vector<PathEnsemble>> PathResults(nThreads);

    PathEnsemble tmp2;
    set<History> L;

    trace_paths(L,G);
    reset_path_ensemble(tmp2, L);

    for (int i=0; i<nThreads; i++) {
        for (int it=0; it<(bins+1); it++) {
            PathResults[i].push_back(tmp2); // after building an empty results
            // vector, use copies of it to fill the big aggregate data vector
        }
    }

    // Initialize vector of threads:
    vector<thread> vThreads(nThreads);

    // Initalize progress monitor vector:
    vector<atomic<int>> progress(nThreads);

    // Dynamical objects now initialized.
    // Run nThreads parallel simulations:
    if ((mode==0)) {
        for (int i=0; i < nThreads; i++) {
            cout << "Spawning thread " << i << "..." <<endl;
            int seed = hostname*nThreads+i;
            cout << "seed " << seed << endl;

            progress[i]=0;
            // Start this simulation:
            vThreads.at(i) = thread(tauleap_beans, &results[i], seed, &progress[i], P);
        }
    } else if (mode==2) {
        for (int i=0; i < nThreads; i++) {
            cout << "Spawning thread " << i << "..." <<endl;
            int seed = hostname*nThreads+i;
            cout << "seed " << seed << endl;

            progress[i]=0;
            // Start simulations:
            vThreads.at(i) = thread(gillespie_beans, &results[i], seed,
            &progress[i], P);
        }
    } else if (mode==3) {
        // lineage-tracking mode
        for (int i=0; i < nThreads; i++) {
            // should spawn threads similarly
            // careful to reference PathResults and not results
            cout << "Spawning thread " << i << "..." <<endl;
            int seed = hostname*nThreads+i;
            cout << "seed " << seed << endl;

            progress[i]=0;
            // Spawn simulations:
            vThreads.at(i) = thread(tauleap_lineages, &PathResults[i], seed,
            &progress[i], P);
        }
    }

    // Report thread progress
    bool Report = 1;  // Set flag to true to enter loop

    time_t starttime = time(0);
    time_t endtime   = starttime+200*60*60;
    time_t now = starttime;
    time_t last_write = starttime;

    while (Report) {
        // all modes MUST be compatible with this section
        cout << "\t|";
        Report = 0; // Reset flag each loop
        double Total = 0.;

        for (auto p = progress.begin(); p != progress.end(); ++p) {
            Total += (*p).load();
            // Test whether to continue reporting:
        //    cout << (*p).load() << ", ";

            Report = (Report || (*p < 100)); 
        }
        //cout << "\r" << flush;
        cout << ((int)(Total/nThreads)) << "%|\r" << flush;

        this_thread::sleep_for(chrono::milliseconds(300));

        // now write to temporary file

        now = time(0);
        if (now - last_write > 60) {
            last_write = now;

            // calculate how many runs have finished
            double Total = 0;
            for (auto p = progress.begin(); p != progress.end(); ++p) {
                // the number of runs that have finished is floor(progress*nThreads/100)
                Total += floor(nThreads*((*p).load())/100);
            }
            if (Total == 0) Total = 1; // avoid division by zero at early times
            double Ztmp = runs*nThreads;//Total;

            // write out temporary results (if it isn't a lineage tracking
            // sim):
            if (mode != 3) {
                ofstream tmpoutput;
                tmpoutput.open("timedoutput.csv");

                for (int it=0; it<bins; it++) {
                    tmpoutput << it << ", ";
                    for (auto p = G.begin(); p != G.end(); ++p) {
                        double freq = 0;

                        for (int i=0; i<nThreads; i++) {
                            freq += (double)results[i][it][*p]/Ztmp;
                        }
                        tmpoutput << freq << ", ";
                    }
                    tmpoutput << endl;
                }

                tmpoutput.close();
            }
        }
    }

    // Wait for all threads to finish
    for (int i=0; i <nThreads; i++) {
        vThreads.at(i).join();
    }

    cout << endl;
    cout << "Simulations complete." << endl;
    cout << "Writing results..." << endl;

    if (PrintKey) {
        // Write out key
        output << "tt, " ;
        if (mode != 3) {
            for (auto p = G.begin();p != G.end(); ++p) {
                output << get<0>(*p);
                output << get<1>(*p);
                output << get<2>(*p) << ", ";
            }
        } else {
            // For the lineages it is slightly more involved:
            // We have already initialised and populated a set of all
            // histories, L in this function:

            // unit test: count paths programmatically:
            int nPaths = 0;
            for (auto p = L.begin(); p != L.end(); ++p) {
                for (auto q = (*p).begin(); q != (*p).end(); ++q) {
                    //for each entry q in lineage p
                    if (q != (*p).begin()) output << "-";

                    output << "(";
                    output << get<0>(*q);
                    output << get<1>(*q);
                    output << get<2>(*q);
                    output << ")";
                }

                output << ",";
            }
            nPaths = L.size();
            output << nPaths;
        }
        output << endl;
    }

    // Calculate probabilities
    // With fitness, each set of results[i] has a variable number 
    // of beans over time

    double Z = nThreads*runs; // we don't use any other values any more

    cout << "Accuracy: +-" << 1.0/Z << endl;

    // Write out results
    for (int it=0; it<bins; it++) {
        // probability = prob that at least one particle is on an end state
        output << it << ", ";

        double zTotal = 0;
        for (int i=0; i<nThreads; i++) {
            zTotal += results[i][it][make_tuple(-1,-1,-1)];
        }
        if ((mode != 5)&&(mode != 3)) {
            for (auto p = G.begin(); p != G.end(); ++p) {
                double freq = 0;
                for (int i=0; i < nThreads; i++) {
                    // total up the results from different threads
                    // Normalize by total beans in this thread at this time
                    if (mode != 2) {
                        freq += (double)results[i][it][*p]/Z;
                    } else {
                        freq += (double)results[i][it][*p]/zTotal;
                    }
                }

                output << freq << ", ";
            }
        } else if (mode == 3) {
            // lineage tracking mode: should refer to PathResults, not
            // results

            // unit test: count paths
            int nPaths = 0;
            for (auto lineage = L.begin(); lineage != L.end(); ++lineage) {
                double freq = 0;
                for (int i=0; i < nThreads; i++) {
                    freq += (double)PathResults[i][it][*lineage]/Z;
                }

                output << freq << ", ";

                nPaths++;
            }
            output << nPaths;
        }

        output << endl;
    }

    // Close output files
    output.close();

    cout << "Finished." << endl;

    return 0;
}
