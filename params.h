/*
* Parameter initialization
*/
#include <cstdint>
using namespace std;

static const long double u        = 1.25e-8; // per year
static const long double rLOH     = 1.36e-4; // per year
static const long double nAPC     = 604;
static const long double nP53     = 73;
static const long double nKRAS    = 20;
static const uint64_t Nbeans = 1e8;
static const uint64_t shiftNbeans = 1;
// multipliers to account for fixation in crypts:
static const long double APC_mult1 = 2.1;
static const long double APC_mult2 = 2.8;
static const long double KRAS_mult = 3.6; // this should be 1.0 when RAS adv. is 0

// Transition matrices:
// A rate matrix gives the transition rates between labelled nodes on the
// underlying graph. Entry [0][1] gives the rate of transition from state 0 to
// 1 of the relevant gene, for example. The natural representation of a rate matrix 
// is as a 2D array.

// APC mutation rate matrix:
const long double rates[5][5] = {
    {(-nAPC*u-rLOH)*APC_mult1,nAPC*u*APC_mult1,rLOH*APC_mult1,0,0},
    {0,(-0.5*nAPC*u-0.5*rLOH)*APC_mult2,0,0.5*rLOH*APC_mult2,0.5*nAPC*u*APC_mult2},
    {0,0,(-0.5*nAPC*u)*APC_mult2,0.5*nAPC*u*APC_mult2,0},
    {0,0,0,0,0},
    {0,0,0,0,0}
};

// TP53 mutation rate matrix:
const long double rates2[5][5] = {
    {(-nP53*u-rLOH),nP53*u,rLOH,0,0},
    {0,(-0.5*nP53*u-0.5*rLOH),0,0.5*rLOH,0.5*nP53*u},
    {0,0,(-0.5*nP53*u),0.5*nP53*u,0},
    {0,0,0,0,0},
    {0,0,0,0,0}
};

// KRAS mutation rate matrix:
const long double rates3[2][2] = {
    {(-nKRAS*u)*KRAS_mult,nKRAS*u*KRAS_mult},
    {0,0}
};

const long double rate_LD = 0.00;
const long double rate_APClost = 0.200;
const long double rate_RASlost = 0.070;//0.000;// when this is 0 KRAS_mult = 1;
const long double rate_BOTHlost= 0.270;//0.200;//
const long double rate_P53lost = 0.00;
