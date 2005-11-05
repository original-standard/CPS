#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_gauge.C
//
// AlgActionGauge is represents the pure gauge contribution to the QCD
// action
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE

AlgActionGauge::AlgActionGauge(AlgMomentum &mom, ActionGaugeArg &g_arg)
  : AlgAction(mom)
{

  cname = "AlgActionGauge(Gclasstype, M*)";
  gauge_arg = &g_arg;
  gluon = g_arg.gluon;

}

AlgActionGauge::~AlgActionGauge() {

}

//!< Heat Bath for the gauge action (i.e., does nothing)
void AlgActionGauge::heatbath() {

}

//!< Calculate gauge contribution to the Hamiltonian
Float AlgActionGauge::energy() {

  char *fname = "energy()";
  Lattice &lat = 
    LatticeFactory::Create(F_CLASS_NONE, gluon);
  Float h = lat.GhamiltonNode();
  LatticeFactory::Destroy();
  return h;

}

//!< evolve method evolves the momentum due to the gauge force
void AlgActionGauge::evolve(Float dt, int steps) 
{

  char *fname = "evolve(Float,int)";
  //!< Create an appropriate lattice
  Lattice &lat = 
    LatticeFactory::Create(F_CLASS_NONE, gluon);

  for (int i=0; i<steps; i++) 
    lat.EvolveMomGforce(mom, dt);

  LatticeFactory::Destroy();

}

//!< Dummy methods
void AlgActionGauge::cost(CgStats *cg_stats_global) {

}

void AlgActionGauge::init() {

}

CPS_END_NAMESPACE