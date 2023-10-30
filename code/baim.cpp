#define HEAD_INFO


#include <sched.h>
#include <stdio.h>
#include <signal.h>
#include <numeric>
#include "head.h"
#include "sfmt/SFMT.h"
#include "argument.h"

#include "infgraph.h"
#include "baim.h"

// static unsigned int rr_num=0;
void OutputSeedSetToFile(vector<int> seed_set, string seedfile) {
  ofstream of;
  of.open(seedfile);
  for (int seed : seed_set) {
    of << seed << " ";
  }
  of << endl;
  of.close();
}

void run_with_parameter(InfGraph &g, const Argument &arg) {
  log_info("algorithm: %d", arg.algo);
  log_info("reuse strategy: %s", g.reuse ? "yes" : "no");
  log_info("use Subsim: %s", g.subsim ? "yes" : "no");
  log_info("lazyInvList opt: %s", g.lazyInvList ? "yes" : "no");
  log_info("use MultiThreads: %s", g.useMultiThread ? "yes" : "no");

  if (arg.algo == ADAP_CA_GREEDY) {
    BAIM::costAwareGreedy(g, arg);
  } else if (arg.algo == NONADAP_CA_GREEDY) {
    log_info("non adaptive cost-aware greedy algorithm");
    BAIM::runBIM(g, arg);
  } else if (arg.algo == MAX_INF_SINGLETON) {
    log_info("the singleton with max influence");
    BAIM::runMaxInfSingleton(g, arg);
  }

  // BAIM::BudgetInfluenceMaximize(g, arg);
  INFO(g.seedSet);
  OutputSeedSetToFile(g.seedSet, arg.seedfile);

  Timer::show(1);
}


void Run(int argn, char **argv) {
  Argument arg;

  for (int i = 0; i < argn; i++) {
    if (argv[i] == string("-help") || argv[i] == string("--help") || argn == 1) {
      cout << "./imm -dataset *** -k ***  -model IC|LT -seedfile *** -time *** -batch "
              "***"
           << endl;
      return;
    }

    if (argv[i] == string("-dataset")) arg.dataset = argv[i + 1];
    if (argv[i] == string("-alpha")) arg.alpha = atof(argv[i + 1]);
    if (argv[i] == string("-beta")) arg.beta = atof(argv[i + 1]);
    if (argv[i] == string("-model")) arg.model = argv[i + 1];
    if (argv[i] == string("-seedfile")) arg.seedfile = argv[i + 1];
    if (argv[i] == string("-time")) arg.time = atoi(argv[i + 1]);
    if (argv[i] == string("-reuse")) arg.reuse = (atoi(argv[i + 1]) == 1);
    if (argv[i] == string("-subsim")) arg.subsim = (atoi(argv[i + 1]) == 1);
    if (argv[i] == string("-budget")) arg.budget = atoi(argv[i + 1]);
    if (argv[i] == string("-algo")) arg.algo = atoi(argv[i + 1]);
    if (argv[i] == string("-mthreads")) arg.mthreads = (atoi(argv[i + 1]) == 1);
    if (argv[i] == string("-rate")) arg.rate = atof(argv[i + 1]);
    if (argv[i] == string("-randomcost")) arg.randomCost = (atoi(argv[i + 1]) == 1);
  }

  ASSERT(arg.dataset != "");
  ASSERT(arg.model == "IC" || arg.model == "LT");
  arg.dataset = arg.dataset + "/";
  string graph_file;
  if (arg.model == "IC")
    graph_file = arg.dataset + "graph_ic.inf";
  else if (arg.model == "LT")
    graph_file = arg.dataset + "graph_lt.inf";
  else
    ASSERT(false);

  InfGraph g(arg.dataset, graph_file);
  g.setConfig(arg);
  g.geneNodeCost(arg.rate, arg.randomCost);
  g.createMultiThreadsEnv();

  if (arg.model == "IC")
    g.setInfuModel(InfGraph::IC);
  else if (arg.model == "LT")
    g.setInfuModel(InfGraph::LT);
  else
    ASSERT(false);

  run_with_parameter(g, arg);
}


int main(int argn, char **argv) {

  __head_version = "v1";
  OutputInfo info(argn, argv);

  Run(argn, argv);
}
