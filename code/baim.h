#pragma once

#include <chrono>
#include <ctime>
#include <ratio>

#define e exp(1)
#define c 2 * (exp(1) - 2)

using namespace std::chrono;

class Math {
 public:
  static double log2(int n) { return log(n) / log(2); }
  static double logcnk(int n, int k) {
    double ans = 0;
    for (int i = n - k + 1; i <= n; i++) {
      ans += log(i);
    }
    for (int i = 1; i <= k; i++) {
      ans -= log(i);
    }
    return ans;
  }
};

class BAIM {
 public:
  static unsigned long long rr_num;
  static unsigned int left_n;

  static void loadSubsimSeedSet(const Argument &arg, int id, vector<int> &seedSet)
  {
    // dataset/graph/
    int len = arg.dataset.length();
    int start_pos = string("dataset/").length();
    string graph_name = arg.dataset.substr(start_pos, len-start_pos - 1);
    cout << graph_name << endl;

    string filename = string("dataset/subsim/") + "seed_" + graph_name + "_subsim_k" + to_string(arg.budget) +
                      "_" + to_string(id);
    cout << filename << endl;
    ifstream cin(filename);
    if (!cin.is_open())
    {
      cout << "error in open " << filename << endl;
      exit(1);
    }

    int node;
    while (cin >> node) 
    {
      seedSet.push_back(node);
    }
    log_info("seed size: %lu\n", seedSet.size());
    return;
  }

  static void maxInfSingleton(InfGraph &g, const Argument &arg, double &upperBound) {
    const double approx = arg.beta;
    double epsilon = 1.0 - arg.beta;

    double left_n = g.NumcurNode;

    double delta_i = 0.01 * epsilon / left_n;
    double epsilon_i = (epsilon - delta_i * left_n) / (1.0 - delta_i * left_n);
    double epsilon_a = epsilon_i / (1.0 - delta_i * left_n);

    double tmp_var = (2.0 + 2.0 * epsilon_a / 3.0) * left_n / epsilon_a / epsilon_a;
    size_t numIter = ceil(log(tmp_var) / log(2)) + 1;
    double a_i = log(2.0 * numIter / delta_i);
    size_t numRBase = log(2.0 / delta_i) + log(left_n);

    for (size_t idx = 1; idx <= numIter; idx++) {
      const auto numR = numRBase << (idx - 1);

      log_debug("iteration: %lu, sample size: %lu", idx, numR);
      g.build_RRsets(numR, g.R1);  // R1
      g.build_RRsets(numR, g.R2);  // R2

      double upperCoverage = g.NumcurNode;
      g.maxInfSingleton_selectNode(g.maxInfSeedSet, g.R1, a_i, upperCoverage);
      double coverage = g.Coverage(g.maxInfSeedSet, g.R2);

      const auto lowerSelect =
          pow2(sqrt(coverage + a_i * 2.0 / 9.0) - sqrt(a_i / 2.0)) - a_i / 18.0;
      const auto currApprox = lowerSelect / upperCoverage;

      log_debug("current approx: %f", currApprox);

      if (currApprox >= approx) {
        double tmp = pow2(sqrt(upperCoverage + a_i / 2.0) + sqrt(a_i / 2.0));
        upperBound = tmp * g.n / g.R2.validRRsetNum;
        log_debug("the cost of maxInf: %f", g.cost[g.maxInfSeedSet[0]]);
        log_debug("expected non adaptive greedy influence: %f, upperBound: %f",
                 coverage * g.n / g.R2.validRRsetNum, upperBound);
        return;
      }
    }
    return;
  }

  static void BIM(InfGraph &g, const Argument &arg, double &lowerBound) {
    const double approx = arg.beta;
    double delta = 1.0 / g.NumcurNode;
    double left_n = g.NumcurNode;
    double epsilon = 1.0 - arg.beta;

    log_info("BIM approximation: %f", approx);

    double tmp_var1 = sqrt(log(6.0 / delta)) + sqrt(log(left_n) + log(6.0 / delta));
    // a multiply factor 100 to increase the numRbase
    size_t numRbase = 100.0 * tmp_var1 * tmp_var1 * 2.0 / epsilon / epsilon;

    uint32_t k_max = g.budget / BASIC_COST;
    double tmp_var2 = sqrt(log(6.0 / delta)) + sqrt(log(left_n) * (double)k_max + log(6.0 / delta));
    size_t numRMax = 2.0 * tmp_var2 * tmp_var2 * left_n / epsilon / epsilon / (double)k_max;

    int numIter = ceil(log2(numRMax / numRbase)) + 1;

    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);

    for (auto idx = 1; idx <= numIter; idx++) {
      const auto numR = numRbase << (idx - 1);

      log_info("iteration: %d, sample size: %lu", idx, numR);
      g.build_RRsets(numR, g.R1);  // R1
      g.build_RRsets(numR, g.R2);  // R2

      double upperCoverage = g.max_cover_lazy_budget(g.R1);
      double coverage = g.Coverage(g.nonAdaSeedSet, g.R2);

      const auto lowerSelect = pow2(sqrt(coverage + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
      const auto upperOPT = pow2(sqrt(upperCoverage + a2 / 2.0) + sqrt(a2 / 2.0));
      const auto currApprox = lowerSelect / upperOPT;

      log_info("current approx: %f", currApprox);

      if (currApprox >= approx) {
        log_info("lower bound: %f", lowerSelect * g.n / g.R2.validRRsetNum);
        log_info("expected non adaptive greedy influence: %f", coverage * g.n / g.R2.validRRsetNum);

        lowerBound = lowerSelect * g.n / g.R2.validRRsetNum;
        return;
      }
    }

    return;
  }

  static void ourCASE(InfGraph &g, const Argument &arg) {
    uint32_t left_n = g.NumcurNode;
    double delta = BASIC_COST / (double)left_n;
    double epsilon_a = (1.0 - delta) / arg.alpha - 1.0;

    double tmp_variable = (2. + 2. * epsilon_a / 3) * left_n / epsilon_a / epsilon_a;
    uint64_t r_max = ceil(tmp_variable * log((double)left_n / delta)) + 1;
    uint64_t sample = ceil(log((double)left_n / delta));
    double a = log(1.0 / delta);

    int currSeed;
    vector<int> batch_set;
    log_debug("r max: %lu", r_max);
    double ratio = 0.0;

    while (sample < r_max) {
      g.build_RRsets((RID_T)sample, g.R1);

      double upperRatio = left_n / BASIC_COST;
      g.CASE_selectNode(currSeed, g.R1, a, upperRatio);

      upperRatio = upperRatio * left_n / g.R1.validRRsetNum;

      g.build_RRsets((RID_T)sample, g.R2);

      batch_set.clear();
      batch_set.push_back(currSeed);
      double cover = g.Coverage(batch_set, g.R2);
      if (g.reuse) {
        int tempSeed;
        g.CASE_selectNode(tempSeed, g.R2, a, upperRatio);
        upperRatio = upperRatio * left_n / g.R2.validRRsetNum;
      }

      log_debug("sample: %lu, g.R1 validNum: %lu, g.R2 validNum: %lu", (uint64_t)sample,
                (uint64_t)g.R1.validRRsetNum, (uint64_t)g.R2.validRRsetNum);
      log_debug("cover ratio: %f", cover / g.R2.validRRsetNum);

      double lowerRatio = sqr(sqrt(cover + 2. * a / 9.) - sqrt(a / 2.)) - a / 18;
      lowerRatio = lowerRatio / g.cost[currSeed];
      lowerRatio = lowerRatio * left_n / g.R2.validRRsetNum - delta * left_n / g.cost[currSeed];

      log_debug("lower ratio: %f, upper ratio:%f", lowerRatio, upperRatio);

      ratio = lowerRatio / upperRatio;

      log_debug("seed cost: %f", g.cost[currSeed]);
      log_debug("lower: %f, upper: %f, ratio: %f", lowerRatio, upperRatio, ratio);
        // log_info("ratio: %f", ratio);
      if (ratio > arg.alpha) {

        // log_info("seed cost: %f, consumed cost: %f", g.cost[currSeed], g.totalCost);

        if (g.totalCost + g.cost[currSeed] > g.budget) {
          log_debug("total seed num: %lu", g.seedSet.size());
          g.needToStop = true;
        } else {
          g.requiredSampleNum += (g.R1.validRRsetNum + g.R2.validRRsetNum);
          g.totalCost += g.cost[currSeed];
          g.seedSet.push_back(currSeed);

          g.realization(batch_set);
          g.avgRatio += ratio;
        }

        return;
      }
      sample = std::min(g.R1.validRRsetNum, g.R2.validRRsetNum);
      sample *= 2;
    }


    g.requiredSampleNum += (g.R1.validRRsetNum + g.R2.validRRsetNum);
    g.totalCost += g.cost[currSeed];
    g.seedSet.push_back(currSeed);

    g.realization(batch_set);
    g.avgRatio += ratio;

    return;
  }

  static void costAwareGreedy(InfGraph &g, const Argument &arg) {
    double total_spread = 0;
    int64_t total_actual_new_sample = 0;
    int64_t total_actual_update_sample = 0;

    double lastUsedTime = 0;
    for (int i = 0; i < arg.time; i++) {
      g.load_possible_world(to_string(i), arg);
      g.init_hyper_graph();
      g.avgRatio = 0;

      left_n = g.n;
      {
        Timer timer(SUMMARY, "total running time");

        while (!g.needToStop) {
          log_trace("seedsize=%lu", g.seedSet.size());
          left_n = g.NumcurNode;
          ourCASE(g, arg);

        }
        log_debug("avg ratio: %f", g.avgRatio / g.seedSet.size());
      }

      total_spread += (g.n - g.NumcurNode);
      total_actual_new_sample += g.actualNewSampleNum;
      total_actual_update_sample += g.actualUpdateNum;
      rr_num += g.requiredSampleNum;

      cout << "SingleSpread " << g.n - g.NumcurNode << endl;
      cout << "SingleRuntime " << Timer::getSpecificInterval(SUMMARY) - lastUsedTime << endl;
      lastUsedTime = Timer::getSpecificInterval(SUMMARY);
      if (arg.reuse) {
        cout << "Single actual sample " << g.actualNewSampleNum << endl;
        cout << "Single update sample " << g.actualUpdateNum << endl;
      }
      cout << "Single required sample " << g.requiredSampleNum << endl;
      cout << "*********************************" << endl;
    }
    cout << "RunningTime(s) " << Timer::getSpecificInterval(SUMMARY) / arg.time << endl;
    cout << "TotalActualNewSample: " << total_actual_new_sample / arg.time << endl;
    cout << "TotalUpdateSample: " << total_actual_update_sample / arg.time << endl;
    cout << "TotalRequiredSample: " << rr_num / arg.time << endl;
    cout << "Total AdaptiveSpread: " << total_spread / arg.time << endl;
  }

  static void runBIM(InfGraph &g, const Argument &arg) {
    double total_spread = 0;

    double lastUsedTime = 0;
    for (int i = 0; i < arg.time; i++) {
      g.load_possible_world(to_string(i), arg);
      g.init_hyper_graph();

      double lowerBound = 0;
      {
        Timer timer(SUMMARY, "total running time");
        BIM(g, arg, lowerBound);
        log_info("seed size: %lu", g.nonAdaSeedSet.size());
      }
      g.realization(g.nonAdaSeedSet);

      total_spread += (g.n - g.NumcurNode);

      cout << "SingleSpread " << g.n - g.NumcurNode << endl;
      cout << "SingleRuntime " << Timer::getSpecificInterval(SUMMARY) - lastUsedTime << endl;
      lastUsedTime = Timer::getSpecificInterval(SUMMARY);
      cout << "*********************************" << endl;
    }
    cout << "RunningTime(s) " << Timer::getSpecificInterval(SUMMARY) / arg.time << endl;
    cout << "Total AdaptiveSpread: " << total_spread / arg.time << endl;
  }

  static void runMaxInfSingleton(InfGraph &g, const Argument &arg) {
    double total_spread = 0;

    double lastUsedTime = 0;
    for (int i = 0; i < arg.time; i++) {
      g.load_possible_world(to_string(i), arg);
      g.init_hyper_graph();

      double upperBound = 0;
      {
        Timer timer(SUMMARY, "total running time");
        maxInfSingleton(g, arg, upperBound);
        g.maxInfSingleton_extendSeedSet(g.R1);
      }
      g.realization(g.maxInfSeedSet);

      total_spread += (g.n - g.NumcurNode);

      cout << "SingleSpread " << g.n - g.NumcurNode << endl;
      cout << "SingleRuntime " << Timer::getSpecificInterval(SUMMARY) - lastUsedTime << endl;
      lastUsedTime = Timer::getSpecificInterval(SUMMARY);
      cout << "*********************************" << endl;
    }
    cout << "RunningTime(s) " << Timer::getSpecificInterval(SUMMARY) / arg.time << endl;
    cout << "Total AdaptiveSpread: " << total_spread / arg.time << endl;
  }

};

unsigned long long BAIM::rr_num = 0;
unsigned int BAIM::left_n = 0;