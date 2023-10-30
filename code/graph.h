#ifndef __GRAPH_H__
#define __GRAPH_H__

#define HEAD_INFO

#include "file_ctrl.h"
#include "head.h"
#include "iclass.h"
#include "log.h"
#include "sfmt/SFMT.h"

using namespace std;
typedef double (*pf)(int, int);
void handle_error(const char* msg);

class Graph {
 public:
  uint32_t n, m;

  vector<vector<int>> gT;
  vector<vector<double>> probT;
  vector<double> cost;
  double budget;
  double costRate = 0.01;

  enum InfluModel { IC, LT };
  InfluModel influModel;
  void setInfuModel(InfluModel p) {
    influModel = p;
    TRACE(influModel == IC);
    TRACE(influModel == LT);
  }

  string folder;
  string graph_file;
  void readNM() {
    ifstream cin((folder + "attribute.txt").c_str());
    ASSERT(!cin == false);
    string s;
    while (cin >> s) {
      if (s.substr(0, 2) == "n=") {
        n = atoi(s.substr(2).c_str());
        continue;
      }
      if (s.substr(0, 2) == "m=") {
        m = (unsigned int)(atoll(s.substr(2).c_str()));
        continue;
      }
      ASSERT(false);
    }
    log_info("n=%u, m=%u", n, m);
    cin.close();
  }

  void readGraph() {
    size_t length;
    int fd = open((graph_file).c_str(), O_RDWR);
    if (fd == -1) handle_error("open");
    struct stat sb;
    int rc = fstat(fd, &sb);
    if (rc == -1) handle_error("fstat");

    length = sb.st_size;
    auto ptr =
        static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));  // byte by byte
    auto f = ptr;

    int gap = 2 * sizeof(int) + sizeof(double);
    // ASSERT(fin != false);
    unsigned int readCnt = 0;
    for (unsigned int i = 0; i < m; i++) {
      readCnt++;
      unsigned int a, b;
      double p;
      memcpy(&a, f, sizeof(int));
      memcpy(&b, f + sizeof(int), sizeof(int));
      memcpy(&p, f + 2 * sizeof(int), sizeof(double));
      f += gap;

      ASSERT(a < n);
      ASSERT(b < n);

      // probT[b].push_back((unsigned int)(p*UI_MAX));
      probT[b].push_back(p);
      gT[b].push_back(a);
    }

    ASSERT(readCnt == m);
    rc = munmap(ptr, length);
    close(fd);
  }

  void geneNodeCost(double costRate, bool randomCost) {
    cost.resize(n);
    if (!randomCost) {
      double maxCost = 0.0;
      for (uint32_t i = 0; i < n; i++) {
        //  cost[i] = BASIC_COST;
        cost[i] = BASIC_COST + gT[i].size() * costRate;
        maxCost = max(cost[i], maxCost);
      }

      log_info("average deg: %f", (double)m / n);
      log_info("total budget: %f", budget);
      log_info("basic cost: %f, max cost: %f, cost per edge: %f", BASIC_COST, maxCost, costRate);
    } else {
      int seed = 12345;  // fixed random seed
      srand(seed);       // set the random seed
      double maxCost = 0.0;
      double minCost = RANDOM_MAX;

      for (uint32_t i = 0; i < n; i++) {
        cost[i] = ((double)rand() / (double)(RAND_MAX)) * (RANDOM_MAX - 1.0) + 1.0;
        maxCost = max(cost[i], maxCost);
        minCost = min(cost[i], minCost);
      }

      log_info("min cost: %f, max cost: %f", minCost, maxCost);
    }
  }

  Graph(string folder, string graph_file) : folder(folder), graph_file(graph_file) {
    readNM();

    string weightFile = folder + "weights.bin";
    string adjFile = folder + "adj.bin";

    if (!seraizlied_graph_exist(weightFile) || !seraizlied_graph_exist(adjFile)) {
      log_info("There exist no weights file and adjacent matrix");
      gT = vector<vector<int>>(n, vector<int>());
      probT = vector<vector<double>>(n, vector<double>());

      readGraph();
      save_file(weightFile, probT);
      save_file(adjFile, gT);
    } else {
      log_info("starting loading graph");
      load_file(weightFile, probT);
      load_file(adjFile, gT);
      log_info("graph is successfully loaded");
    }
  }
};

double sqr(double t) { return t * t; }

void handle_error(const char* msg) {
  perror(msg);
  exit(255);
}

#endif