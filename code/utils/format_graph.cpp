#include <stdlib.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define ASSERT(v)                                                      \
  {                                                                    \
    if (!(v)) {                                                        \
      cerr << "ASSERT FAIL @ " << __FILE__ << ":" << __LINE__ << endl; \
      exit(1);                                                         \
    }                                                                  \
  }

void readNM(string &folder, int &n, int64_t &m) {
  ifstream cin(folder + "attribute.txt");
  ASSERT(!cin == false);
  string s;
  while (cin >> s) {
    if (s.substr(0, 2) == "n=") {
      n = atoi(s.substr(2).c_str());
      continue;
    }
    if (s.substr(0, 2) == "m=") {
      m = atoll(s.substr(2).c_str());
      continue;
    }
    ASSERT(false);
  }
  cin.close();
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cout << "./format_graph graph_folder" << endl;
    return 0;
  }

  string folder = argv[1];
  folder += "/";
  string graph_file = folder + "graph.txt";
  ifstream input(graph_file);

  if (!input.is_open()) {
    cout << "error in opening file " << graph_file << endl;
    return 0;
  }

  int n;
  int64_t m;
  readNM(folder, n, m);
  cout << "n = " << n << ", m = " << m << endl;

  vector<vector<int>> graph;
  graph.resize(n);
  int s, t;

  long long verify_m = 0;
  vector<int> indeg(n, 0);

  while (input >> s >> t) {
    assert(s < n);
    assert(t < n);
    graph[s].push_back(t);
    indeg[t]++;
    verify_m++;
  }
  assert(verify_m == m);

  vector<double> prob(n);
  for (int node = 0; node < n; node++) {
    prob[node] = 1.0 / indeg[node];
  }

  string output = folder + "graph_ic.inf";
  FILE *ofile = fopen(output.c_str(), "w");
  if (!ofile) {
    cout << "error in opening file " << output << endl;
    return 0;
  }

  for (long long i = 0; i < (long long)graph.size(); i++) {
    auto nbrs = graph[i];
    for (int j = 0; j < (int)nbrs.size(); j++) {
      int target = nbrs[j];
      fwrite(&i, sizeof(int), 1, ofile);
      fwrite(&target, sizeof(int), 1, ofile);
      double p = prob[target];
      fwrite(&p, sizeof(double), 1, ofile);
    }
  }
  fclose(ofile);
  return 0;
}
