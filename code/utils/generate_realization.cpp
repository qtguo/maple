#include <stdio.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>
// for mmap:
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

// for cusomer
#include <time.h>
#include <unistd.h>  //close open
#include <unistd.h>

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "../sfmt/SFMT.h"

// typedef unsigned long long int;

#define ASSERT(v)                                                      \
  {                                                                    \
    if (!(v)) {                                                        \
      cerr << "ASSERT FAIL @ " << __FILE__ << ":" << __LINE__ << endl; \
      exit(1);                                                         \
    }                                                                  \
  }

using namespace std;

void handle_error(const char* msg) {
  perror(msg);
  exit(255);
}

void readNM(string& folder, int& n, int64_t& m) {
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

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "Usage: ./gene_real dataset_folder number" << endl;
    exit(0);
  }

  int n;
  int64_t m;

  string folder = argv[1];
  folder = folder + "/";
  int num = atoi(argv[2]);

  sfmt_t sfmtSeed;
  // srand(95082);
  srand(time(NULL));
  // sfmt_init_gen_rand(&sfmtSeed, 95082);
  sfmt_init_gen_rand(&sfmtSeed, rand());

  readNM(folder, n, m);
  vector<vector<int>> in_edge(n, vector<int>());

  size_t length;
  string graph_file = folder + "graph_ic.inf";
  int fd = open((graph_file).c_str(), O_RDWR);
  if (fd == -1) handle_error("open");
  struct stat sb;
  int rc = fstat(fd, &sb);
  if (rc == -1) handle_error("fstat");

  length = sb.st_size;
  auto ptr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));  // byte by
                                                                                      // byte
  auto f = ptr;

  int gap = 2 * sizeof(int) + sizeof(double);
  // ASSERT(fin != false);

  cout << n << " " << m << endl;
  for (int64_t i = 0; i < m; i++) {
    int a, b;
    double p;
    // int c = fscanf(fin, "%d%d%lf", &a, &b, &p);
    memcpy(&a, f, sizeof(int));
    memcpy(&b, f + sizeof(int), sizeof(int));
    memcpy(&p, f + 2 * sizeof(int), sizeof(double));
    f += gap;

    ASSERT(a < n);
    ASSERT(b < n);

    in_edge[b].push_back(a);
  }

  rc = munmap(ptr, length);
  close(fd);

  // IC
  for (int k = 0; k < num; k++) {
    string index = to_string(k);
    string outfile = folder + "realization_" + index;
    ofstream output(outfile);

    for (int i = 0; i < n; i++) {
      for (size_t j = 0; j < in_edge[i].size(); j++) {
        // if ((double)rand() / RAND_MAX < 1.0 / in_edge[i].size())output << in_edge[i][j] << " " <<
        // i << endl;
        if (sfmt_genrand_real1(&sfmtSeed) < 1.0 / in_edge[i].size())
          output << in_edge[i][j] << " " << i << endl;
        // if ((double)rand() / RAND_MAX < probT[i][j])output << in_edge[i][j] << " " << i << endl;
        // if ((double)rand() / RAND_MAX < probT[i][j])output << in_edge[i][j] << " " << i << endl;
      }
    }

    output.close();
  }

  return 0;
}
