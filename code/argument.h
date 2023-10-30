#pragma once
#include <stdio.h>
#include <string>

using namespace std;
class Argument {
 public:
  string dataset;
  double alpha = 0.5;
  double beta = 0.8;
  double rate = 0.01;

  string model;
  string seedfile;
  int time;  // this is the number of times needed to measure the result.
  bool reuse = true;
  bool subsim = true;
  uint32_t budget = 100;
  int algo = 1;
  bool mthreads = false;
  bool randomCost = false;
};