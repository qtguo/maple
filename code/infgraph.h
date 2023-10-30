#pragma once

#include <chrono>
#include <ctime>
#include <future>
#include <numeric>
#include <queue> //priority_queue
#include <set>
#include <unordered_set>
#include <utility> // pair

#include "argument.h"
#include "graph.h"
#include "iclass.h"
using namespace std::chrono;

class InfGraph : public Graph
{
private:
  vector<bool> visit;
  vector<int> visit_mark;
  vector<int> visit_pos;
  vector<int> curNode;
  vector<int> curNodeIdx;

public:
  unsigned int NumcurNode;
  bool reuse;
  bool subsim;
  bool lazyInvList = true;
  int64_t requiredSampleNum;
  int64_t actualNewSampleNum;
  int64_t actualUpdateNum;
  double avgRatio;

  bool useMultiThread = true;
  vector<vector<bool>> visitFlagThds;
  vector<vector<vector<pair<int, RID_T>>>> RRsetsBuff;

  struct rrCollection
  {
    bool finishUpdate;
    RID_T validRRsetNum;
    // unordered_set<int> hittedRRsets;
    vector<bool> hittedFlagVec;
    vector<RID_T> hittedIdVec;

    vector<vector<int>> hedges;
    vector<vector<int>> ePos;
    vector<vector<RID_T>> invertedList;
    vector<int> deduction;

    double newRRTime;
    double updateTime;
    double buildInvListTime;
  };

  rrCollection R1;
  rrCollection R2;

  // without reuse
  rrCollection R3;

  vector<bool> active_set;
  vector<vector<int>> PO;

  sfmt_t sfmtSeed;
  vector<sfmt_t> sfmtSeedVec;
  vector<int> seedSet;

  vector<int> nonAdaSeedSet;
  vector<int> maxInfSeedSet;

  double totalCost;
  bool needToStop;

  InfGraph(string folder, string graph_file) : Graph(folder, graph_file)
  {
    srand(time(NULL));
    sfmt_init_gen_rand(&sfmtSeed, rand());

    visit = vector<bool>(n);
    visit_mark = vector<int>(n);
    visit_pos = vector<int>(n);

    curNode = vector<int>(n);
    curNodeIdx = vector<int>(n);
    iota(curNode.begin(), curNode.end(), 0);
    iota(curNodeIdx.begin(), curNodeIdx.end(), 0);
    NumcurNode = n;

    active_set = vector<bool>(n, false);

    R1.invertedList.resize(n);
    R1.deduction.resize(n);
    R2.invertedList.resize(n);
    R2.deduction.resize(n);

    R3.invertedList.resize(n);
    R3.deduction.resize(n);

    PO.resize(n, vector<int>());
    actualNewSampleNum = 0;
    actualUpdateNum = 0;
    requiredSampleNum = 0;

  }

  void createMultiThreadsEnv()
  {
    if (useMultiThread) {
      sfmtSeedVec.resize(THD_NUM);
      for (int tid = 0; tid < THD_NUM; tid++) {
        visitFlagThds.emplace_back(vector<bool>(n, false));
        vector<vector<pair<int, RID_T>>> tmp(THD_NUM);
        RRsetsBuff.emplace_back(std::move(tmp));
        sfmt_init_gen_rand(&sfmtSeedVec[tid], sfmt_genrand_uint32(&sfmtSeed));
      }
      log_info("there are %d threads", (int)visitFlagThds.size());
    }
  }

  void setConfig(Argument &arg)
  {
    reuse = arg.reuse;
    subsim = arg.subsim;
    budget = arg.budget;
    useMultiThread = arg.mthreads;
  }

  void initRRcollection(rrCollection &R)
  {
    R.finishUpdate = true;
    R.validRRsetNum = 0;

    clearRRsets(R);

    if (reuse)
    {
      fill(R.deduction.begin(), R.deduction.end(), 0);
      R.hittedFlagVec.clear();
      R.hittedIdVec.clear();
    }
  }

  void init_hyper_graph()
  {

    initRRcollection(R1);
    initRRcollection(R2);
    initRRcollection(R3);
    actualNewSampleNum = 0;
    actualUpdateNum = 0;
    requiredSampleNum = 0;

    seedSet.clear();
    fill(active_set.begin(), active_set.end(), false);

    iota(curNode.begin(), curNode.end(), 0);
    iota(curNodeIdx.begin(), curNodeIdx.end(), 0);
    NumcurNode = n;

    totalCost = 0;
    needToStop = false;
  }

  // initialize PO each time it loads a possible world.
  void load_possible_world(string index, const Argument &arg)
  {
    // initialization
    for (unsigned int i = 0; i < PO.size(); ++i)
      PO[i].clear();

    string file_name = arg.dataset + "realization_" + index;
    cout << file_name << endl;

    if (influModel == LT)
      file_name += "_lt";

    size_t length;
    auto ptr = map_file(file_name.c_str(), length);
    auto f = ptr;
    auto tep = f;
    char buffer[256];

    auto l = f + length;

    unsigned int t1, t2;

    while (f && f != l)
    {
      if ((f = static_cast<char *>(memchr(f, '\n', l - f))))
      {
        memset(buffer, 0, sizeof(buffer));
        memcpy(buffer, tep, f - tep);
        f++;
        tep = f;
        stringstream s;
        s << buffer;
        s >> t1 >> t2;
        ASSERT(t1 < n);
        ASSERT(t2 < n);

        PO[t1].push_back(t2);
      }
    }
    munmap(ptr, length);
  }

  char *map_file(const char *fname, size_t &length)
  {
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
      handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
      handle_error("fstat");

    length = sb.st_size;

    char *addr = static_cast<char *>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
      handle_error("mmap");

    close(fd);
    return addr;
  }

  void singleThreadSampleRRSets(int tid, RID_T startId, RID_T endId, rrCollection &R)
  {
    auto &visit_flag = visitFlagThds[tid];
    auto &randomSeed = sfmtSeedVec[tid];

    const int seg = n / THD_NUM + 1;
    auto &threadBuff = RRsetsBuff[tid];

    for (RID_T id = startId; id < endId; id++)
    {
      auto &hedge = R.hedges[id];

      const auto Idx = sfmt_genrand_uint32(&randomSeed) % NumcurNode;
      auto uStart = curNode[Idx];

      unsigned int n_visit_mark = 0, curIdx = 0;
      hedge.push_back(uStart);
      n_visit_mark++;
      visit_flag[uStart] = true;
      int segId = uStart / seg;
      threadBuff[segId].emplace_back(make_pair(uStart, id));

      if (influModel == IC)
      {
        while (curIdx < n_visit_mark)
        {
          int i = hedge[curIdx++];

          if (gT[i].size() == 0)
          {
            continue;
          }
          int startPos = 0;
          int endPos = (int)gT[i].size();

          if (endPos == 1)
          {
            int v = gT[i][0];

            if (active_set[v] || visit_flag[v])
            {
              continue;
            }

            visit_flag[v] = true;
            hedge.push_back(v);
            n_visit_mark++;

            int segId = v / seg;
            threadBuff[segId].emplace_back(make_pair(v, id));
          }
          else
          {
            double p = probT[i][0];
            const double log_prob = std::log(1 - p);
            while (startPos < endPos)
            {
              double prob = sfmt_genrand_real1(&randomSeed);
              startPos += int(std::log(prob) / log_prob);

              if (startPos < 0 || startPos >= endPos)
              {
                break;
              }

              const auto v = gT[i][startPos];
              startPos++;

              if (active_set[v] || visit_flag[v])
              {
                continue;
              }

              visit_flag[v] = true;
              hedge.push_back(v);
              n_visit_mark++;

              int segId = v / seg;
              threadBuff[segId].emplace_back(make_pair(v, id));
            }
          }
        }
      }
      else if (influModel == LT)
      {
        while (curIdx < n_visit_mark)
        {
          int expand = hedge[curIdx++];

          if (gT[expand].size() == 0)
            continue;

          int index = sfmt_genrand_uint32(&randomSeed) % gT[expand].size();
          int v = gT[expand][index];
          if (active_set[v] || visit_flag[v])
            continue;

          visit_flag[v] = true;
          hedge.push_back(v);
          n_visit_mark++;
          int segId = v / seg;
          threadBuff[segId].emplace_back(make_pair(v, id));
        }
      }
      else
        ASSERT(false);

      for (unsigned int i = 0; i < n_visit_mark; i++)
        visit_flag[hedge[i]] = false;
    }
  }

  void singleThreadInvList(int block, rrCollection &R)
  {
    for (int tid = 0; tid < THD_NUM; tid++)
    {
      auto &buff = RRsetsBuff[tid];
      for (auto &elem : buff[block])
      {
        R.invertedList[elem.first].push_back(elem.second);
      }
      buff[block].clear();
    }
  }

  void multiThreadSampleRRsets(RID_T num, rrCollection &R)
  {
    if (num < 1024)
    {
      {
        Timer timer(VANILLA_RRSETS, "vanilla");
        return subsimSampleRRsets(num, R);
      }
    }

    RID_T baseId = R.hedges.size();
    RID_T endId = baseId + num;
    {
      Timer timer(MULTITHD_PREP, "threadPreprocess");
      for (RID_T id = baseId; id < endId; id++)
      {
        R.validRRsetNum++;
        R.hedges.emplace_back(vector<int>());
        // R.ePos.emplace_back(vector<int>());
      }
    }

    {
      Timer timer(MULTITHD_SAMPLE, "multiple sample");
      vector<future<void>> futures(THD_NUM);
      RID_T avg_thread_num = num / THD_NUM + 1;
      for (int tid = 0; tid < THD_NUM; tid++)
      {
        RID_T startThd = tid * avg_thread_num + baseId;
        RID_T endThd = startThd + avg_thread_num;
        endThd = min(endId, endThd);
        futures[tid] = std::async(std::launch::async, &InfGraph::singleThreadSampleRRSets, this,
                                  tid, startThd, endThd, std::ref(R));
      }

      std::for_each(futures.begin(), futures.end(), std::mem_fn(&std::future<void>::wait));
    }

    {
      Timer timer(MULTITHD_INVLIST, "build inverted list");
      vector<future<void>> futures(THD_NUM);
      for (int block = 0; block < THD_NUM; block++)
      {
        futures[block] = std::async(std::launch::async, &InfGraph::singleThreadInvList, this,
                                    block,
                                    std::ref(R));
      }

      std::for_each(futures.begin(), futures.end(), std::mem_fn(&std::future<void>::wait));
    }
  }

  void vanillaSampleRRsets(RID_T num, rrCollection &R)
  {
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    RID_T id = R.hedges.size();
    RID_T endID = id + num;

    for (; id < endID; id++)
    {
      R.validRRsetNum++;
      const auto Idx = sfmt_genrand_uint32(&sfmtSeed) % NumcurNode;
      auto uStart = curNode[Idx];

      unsigned int n_visit_mark = 0, curIdx = 0;
      visit_mark[n_visit_mark++] = uStart;
      visit[uStart] = true;
      R.invertedList[uStart].push_back(id);

      if (influModel == IC)
      {
        while (curIdx < n_visit_mark)
        {
          int i = visit_mark[curIdx++];
          for (int j = 0; j < (int)gT[i].size(); j++)
          {
            int v = gT[i][j];
            if (active_set[v] || visit[v])
              continue;
            double randDouble = sfmt_genrand_real1(&sfmtSeed);
            if (randDouble > probT[i][j])
              continue;

            visit[v] = true;
            visit_mark[n_visit_mark++] = v;
            R.invertedList[v].push_back(id);
            //   hyperG[v].insert(hyperId);
          }
          //   pos.push_back(n_visit_mark - 1);
        }
      }
      else if (influModel == LT)
      {
        while (curIdx < n_visit_mark)
        {
          int expand = visit_mark[curIdx++];

          if (gT[expand].size() == 0)
            continue;

          int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
          int v = gT[expand][index];
          if (active_set[v] || visit[v])
            continue;

          visit[v] = true;
          visit_mark[n_visit_mark++] = v;
          //   hyperG[v].push_back(hyperId);
          // hyperG[v].insert(hyperId);
        }
      }
      else
        ASSERT(false);

      for (unsigned int i = 0; i < n_visit_mark; i++)
        visit[visit_mark[i]] = false;

      R.hedges.emplace_back(
          std::move(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark)));
    }
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
    R.newRRTime += interval.count();
  }



  void subsimSampleRRsets(RID_T num, rrCollection &R)
  {
     if (reuse){
      subsimSampleRRsetsWithReuse(num, R);
     }
     else
     {
      subsimSampleRRsetsWithoutReuse(num, R);
     }
  }

    void subsimSampleRRsetsWithoutReuse(RID_T num, rrCollection &R)
  {
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    RID_T id = R.hedges.size();
    RID_T endID = id + num;

    for (; id < endID; id++)
    {
      R.validRRsetNum++;
      const auto Idx = sfmt_genrand_uint32(&sfmtSeed) % NumcurNode;
      auto uStart = curNode[Idx];
      int n_visit_mark = 0, curIdx = 0;

      visit_mark[n_visit_mark++] = uStart;
      visit[uStart] = true;
      R.invertedList[uStart].push_back(id);

      if (influModel == IC)
      {
        while (curIdx < n_visit_mark)
        {
          int i = visit_mark[curIdx++];
          if (gT[i].size() == 0)
          {
            continue;
          }
          int startPos = 0;
          int endPos = (int)gT[i].size();

          if (endPos == 1)
          {
            int v = gT[i][0];

            if (active_set[v] || visit[v])
            {
              continue;
            }

            visit[v] = true;
            visit_mark[n_visit_mark++] = v;
            R.invertedList[v].push_back(id);
          }
          else
          {
            double p = probT[i][0];
            const double log_prob = std::log(1 - p);
            while (startPos < endPos)
            {
              double prob = sfmt_genrand_real1(&sfmtSeed);
              startPos += int(std::log(prob) / log_prob);

              if (startPos < 0 || startPos >= endPos)
              {
                break;
              }

              const auto nbrId = gT[i][startPos];
              startPos++;

              if (active_set[nbrId] || visit[nbrId])
              {
                continue;
              }

              visit_mark[n_visit_mark++] = nbrId;
              visit[nbrId] = true;
              R.invertedList[nbrId].push_back(id);
            }
          }
        }
      }
      else if (influModel == LT)
      {
        while (curIdx < n_visit_mark)
        {
          int expand = visit_mark[curIdx++];

          if (gT[expand].size() == 0)
            continue;

          int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
          int v = gT[expand][index];
          if (active_set[v] || visit[v])
            continue;

          visit[v] = true;
          visit_mark[n_visit_mark++] = v;
          R.invertedList[v].push_back(id);
          // hyperG[v].insert(hyperId);
        }
      }
      else
        ASSERT(false);

      for (int i = 0; i < n_visit_mark; i++)
        visit[visit_mark[i]] = false;

      R.hedges.emplace_back(
          std::move(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark)));
    }
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
    R.newRRTime += (double)interval.count();
  }
  

void subsimSampleRRsetsWithReuse(RID_T num, rrCollection &R)
  {
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    RID_T id = R.hedges.size();
    RID_T endID = id + num;

    for (; id < endID; id++)
    {
      R.validRRsetNum++;
      const auto Idx = sfmt_genrand_uint32(&sfmtSeed) % NumcurNode;
      auto uStart = curNode[Idx];
      int n_visit_mark = 0, curIdx = 0;
      visit_pos.clear();

      visit_mark[n_visit_mark++] = uStart;
      visit_pos.push_back(-1);
      visit[uStart] = true;
      R.invertedList[uStart].push_back(id);

      if (influModel == IC)
      {

        while (curIdx < n_visit_mark)
        {
          int frontier = curIdx;
          curIdx++;
          int i = visit_mark[frontier];

          if (gT[i].size() == 0)
          {
            continue;
          }
          int startPos = 0;
          int endPos = (int)gT[i].size();

          if (endPos == 1)
          {
            int v = gT[i][0];

            if (active_set[v] || visit[v])
            {
              continue;
            }

            visit[v] = true;
            visit_mark[n_visit_mark++] = v;
            R.invertedList[v].push_back(id);
            visit_pos.push_back(frontier);
          }
          else
          {
            double p = probT[i][0];
            const double log_prob = std::log(1 - p);
            while (startPos < endPos)
            {
              double prob = sfmt_genrand_real1(&sfmtSeed);
              startPos += int(std::log(prob) / log_prob);

              if (startPos < 0 || startPos >= endPos)
              {
                break;
              }

              const auto nbrId = gT[i][startPos];
              startPos++;

              if (active_set[nbrId] || visit[nbrId])
              {
                continue;
              }

              visit_mark[n_visit_mark++] = nbrId;
              visit[nbrId] = true;
              R.invertedList[nbrId].push_back(id);
              visit_pos.push_back(frontier);
            }
          }
        }
      }
      else if (influModel == LT)
      {
        while (curIdx < n_visit_mark)
        {
          int expand = visit_mark[curIdx++];

          if (gT[expand].size() == 0)
            continue;

          int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
          int v = gT[expand][index];
          if (active_set[v] || visit[v])
            continue;

          visit[v] = true;
          visit_mark[n_visit_mark++] = v;
          R.invertedList[v].push_back(id);
          // hyperG[v].insert(hyperId);
        }
      }
      else
        ASSERT(false);

      for (int i = 0; i < n_visit_mark; i++)
        visit[visit_mark[i]] = false;

      R.hedges.emplace_back(
          std::move(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark)));
      R.ePos.emplace_back(
          std::move(vector<int>(visit_pos.begin(), visit_pos.begin() + n_visit_mark)));
    }
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
    R.newRRTime += (double)interval.count();
  }


  void updateRRsets(rrCollection &R)
  {
    if (R.finishUpdate)
    {
      return;
    }

    RID_T prevSize = (RID_T)R.hedges.size();
    // log_info("number of hitted: %lu, valid: %lu, ratio: %f", R.hittedIdVec.size(), R.validRRsetNum, (double)R.hittedIdVec.size()/R.validRRsetNum);
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    for (auto &rid : R.hittedIdVec)
    {
      vector<int> rrset;
      vector<int> pos;

      R.hedges[rid].swap(rrset);
      R.ePos[rid].swap(pos);

      assert(rrset.size() != 0);

      if (lazyInvList)
      {
        // the RR set id is kept in the inverted list if lazyInvList = true
        for (auto node : rrset)
        {
          R.deduction[node]++;
        }
      }

      if (!active_set[rrset[0]])
      {
        actualUpdateNum++;
        correctSingleRRsetSubsim(rrset, pos);
        assert(rrset.size() == pos.size());
        R.hedges.emplace_back(std::move(rrset));
        R.ePos.emplace_back(std::move(pos));
      }
      else
      {
        actualNewSampleNum++;
        // the target node is activated
        R.validRRsetNum--;
        subsimSampleRRsets(1, R);
      }
    }
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
    R.updateTime += (double)interval.count();

    if (lazyInvList)
    {

      // if lazyInvList = true, we only append the newly added RR set Id to the inverted list
      for (RID_T rid = prevSize; rid < (RID_T)R.hedges.size(); rid++)
      {
        for (auto &node : R.hedges[rid])
        {
          R.invertedList[node].push_back(rid);
        }
      }
    }
    else
    {
      // if lazyInvList = false, we build the inverted list from scratch
      for (auto &list : R.invertedList)
      {
        list.clear();
      }

      RID_T rid = 0;
      for (; rid < (RID_T)R.hedges.size(); rid++)
      {

        if (R.hedges[rid].size() == 0)
          continue;

        for (auto &node : R.hedges[rid])
        {
          R.invertedList[node].push_back(rid);
        }
      }
    }

    R.finishUpdate = true;

    high_resolution_clock::time_point endTime2 = high_resolution_clock::now();
    duration<double> interval2 = duration_cast<duration<double>>(endTime2 - endTime);
    R.buildInvListTime += (double)interval2.count();
  }

  void build_RRsets(RID_T totalNum, rrCollection &R)
  {
    if (reuse)
    {
      updateRRsets(R);
      auto counter = R.validRRsetNum;
      if (totalNum < counter)
      {
        return;
      }

      subsimSampleRRsets(totalNum - counter, R);
      actualNewSampleNum += totalNum - counter;

      assert(R.validRRsetNum == totalNum);
    }
    else
    {
      int incr = totalNum - R.validRRsetNum;
      assert(incr > 0);
      if (subsim && !useMultiThread)
      {
        subsimSampleRRsets(incr, R);
      }
      else
      {
        if (useMultiThread)
        {
          multiThreadSampleRRsets(incr, R);
        }
        else
        {
          vanillaSampleRRsets(incr, R);
        }
      }
    }
  }


void correctSingleRRsetSubsim(vector<int> &orig, vector<int> &pos)
  {
    int index = 0;

    for (; index < (int)orig.size(); index++)
    {
      auto node = orig[index];
      if (active_set[node])
      {
        break;
      }
    }
    assert(index < (int)orig.size());
    assert(index != 0);

    int starting_node = orig[index];

    orig.resize(index);

    unsigned int n_visit_mark = 0, curIdx = 0;

    {
      //   int pidx = 0;
      for (int idx = 0; idx < index; idx++)
      {
        auto v = orig[idx];
        visit[v] = true;
        visit_mark[n_visit_mark++] = v;
      }

      curIdx = pos[index];
      pos.resize(index);


      int v = visit_mark[curIdx];
      bool found = false;
      auto &nbrs = gT[v];
      for (int idx = 0; idx < (int)nbrs.size(); idx++)
      {
        int nbr = nbrs[idx];
        if (starting_node == nbr)
        {
          found = true;
        }

        if (!found)
          continue;

        if (active_set[nbr] || visit[nbr])
        {
          continue;
        }
        double randDouble = sfmt_genrand_real1(&sfmtSeed);
        if (randDouble > probT[v][idx])
          continue;

        visit[nbr] = true;
        visit_mark[n_visit_mark++] = nbr;
        pos.push_back(curIdx);
      }

      curIdx++;
    }

    while (curIdx < n_visit_mark)
    {
      int frontier = curIdx;
      curIdx++;
      int i = visit_mark[frontier];
      if (gT[i].size() == 0)
      {
        continue;
      }
      int startPos = 0;
      int endPos = (int)gT[i].size();

      if (endPos == 1)
      {
        int v = gT[i][0];

        if (active_set[v] || visit[v])
        {
          continue;
        }

        visit[v] = true;
        visit_mark[n_visit_mark++] = v;
        pos.push_back(frontier);
      }
      else
      {
        double p = probT[i][0];
        const double log_prob = std::log(1 - p);
        while (startPos < endPos)
        {
          double prob = sfmt_genrand_real1(&sfmtSeed);
          startPos += int(std::log(prob) / log_prob);

          if (startPos < 0 || startPos >= endPos)
          {
            break;
          }

          const auto nbrId = gT[i][startPos];
          startPos++;

          if (active_set[nbrId] || visit[nbrId])
          {
            continue;
          }

          visit_mark[n_visit_mark++] = nbrId;
          visit[nbrId] = true;
          pos.push_back(frontier);
        }
      }
    }

    vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark).swap(orig);
    for (unsigned int i = 0; i < n_visit_mark; i++)
      visit[visit_mark[i]] = false;
  }


 void correctSingleRRset(vector<int> &orig, vector<int> &pos)
  {
    int index = 0;

    // vector<int> bak1(orig.begin(), orig.end());
    // vector<int> bak2(pos.begin(), pos.end());
    for (; index < (int)orig.size(); index++)
    {
      auto node = orig[index];
      if (!active_set[node])
      {
        continue;
      }
      else
      {
        break;
      }
    }
    assert(index < (int)orig.size());

    assert(index != 0);
    int starting_node = orig[index];
    orig.resize(index);

    unsigned int n_visit_mark = 0, curIdx = 0;
    {
      //   int pidx = 0;

      const int loc = orig.size();
      for (int idx = 0; idx < loc; idx++)
      {
        auto v = orig[idx];
        visit[v] = true;
        visit_mark[n_visit_mark++] = v;
      }

      //找到第一个点，使得ending point包含了污染点
      for (int idx = 0; idx < loc; idx++)
      {
        if (pos[idx] >= (loc - 1))
        {
          curIdx = idx;
          break;
        }
      }

      assert(pos[curIdx] >= (loc - 1));
      pos.resize(curIdx);

      int v = visit_mark[curIdx];
      bool found = false;
      auto &nbrs = gT[v];
      for (int idx = 0; idx < (int)nbrs.size(); idx++)
      {
        int nbr = nbrs[idx];
        if (starting_node == nbr)
        {
          found = true;
        }

        if (!found)
          continue;

        if (active_set[nbr] || visit[nbr])
        {
          continue;
        }
        double randDouble = sfmt_genrand_real1(&sfmtSeed);
        if (randDouble > probT[v][idx])
          continue;

        visit[nbr] = true;
        visit_mark[n_visit_mark++] = nbr;
        // hyperG[v].push_back(hyperId);
      }
      pos.push_back(n_visit_mark - 1);
      curIdx++;
    }

    while (curIdx < n_visit_mark)
    {
      int i = visit_mark[curIdx++];
      for (int j = 0; j < (int)gT[i].size(); j++)
      {
        int v = gT[i][j];
        if (active_set[v] || visit[v])
          continue;
        double randDouble = sfmt_genrand_real1(&sfmtSeed);
        if (randDouble > probT[i][j])
          continue;

        visit[v] = true;
        visit_mark[n_visit_mark++] = v;
        // hyperG[v].push_back(hyperId);
      }
      pos.push_back(n_visit_mark - 1);
    }

    vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark).swap(orig);
    for (unsigned int i = 0; i < n_visit_mark; i++)
      visit[visit_mark[i]] = false;
  }

  void getHittedRRsets(vector<int> &realization, rrCollection &R)
  {

    R.hittedFlagVec.resize(R.hedges.size());

    R.hittedIdVec.clear();

    auto &hittedFlagVec = R.hittedFlagVec;
    auto &IdVec = R.hittedIdVec;
    for (auto infNode : realization)
    {
      for (auto &idx : R.invertedList[infNode])
      {
        if (hittedFlagVec[idx])
          continue;
        hittedFlagVec[idx] = true;
        IdVec.push_back(idx);
      }
    }
    R.finishUpdate = false;
  }

  void clearRRsets(rrCollection &R)
  {
    for (auto &list : R.invertedList)
      list.clear();
    for (auto &rrset : R.hedges)
      vector<int>().swap(rrset);
    R.hedges.clear();

    for (auto &pos : R.ePos)
      vector<int>().swap(pos);
    R.ePos.clear();

    R.validRRsetNum = 0;
  }

  void realization(vector<int> &batch_set)
  {
    // printf("hello\n");
    unsigned int n_visit_mark = 0, curIdx = 0;
    vector<int> new_observation(batch_set.begin(), batch_set.end());
    for (auto seed : batch_set)
    {
      active_set[seed] = true;
      visit_mark[n_visit_mark++] = seed;
      --NumcurNode;
      curNodeIdx[curNode[NumcurNode]] = curNodeIdx[seed];
      curNode[curNodeIdx[seed]] = curNode[NumcurNode];
    }

    while (curIdx < n_visit_mark)
    {
      int expand = visit_mark[curIdx++];
      for (auto v : PO[expand])
      {
        if (active_set[v])
          continue;
        visit_mark[n_visit_mark++] = v;
        active_set[v] = true;
        --NumcurNode;
        curNodeIdx[curNode[NumcurNode]] = curNodeIdx[v];
        curNode[curNodeIdx[v]] = curNode[NumcurNode];

        new_observation.push_back(v);
      }
    }

    // remove all RR sets
    // save_rrsets_1(new_observation, hyperGT_2, hyperPos_2, hyperG_2);
    // save_rrsets_2(new_observation, hyperGT, hyperPos, hyperG);
    if (reuse)
    {
      getHittedRRsets(new_observation, R1);
      getHittedRRsets(new_observation, R2);
      //   clearRRsets(R1);
      clearRRsets(R3);
    }
    else
    {
      clearRRsets(R1);
      clearRRsets(R2);
      clearRRsets(R3);
    }
  }

/* continue to select seed node after it has found the node with largest influence*/
  void maxInfSingleton_extendSeedSet(rrCollection &R)
  {
    vector<int> coverage(n, 0);
    size_t maxDeg = 0;
    for (auto i = n; i--;) {
      // const auto deg = node_influence_1[i];
      const auto deg = R.invertedList[i].size();
      coverage[i] = deg - R.deduction[i];
      if (deg > maxDeg) maxDeg = deg;
    }

    assert(maxInfSeedSet.size() == 1);
    // check if an edge is removed
    std::vector<bool> edgeMark(R.hedges.size(), false);

    coverage[maxInfSeedSet[0]] = 0;
    for (auto edgeIdx : R.invertedList[maxInfSeedSet[0]]) {
      if (R.hedges[edgeIdx].size() == 0 || edgeMark[edgeIdx]) continue;
      edgeMark[edgeIdx] = true;
      for (auto nodeIdx : R.hedges[edgeIdx]) {
        if (coverage[nodeIdx] == 0) continue;  // This node is seed, skip
        coverage[nodeIdx]--;
      }
    }

    vector<vector<int>> degMap(maxDeg + 1);  // degMap: map degree to the nodes with this degree
    for (auto i = n; i--;) {
    //   if (coverage[i] == 0) continue;
      degMap[coverage[i]].push_back(i);
    }


    vector<int> sortedNode(n);  // sortedNode: record the sorted nodes in ascending order of degree
    vector<int> nodePosition(n);  // nodePosition: record the position of each node in the sortedNode
    vector<int> degreePosition(maxDeg + 2);  // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;
    for (auto& nodes : degMap) {
      degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
      idxDegree++;
      for (auto& node : nodes) {
        nodePosition[node] = idxSort;
        sortedNode[idxSort++] = node;
      }
    }

    /*
     * sortedNode: position -> node
     * nodePosition: node -> position
     * degreePosition: degree -> position (start position of this degree)
     * coverage: node -> degree
     * e.g., swap the position of a node with the start position of its degree
     * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
     */
    double total_cost = cost[maxInfSeedSet[0]];
    for (auto k = n; k--;) {
      const auto seed = sortedNode.back();
      sortedNode.pop_back();

      log_debug("potential seed: %d， cost: %f", seed, cost[seed]);
      if (total_cost + cost[seed] > budget) {
        log_debug("total cost: %f, seed size: %lu", total_cost, maxInfSeedSet.size());

        // avoid the situation that we encounter a node with large cost 
        return;
      }
      total_cost += cost[seed];
      maxInfSeedSet.push_back(seed);

      coverage[seed] = 0;
      for (auto edgeIdx : R.invertedList[seed]) {
        if (R.hedges[edgeIdx].size() == 0 || edgeMark[edgeIdx]) continue;
        edgeMark[edgeIdx] = true;
        for (auto nodeIdx : R.hedges[edgeIdx]) {
          if (coverage[nodeIdx] == 0) continue;           // This node is seed, skip
          const auto currPos = nodePosition[nodeIdx];     // The current position
          const auto currDeg = coverage[nodeIdx];         // The current degree
          const auto startPos = degreePosition[currDeg];  // The start position of this degree
          const auto startNode = sortedNode[startPos];    // The node with the start position
          // Swap this node to the start position with the same degree, and update their positions
          // in nodePosition
          std::swap(sortedNode[currPos], sortedNode[startPos]);
          nodePosition[nodeIdx] = startPos;
          nodePosition[startNode] = currPos;
          // Increase the start position of this degree by 1, and decrease the degree of this node
          // by 1
          degreePosition[currDeg]++;
          coverage[nodeIdx]--;
        }
      }
    }

    return;
  }

  double max_cover_lazy_budget(rrCollection &R)
  {
    // mode: optimization mode.
    // 0->no optimization,
    // 1->optimization with the upper bound in the last round,
    // 2->optimization with minimum upper bound among all rounds [Default].

    // cost influence most
    std::vector<double> coverage_cost(n, 0);
    int maxDeg_cost = 0;
    for (auto i = n; i--;)
    {
      const int deg_cost = R.invertedList[i].size() - R.deduction[i];
      coverage_cost[i] = (double)deg_cost / cost[i];
    }

    // Normalization
    for (auto i = n; i--;)
    {
      int deg_cost = int(coverage_cost[i]) + 1;
      if (deg_cost > maxDeg_cost)
      {
        maxDeg_cost = deg_cost;
      }
    }

    // int(u.inf/u.cost) -> sort(pair<u.id, u.inf/u.cost>)
    std::vector<std::vector<std::pair<RID_T, double>>> degMap_cost(maxDeg_cost + 1);
    for (auto i = n; i--;)
    {
      int deg_cost = int(coverage_cost[i]) + 1;
      degMap_cost[deg_cost].push_back(std::make_pair(i, coverage_cost[i]));
    }

    size_t sumInf_cost = 0;
    double sumBudget_cost = 0;
    double minUpCoverage = R.hedges.size();

    // check if an edge is removed
    std::vector<bool> edgeMark_cost(R.hedges.size(), false);
    std::vector<bool> node_sel_cost(n, false);

    nonAdaSeedSet.clear();

    bool exceed_budget = false;

    for (size_t deg = maxDeg_cost; deg > 0; deg--) // Enusre deg > 0
    {
      if (exceed_budget)
      {
        break;
      }

      std::sort(degMap_cost[deg].begin(), degMap_cost[deg].end(), [](auto &left, auto &right)
                { return std::get<1>(left) < std::get<1>(right); });
      auto &vecNode = degMap_cost[deg];
      for (auto idx = vecNode.size(); idx--;)
      {
        auto argmaxIdx = std::get<0>(vecNode[idx]);
        const size_t currDeg = int(coverage_cost[argmaxIdx]) + 1;
        if (deg > currDeg)
        {
          degMap_cost[currDeg].push_back(std::make_pair(argmaxIdx, coverage_cost[argmaxIdx]));
          degMap_cost[deg].pop_back();
          continue;
        }
        if (sumBudget_cost + cost[argmaxIdx] <= budget)
        {
          // mode
          if (2 == 2)
          {
            // Find upper bound
            auto degBound = deg;
            vector<uint32_t> vecBound;
            double vecBound_Budget = 0;
            // Initialize vecBound
            auto idxBound = idx + 1;
            bool is_stop = false;
  
            if (is_stop == false)
            {
              while (vecBound_Budget <= budget && idxBound--)
              {
                if (vecBound_Budget + cost[std::get<0>(degMap_cost[degBound][idxBound])] <= budget)
                {
                  vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                     cost[std::get<0>(degMap_cost[degBound][idxBound])]);
                  vecBound_Budget += cost[std::get<0>(degMap_cost[degBound][idxBound])];
                }
                else
                {
                  vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                     (budget - vecBound_Budget));
                  vecBound_Budget = budget;
                  is_stop = true;
                  break;
                }
              }
            }
            if (is_stop == false)
            {
              while (vecBound_Budget <= budget && --degBound)
              {
                std::sort(degMap_cost[degBound].begin(), degMap_cost[degBound].end(), [](auto &left, auto &right)
                          { return std::get<1>(left) < std::get<1>(right); });
                idxBound = degMap_cost[degBound].size();
                while (vecBound_Budget <= budget && idxBound--)
                {
                  {
                    if (vecBound_Budget + cost[std::get<0>(degMap_cost[degBound][idxBound])] <= budget)
                    {
                      vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                         cost[std::get<0>(degMap_cost[degBound][idxBound])]);
                      vecBound_Budget += cost[std::get<0>(degMap_cost[degBound][idxBound])];
                    }
                    else
                    {
                      vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                         (budget - vecBound_Budget));
                      vecBound_Budget = budget;
                      is_stop = true;
                      break;
                    }
                  }
                }
                if (is_stop == true)
                {
                  break;
                }
              }
            }
            std::sort(vecBound.begin(), vecBound.end());
            // make_min_heap(vecBound);

            double UpCoverage = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf_cost) ;
            minUpCoverage = min(minUpCoverage, UpCoverage);
          }
          sumInf_cost = sumInf_cost + coverage_cost[argmaxIdx] * cost[argmaxIdx];

          nonAdaSeedSet.push_back(argmaxIdx);
          sumBudget_cost = sumBudget_cost + cost[argmaxIdx];
          coverage_cost[argmaxIdx] = 0;
          node_sel_cost[argmaxIdx] = true;

          // min_cost = min_vector_nonzero_mask(cost, node_sel_cost);
          for (auto edgeIdx : R.invertedList[argmaxIdx])
          {
            if (edgeMark_cost[edgeIdx])
            {
              continue;
            }
            edgeMark_cost[edgeIdx] = true;
            for (auto nodeIdx : R.hedges[edgeIdx])
            {
              if (coverage_cost[nodeIdx] == 0)
              {
                continue; // This node is seed, skip
              }
              coverage_cost[nodeIdx] = (double)(coverage_cost[nodeIdx] * cost[nodeIdx] - 1) /
                                       cost[nodeIdx];
            }
          }
          degMap_cost[deg].pop_back();
          for (auto degup = deg; degup > 0; degup--)
          {
            auto &vecNodeup = degMap_cost[degup];
            for (auto idxup = vecNodeup.size(); idxup--;)
            {
              std::get<1>(degMap_cost[degup][idxup]) = coverage_cost[std::get<0>(degMap_cost[degup][idxup])];
            }
          }
          std::sort(degMap_cost[deg].begin(), degMap_cost[deg].end(), [](auto &left, auto &right)
                    { return std::get<1>(left) < std::get<1>(right); });
        }
        else
        {
          log_debug("exceeding budget");
          exceed_budget = true;
          break;
        }
      }
      degMap_cost.pop_back();
    }

    log_debug("budget used: %f, size: %lu", sumBudget_cost, nonAdaSeedSet.size());
    log_debug("coverage: %lu, upperCoverage: %f", sumInf_cost, minUpCoverage);

    return minUpCoverage;
  }

  void CASE_selectNode(int &seed, rrCollection &R, double a, double &upperRatio)
  {
    double maxRatio = 0;
    double maxUpperRatio = 0;
    uint32_t maxRatioNode = 0;
    uint32_t coverage = 0;

    for (auto i = n; i--;)
    {
      if (active_set[i] || cost[i] > budget)
        continue;
      coverage = R.invertedList[i].size() - R.deduction[i];
      double ratio = (double)coverage / cost[i];
      // double upperRatio = ratio;

    //   double upperRatio = sqr(sqrt(coverage + a / 2.) + sqrt(a / 2.)) / cost[i];
      double upperRatio = (double)coverage / cost[i];
      if (maxRatio < ratio) {
        maxRatio = ratio;
        maxRatioNode = i;
      }

      maxUpperRatio = max(maxUpperRatio, upperRatio);
    }

    seed = maxRatioNode;
    upperRatio = maxUpperRatio;

    return;
  }

  void maxInfSingleton_selectNode(vector<int> &batch_set, rrCollection &R, double a, double &upperInf)
  {
    double maxInf = 0;
    uint32_t maxInfNode = 0;
    uint32_t coverage = 0;

    for (auto i = n; i--;)
    {
      if (active_set[i] || cost[i] > budget)
        continue;

      // const auto deg = node_influence_1[i];
      coverage = R.invertedList[i].size() - R.deduction[i];
      double inf = (double)coverage;

      if (maxInf < inf)
      {
        maxInf = inf;
        maxInfNode = i;
      }
    }

    batch_set.clear();
    batch_set.push_back(maxInfNode);

    // upperInf = sqr(sqrt(maxInf + a / 2.) + sqrt(a / 2.));
    upperInf = maxInf;
    return;
  }

  double build_seedset(int targetSize, vector<int> &batch_set, rrCollection &R)
  {
    vector<int> coverage(n, 0);
    size_t maxDeg = 0;
    for (auto i = n; i--;)
    {
      // const auto deg = node_influence_1[i];
      const auto deg = R.invertedList[i].size();
      coverage[i] = deg - R.deduction[i];
      if (deg > maxDeg)
        maxDeg = deg;
    }

    vector<vector<int>> degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
    for (auto i = n; i--;)
    {
      // if (coverage[i] == 0) continue;
      degMap[coverage[i]].push_back(i);
    }
    vector<int> sortedNode(n); // sortedNode: record the sorted nodes in ascending order of degree
    vector<int> nodePosition(
        n); // nodePosition: record the position of each node in the sortedNode
    vector<int> degreePosition(
        maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;
    for (auto &nodes : degMap)
    {
      degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
      idxDegree++;
      for (auto &node : nodes)
      {
        nodePosition[node] = idxSort;
        sortedNode[idxSort++] = node;
      }
    }
    // check if an edge is removed
    std::vector<bool> edgeMark(R.hedges.size(), false);
    // record the total of top-k marginal gains
    int sumTopk = 0;
    for (auto deg = maxDeg + 1; deg--;)
    {
      if (degreePosition[deg] <= (int)(n - targetSize))
      {
        sumTopk += deg * (degreePosition[deg + 1] - (int)(n - targetSize));
        break;
      }
      sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }
    double boundMin = 1.0 * sumTopk;
    double sumInf = 0;

    batch_set.clear();

    /*
     * sortedNode: position -> node
     * nodePosition: node -> position
     * degreePosition: degree -> position (start position of this degree)
     * coverage: node -> degree
     * e.g., swap the position of a node with the start position of its degree
     * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
     */
    for (auto k = targetSize; k--;)
    {
      const auto seed = sortedNode.back();
      sortedNode.pop_back();
      int newNumV = sortedNode.size();
      sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
      sumInf += coverage[seed];
      batch_set.push_back(seed);
      coverage[seed] = 0;
      for (auto edgeIdx : R.invertedList[seed])
      {
        if (R.hedges[edgeIdx].size() == 0 || edgeMark[edgeIdx])
          continue;
        edgeMark[edgeIdx] = true;
        for (auto nodeIdx : R.hedges[edgeIdx])
        {
          if (coverage[nodeIdx] == 0)
            continue;                                    // This node is seed, skip
          const auto currPos = nodePosition[nodeIdx];    // The current position
          const auto currDeg = coverage[nodeIdx];        // The current degree
          const auto startPos = degreePosition[currDeg]; // The start position of this degree
          const auto startNode = sortedNode[startPos];   // The node with the start position
          // Swap this node to the start position with the same degree, and update their positions
          // in nodePosition
          std::swap(sortedNode[currPos], sortedNode[startPos]);
          nodePosition[nodeIdx] = startPos;
          nodePosition[startNode] = currPos;
          // Increase the start position of this degree by 1, and decrease the degree of this node
          // by 1
          degreePosition[currDeg]++;
          coverage[nodeIdx]--;
          // If the start position of this degree is in top-k, reduce topk by 1
          if (startPos >= newNumV - targetSize)
            sumTopk--;
        }
      }
      const auto boundLast = 1.0 * (sumInf + sumTopk);
      if (boundMin > boundLast)
        boundMin = boundLast;
    }
    // cout << batch_set.size() << ": boundMin: " << boundMin << endl;// " boundLast: " << boundLast
    // << endl;// << " sumTopk: " << sumTopk << " sum: " << sumInf + sumTopk << endl;
    return 1.0 * boundMin;
  }

  double Coverage(vector<int> &batch_set, rrCollection &R)
  {

    vector<bool> overlap(R.hedges.size(), false);
    RID_T flagSize = R.hittedFlagVec.size();

    if (lazyInvList)
    {
      for (auto seed : batch_set)
      {
        for (auto rrIdx : R.invertedList[seed])
        {

          if (rrIdx < flagSize && R.hittedFlagVec[rrIdx])
            continue;

          overlap[rrIdx] = true; // rrIdx is the index of the rr set that contains seedIdx.
        }
      }
    }
    else
    {
      for (auto seed : batch_set)
      {
        for (auto rrIdx : R.invertedList[seed])
        {
          overlap[rrIdx] = true; // rrIdx is the index of the rr set that contains seedIdx.
        }
      }
    }

    // assert(batch_set.size() == 1);
    uint32_t coverage = count(overlap.begin(), overlap.end(), true);
    log_debug("coverage: %u", coverage);

    return coverage;
  }

  inline double randomNumber()
  {
    return sfmt_genrand_real1(&sfmtSeed);
  }
};