#ifndef __HEAD_H__
#define __HEAD_H__

#ifdef DEBUG_INFO
this_is_deprecated
#endif
#ifdef DEBUG_TRACE
this_is_deprecated
#endif


#if defined(WIN32)
#elif defined(__CYGWIN__) // cygwin
#include <sys/time.h>
#else //linux
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif



#include <iostream>
#include <set>
#include <list>
#include <sstream>
#include <cmath>
#include <queue>
#include <fstream>
#include <string>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <climits>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <deque>
#include <chrono>
#include <ctime>
#include <ratio>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unordered_set>
#include <set>
#include <assert.h>

using namespace std;
// typedef unsigned int uint;
// typedef unsigned char uint8_t;
// typedef long long int64_t;
// typedef unsigned long long uint64_t;
// typedef pair<int,int> ipair;
// typedef pair<double,double> dpair;
#define MP make_pairint64_t
#define F first
#define S second

#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (2393.910e6)
#endif

typedef uint64_t RID_T;

// typedef char int8_t;
// typedef unsigned char uint8_t;
// typedef long long int64_t;
// typedef unsigned long long uint64_t;

typedef pair<int64_t, int64_t> offset_t;
#define THD_NUM 32

#define BASIC_COST 1.0f
#define RANDOM_MAX 10.0f

inline double pow2(const double t)
{
    return t * t;
}

//LINUX
uint64_t rdtsc(void)
{
    unsigned a, d;
    //asm("cpuid");
    asm volatile("rdtsc" : "=a" (a), "=d" (d));
    return (((uint64_t)a) | (((uint64_t)d) << 32));
}


string nowStr(){
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char str[100];
    sprintf(str,"%d_%d_%d_%d_%d_%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    return string(str);
}


map<string, timeval> __head_h_start;
void timer_init(string arg="default"){
    timeval ts;
    gettimeofday(&ts,NULL);
    __head_h_start[arg]=ts;
}
int64_t timer_elapse(string arg="default"){ // unit ms
    cout << "fuck " <<endl;
    struct timeval now;
    gettimeofday(&now,NULL);
    int64_t sec=now.tv_sec-__head_h_start[arg].tv_sec;
    int64_t usec=now.tv_usec-__head_h_start[arg].tv_usec;
    return sec * 1000 + usec/1000;
}
// #define SIZE(t) (int)(t.size())
// #define ALL(t) (t).begin(), (t).end()
// #define FOR(i,n) for(int (i)=0; (i)<((int)(n)); (i)++)
// #ifdef WIN32
// #define FORE(i, x) for (typeid((x).begin()) i = (x).begin(); (i) != (x).end(); (i)++)
// #else
// #define FORE(i, x) for (__typeof((x).begin()) i = (x).begin(); (i) != (x).end(); (i)++)
// #endif


// static inline string &ltrim(string &s) { s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace)))); return s; }
// static inline string &rtrim(string &s) { s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end()); return s; }
// static inline string &trim(string &s) { return ltrim(rtrim(s)); }
// string __n_variable(string t, int n){ t=t+','; int i=0; if(n) for(; i<SIZE(t)&&n; i++) if(t[i]==',') n--; n=i; for(;t[i]!=',';i++); t=t.substr(n, i-n); trim(t);if(t[0]=='"') return ""; return t+"="; }
// #define __expand_nv(x) __n_variable(t, x)<< t##x << " "
// template<class T0>
// void ___debug(string t,deque<T0> t0,ostream&os){os<<__n_variable(t,0);FOR(i, SIZE(t0))os<<t0[i]<<" ";}
// template<class T0>
// void ___debug(string t,set<T0> t0,ostream&os){os<<__n_variable(t,0);FORE(i,t0)os<<*i<<" ";}
// template<class T0>
// void ___debug(string t,vector<T0> t0,ostream&os){os<<__n_variable(t,0);FOR(i, SIZE(t0))os<<t0[i]<<" ";}
// template<class T0,class T1>
// void ___debug(string t,vector<pair<T0,T1> > t0,ostream&os){os<<__n_variable(t,0);FOR(i, SIZE(t0))os<<t0[i].F<<","<<t0[i].S<<" ";}
// template<class T0>
// void ___debug(string t,T0 t0,ostream&os){os<<__expand_nv(0);}
// template<class T0,class T1>
// void ___debug(string t,T0 t0,T1 t1,ostream&os){os<<__expand_nv(0)<<__expand_nv(1);}
// template<class T0,class T1,class T2>
// void ___debug(string t,T0 t0,T1 t1,T2 t2,ostream&os){os<<__expand_nv(0)<<__expand_nv(1)<<__expand_nv(2);}
// template<class T0,class T1,class T2,class T3>
// void ___debug(string t,T0 t0,T1 t1,T2 t2,T3 t3,ostream&os){os<<__expand_nv(0)<<__expand_nv(1)<<__expand_nv(2)<<__expand_nv(3);}
// template<class T0,class T1,class T2,class T3,class T4>
// void ___debug(string t,T0 t0,T1 t1,T2 t2,T3 t3, T4 t4,ostream&os){os<<__expand_nv(0)<<__expand_nv(1)<<__expand_nv(2)<<__expand_nv(3)<<__expand_nv(4);}


//#define DO_ONCE


#define RUN_TIME(...) { int64_t t=rdtsc();  __VA_ARGS__; t=rdtsc()-t; cout<<  #__VA_ARGS__ << " : " << t/TIMES_PER_SEC <<"s"<<endl;  }

#define TRACE(...) ;
#define IF_TRACE(args) ;
#define TRACE_LINE(...) ;
#define TRACE_SKIP(a, ...) ;
#define TRACE_LINE_SKIP(a, ...) ;
#define TRACE_LINE_END(...) ;
#define TRACE_LOG(...) ;




#define ASSERTT(v, ...) ;
#define ASSERT(v ) ;
#define INFO(...) ;

enum { ADAP_CA_GREEDY = 1, NONADAP_CA_GREEDY = 2, MAX_INF_SINGLETON = 3};

enum {SUMMARY, VANILLA_RRSETS, SUBSIM_RRSETS, MULTITHD_PREP, MULTITHD_INVLIST, MULTITHD_SAMPLE};
class Timer
{
    public:
        static vector<double> timeUsed;
        static vector<string> timeUsedDesc;
        static vector<chrono::high_resolution_clock::time_point> procPauseStart;
        static vector<double> procPauseInterval;


        chrono::duration<double> interval;
        int id;
        //uint64_t startTime;
        chrono::high_resolution_clock::time_point startTime;
        chrono::high_resolution_clock::time_point endTime;
        bool showOnDestroy;
        Timer(int id, string desc="", bool showOnDestroy=false)
        {
            this->id=id;
            while((int)timeUsed.size() <= id)
            {
                timeUsed.push_back(0);
                timeUsedDesc.push_back("");
            }
            timeUsedDesc[id]=desc;
            //startTime=rdtsc();
            startTime=chrono::high_resolution_clock::now();

            this->showOnDestroy=showOnDestroy;
        }
        ~Timer()
        {
            endTime=chrono::high_resolution_clock::now();
            interval= chrono::duration_cast<chrono::duration<double>>(endTime-startTime);

            

            //cout<<"id="<<id<<" Running time="<<interval.count()<<endl;
            //timeUsed[id]+=(rdtsc()-startTime);
            timeUsed[id]+=((double)interval.count() - getPauseInterval(startTime));
            //cout<<"id="<<id<<" "<<timeUsed[id]<<endl;
        
        }

        double getPauseInterval(chrono::high_resolution_clock::time_point timerStart)
        {
            double pauseInterval = 0;
            assert(procPauseStart.size() == procPauseInterval.size());
            for (size_t index = 0; index < procPauseStart.size(); index++)
            {
                if (timerStart < procPauseStart[index])
                {
                    pauseInterval+= procPauseInterval[index];
                }
            }

            return pauseInterval;
        }

        static double getSpecificInterval(int id) { return timeUsed[id]; }

        static void pauseHandler()
        {
            procPauseStart.push_back(chrono::high_resolution_clock::now());
        }

        static void resumeHandler()
        {
            if (procPauseStart.size() == procPauseInterval.size())
            {
              printf("pause vector size: %lu\n", procPauseStart.size());
              return;
            }
            
            int index = procPauseStart.size() -1;
            assert(index >= 0);
            auto endPauseTime =  chrono::high_resolution_clock::now();
            auto pauseInterval= chrono::duration_cast<chrono::duration<double>>(endPauseTime-procPauseStart[index]);
            procPauseInterval.push_back((double)pauseInterval.count());
            printf("current interval: %f\n", procPauseInterval[index]);
            assert(procPauseStart.size() == procPauseInterval.size());
        }


        static void show(int time, bool debug=false)
        {
            //INFO("########## Timer ##########");
 
            for (int i=0;i<(int)timeUsed.size();i++) 
            {
            if (timeUsed[i]>0)
            {
                char str[100];
                //sprintf(str,"%.6lf",timeUsed[i]/TIMES_PER_SEC );
                sprintf(str,"%.6lf", timeUsed[i]/time);
                string s=str;
                if ((int)s.size()<15) s=" "+s;
                char t[100];
                memset(t, 0, sizeof t);
                sprintf(t,"Spend %s seconds on %s",s.c_str(), timeUsedDesc[i].c_str());
                cout<< t << endl;
            }
            }
        }
        static void clearAll()
        {
            timeUsed.clear();
            timeUsedDesc.clear();
        }
};
vector<double> Timer::timeUsed;
vector<string> Timer::timeUsedDesc;
vector<chrono::high_resolution_clock::time_point> Timer::procPauseStart;
vector<double> Timer::procPauseInterval;




string __head_version = "";
class OutputInfo
{
    public:
    OutputInfo(int argn, char ** argv){
        //cout<<"\e\[0;32mProgram Start at: " << currentTimestampStr()<<"\e[0m"<<endl;
        //cout<<"\e\[0;32mProgram version: " << __head_version << "\e[0m"<<endl;
        /*
        cout<<"Arguments: ";
        for(int i=0; i<argn; i++){
            cout<<argv[i]<<" ";
        }
        cout<<endl;
        cout<<"--------------------------------------------------------------------------------" <<endl;
        */
    }
    ~OutputInfo(){
        //cout<<"\e\[0;31mProgram Terminate at: " << currentTimestampStr()<< "\e[0m"<<endl;
        //cout<<"\e\[0;31mProgram version: " << __head_version << "\e[0m"<<endl;
    }
};

#endif //__HEAD_H__
