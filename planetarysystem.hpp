/* 
 * File:   planetarysystem.hpp
 * Author: user
 *
 * Created on 29. August 2008, 16:34
 */

#ifndef _PLANETARYSYSTEM_HPP
#define	_PLANETARYSYSTEM_HPP

//----------------------------------------------------------------------------------------
#include <string>
#include <vector>
#include <sys/time.h>

//----------------------------------------------------------------------------------------
#include "simtypes.hpp"

//----------------------------------------------------------------------------------------

class PlanetarySystem
{
public:
    
    PlanetarySystem();
    
    void InitTracerGrid(int num, myfloat xmin, myfloat xmax, myfloat ymin, myfloat ymax);
    void InitTracerRandom(int num, myfloat rad);
    void InitTracerResonance();
    void InitFromFile(const char *szFile);
    void SetTracerTo(myfloat x, myfloat y);

    void Move();
    int GetPlanetDataArrays(myfloat **x, myfloat **y, myfloat **rad, myfloat **mass, int **col);
    int GetTracerDaraArrays(myfloat **x, myfloat **y, myfloat **vx, myfloat **vy);
    int GetStatData(myfloat **exc);
    std::size_t GetNumPlanets() const;
    std::size_t GetNumActive() const;
    std::size_t GetNumSteps() const;
    void GetOrbitalVelocity(int planet, myfloat x, myfloat y, myfloat &vx, myfloat &vy);
    myfloat GetFPS() const;
    myfloat GetTime() const;
    void SetCenter(int idx);
    void DoTracerReset(bool bStat);

private:

    void ParseDefine(std::ifstream &ifs);
    void ParseParam(std::ifstream &ifs);
    void MoveVerlet();
    void MovePlanets();
    void CheckBenchExit();

    void ResetTracer(std::size_t i);
    void InitArrays(int num);
    
    static const myfloat DEG_TO_RAD;
    static const myfloat GAMMA_;   // m³/(kg*s²)
    
    std::size_t m_nNumPlanets;
    std::size_t m_nNumTracer;
    std::size_t m_nNumActiveTracer;
    bool m_bDoReset;
    myfloat m_xmin;
    myfloat m_xmax;
    myfloat m_ymin;
    myfloat m_ymax;
    myfloat m_dt;
    myfloat m_fps;
    myfloat m_time;
    int m_center;
    int m_frames;
    struct timeval m_tv1, m_tv2;
    
    // planet data
    std::vector<myfloat> m_plPos[2];
    std::vector<myfloat> m_plPosOld[2]; 
    std::vector<myfloat> m_plVel[2];
    std::vector<myfloat> m_plAcc[2];
    std::vector<myfloat> m_plRad;       // Planets radius in km
    std::vector<myfloat> m_plMass;      // Planetary Mass in kg
    std::vector<std::string> m_plName;  // The planets name
    std::vector<int> m_plColor;
    
    // tracer positions and velocities
    std::vector<myfloat> m_trPos[2];
    std::vector<myfloat> m_trPosOld[2];
    std::vector<myfloat> m_trAcc[2];
    std::vector<myfloat> m_trMinDist;
    std::vector<myfloat> m_trMaxDist;
    std::vector<myfloat> m_trExc;
    std::vector<int> m_trColor;
};

#endif	/* _PLANETARYSYSTEM_HPP */

