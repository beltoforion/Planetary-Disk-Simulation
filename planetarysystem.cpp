#include "planetarysystem.hpp"
#include "planetrender.hpp"

#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sys/time.h>
#include <omp.h>
#include <fstream>
#include <stdexcept>

const myfloat PlanetarySystem::DEG_TO_RAD = 0.0174532925;
const myfloat PlanetarySystem::GAMMA_ = 6.67428e-20; // km³/(kg*s²)

//----------------------------------------------------------------------------------------

PlanetarySystem::PlanetarySystem()
: m_nNumPlanets(0)
, m_nNumTracer(0)
, m_nNumActiveTracer(0)
, m_bDoReset(true)
, m_xmin(0)
, m_xmax(0)
, m_ymin(0)
, m_ymax(0)
, m_dt(0)
, m_fps(0)
, m_time(0)
, m_center(0)
, m_frames(0)
{
  gettimeofday(&m_tv1, NULL);
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::SetCenter(int idx)
{
  m_center = idx;
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::InitArrays(int num)
{
  m_nNumTracer = num;
  m_nNumActiveTracer = num;

  for (std::size_t i = 0; i < 2; ++i)
  {
    m_trPos[i].resize(num, 0);
    m_trPosOld[i].resize(num, 0);
    m_trAcc[i].resize(num, 0);
  }
  m_trColor.resize(num, 0);
  m_trMinDist.resize(num, 9e29);
  m_trMaxDist.resize(num, 0);
  m_trExc.resize(num, 0);
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::DoTracerReset(bool bStat)
{
  m_bDoReset = bStat;
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::ParseParam(std::ifstream &ifs)
{
  const char *szKeys[] = {"DELTA_T",
                          "BENCHMARK_EXIT",
                          "END",
                          NULL
  };
  char szLine[1024];

  std::cout << "Parsing simulation parameters\n";

  while (!ifs.getline(szLine, sizeof (szLine)).eof())
  {
    if (szLine[0] == '#')
      continue;

    for (int i = 0; szKeys[i]; ++i)
    {
      if (!strstr(szLine, szKeys[i]))
        continue;

      switch (i)
      {
        case 0: // DELTA_T
        {
          double buf(0);
          sscanf(szLine, "%*[^=]=%lf", &buf);
          m_dt = buf;
        }
          break;

        case 2: // END
          std::cout << "  DELTA_T = " << m_dt << "\n";
          std::cout << std::endl;
          return;
      } // switch keyword
    } // for all keywords
  } // while read line
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::ParseDefine(std::ifstream &ifs)
{
  const char *szKeys[] = {"NAME",
                          "RAD",
                          "ORBIT",
                          "POSITION",
                          "VELOCITY",
                          "MASS",
                          "COLOR",
                          "END",
                          NULL
  };

  char szLine[1024], szBuf[1024];
  int stat(0);
  std::size_t idx = m_nNumPlanets;

  std::cout << "Parsing planet/star/moon definition\n";

  m_plPos[0].push_back(0);
  m_plPos[1].push_back(0);
  m_plPosOld[0].push_back(0);
  m_plPosOld[1].push_back(0);
  m_plVel[0].push_back(0);
  m_plVel[1].push_back(0);
  m_plAcc[0].push_back(0);
  m_plAcc[1].push_back(0);
  m_plMass.push_back(0);
  m_plName.push_back("");
  m_plColor.push_back(0);
  m_plRad.push_back(0);

  while (!ifs.getline(szLine, sizeof (szLine)).eof())
  {
    if (szLine[0] == '#')
      continue;

    for (int i = 0; szKeys[i]; ++i)
    {
      if (!strstr(szLine, szKeys[i]))
        continue;

      switch (i)
      {
        case 0: // NAME
          stat = sscanf(szLine, "%*[^=]=%s", szBuf);
          if (stat != 1)
            throw std::runtime_error("Parsing error \"NAME\" line.");
          m_plName.back() = szBuf;
          break;
        case 1: // RAD
        {
          double rad(0);
          stat = sscanf(szLine, "%*[^=]=%lf", &rad);
          if (stat != 1)
            throw std::runtime_error("Parsing error \"RAD\" line.");
          m_plRad.back() = rad;
        }
          break;
        case 2: // ORBIT
        {
          double elon(0), dist(0);
          std::size_t k(0);

          sscanf(szLine, "%*[^=]=%*[ ]%[a-zA-Z],%lf,%lf", szBuf, &elon, &dist);
          std::cout << "  ORBIT = \"" << szBuf << "\", " << elon << ", " << dist << " km\n";
          myfloat x(dist * sin(elon * DEG_TO_RAD)),
                  y(dist * cos(elon * DEG_TO_RAD)), vx(0), vy(0);

          for (k = 0; k < m_plName.size(); ++k)
          {
            if (!strcmp(m_plName[k].c_str(), szBuf))
            {
              GetOrbitalVelocity(k, x, y, vx, vy);
              m_plPos[0].back() = x;
              m_plPos[1].back() = y;
              m_plVel[0].back() = vx;
              m_plVel[1].back() = vy;
            }
          }

          if (k == m_plName.size() + 1)
            throw std::runtime_error("DEFINE/ORBIT: invalid orbit reference.");

          break;
        }
        case 3: // POSITION
        {
          double x(0), y(0);
          sscanf(szLine, "%*[^=]=%*[ ]%lf,%lf", &x, &y);
          m_plPos[0].back() = x;
          m_plPos[1].back() = y;
        }
          break;
        case 4: // VELOCITY
        {
          double x(0), y(0);
          sscanf(szLine, "%*[^=]=%*[ ]%lf,%lf", &x, &y);
          m_plVel[0].back() = x;
          m_plVel[1].back() = y;
        }
          break;
        case 5: // MASS
        {
          double mass(0);
          sscanf(szLine, "%*[^=]=%*[ ]%lf", &mass);
          m_plMass.back() = mass;
        }
          break;
        case 6: // COLOR
          sscanf(szLine, "%*[^=]=%*[ ]%x", &m_plColor.back());
          break;
        case 7: // END
          std::cout << "  NAME = \"" << m_plName.back() << "\"\n";
          std::cout << "  RAD = " << m_plRad.back() << " km\n";
          std::cout << "  POS = " << m_plPos[0].back() << ", " << m_plPos[1].back() << " km\n";
          std::cout << "  VEL = " << m_plVel[0].back() << ", " << m_plVel[1].back() << " km/s\n";
          std::cout << "  MASS = " << m_plMass.back() << " kg\n";
          std::cout << "  COLOR = " << std::hex << m_plColor.back() << "\n";
          std::cout << std::endl;
          m_nNumPlanets = m_plMass.size();
          return;
      }
    } // for all keywords
  } // while read line
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::InitFromFile(const char *szFile)
{
  assert(szFile);

  char szBuf[1024];
  const char *szKeys[] = {"DEFINE",
                          "PARAM",
                          "TRACER",
                          NULL
  };

  std::ifstream ifs(szFile);
  while (!ifs.getline(szBuf, sizeof (szBuf)).eof())
  {
    if (szBuf[0] == '#')
      continue;
    
    for (int i = 0; szKeys[i]; ++i)
    {
      if (!strstr(szBuf, szKeys[i]))
        continue;

      switch (i)
      {
        case 0: // DEFINE
          ParseDefine(ifs);
          break;

        case 1: // PARAM
          ParseParam(ifs);
          break;

        case 2: // TRACER
          //ParseTracer(ifs);
          break;
      }
    } // for all keywords
  } // while not eof

  // Calculate old positions
  for (std::size_t i = 0; i < m_nNumPlanets; ++i)
  {
    m_plPosOld[0][i] = m_plPos[0][i] - m_dt * m_plVel[0][i];
    m_plPosOld[1][i] = m_plPos[1][i] - m_dt * m_plVel[1][i];
  }
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::InitTracerGrid(int num, myfloat xmin, myfloat xmax, myfloat ymin, myfloat ymax)
{
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::InitTracerRandom(int num,
                                       myfloat rad)
{
  InitArrays(num);

  myfloat dist, elon;
  myfloat vel[2];
  myfloat sun_spc = m_plRad[0] * 20;
  rad = rad - sun_spc; // do not place directly beside sun

  for (int i = 0; i < num; ++i)
  {
    dist = sun_spc + rad * ((myfloat) rand() / RAND_MAX);
    elon = 2 * M_PI * ((myfloat) rand() / RAND_MAX);
    m_trPos[0][i] = m_plPos[0][0] + dist * cos(elon);
    m_trPos[1][i] = m_plPos[1][0] + dist * sin(elon);

    GetOrbitalVelocity(0, m_trPos[0][i], m_trPos[1][i], vel[0], vel[1]);
    m_trPosOld[0][i] = m_trPos[0][i] - m_dt * vel[0];
    m_trPosOld[1][i] = m_trPos[1][i] - m_dt * vel[1];

    m_trColor[i] = 0xffffffff;
  }

  m_xmin = m_plPos[0][0] - rad;
  m_xmax = m_plPos[0][0] + rad;
  m_ymin = m_plPos[1][0] - rad;
  m_ymax = m_plPos[1][0] + rad;
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::InitTracerResonance()
{
  // 3 rings of dust 
  // before, at and behint the 2:1 resonance of jupiter

  for (std::size_t i = 0; i < 360; ++i)
  {
  }

}

//----------------------------------------------------------------------------------------

std::size_t PlanetarySystem::GetNumPlanets() const
{
  return m_nNumPlanets;
}

//----------------------------------------------------------------------------------------

std::size_t PlanetarySystem::GetNumActive() const
{
  return m_nNumActiveTracer;
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::ResetTracer(std::size_t i)
{
  // copy the last active tracer to position i
  if (!m_bDoReset)
  {
    m_nNumActiveTracer--;
    m_trPos[0][i] = m_trPos[0][m_nNumActiveTracer];
    m_trPos[1][i] = m_trPos[1][m_nNumActiveTracer];
    m_trPosOld[0][i] = m_trPosOld[0][m_nNumActiveTracer];
    m_trPosOld[1][i] = m_trPosOld[1][m_nNumActiveTracer];
    m_trAcc[0][i] = m_trAcc[0][m_nNumActiveTracer];
    m_trAcc[1][i] = m_trAcc[1][m_nNumActiveTracer];
    m_trColor[i] = m_trColor[m_nNumActiveTracer];
    m_trColor[i] = m_trColor[m_nNumActiveTracer];

    m_trMinDist[i] = m_trMinDist[m_nNumActiveTracer];
    m_trMaxDist[i] = m_trMaxDist[m_nNumActiveTracer];
    m_trExc[i] = m_trExc[m_nNumActiveTracer];
  }
  else
  {
    myfloat xrng = m_xmax - m_xmin;
    myfloat yrng = m_ymax - m_ymin;
    myfloat vel[2];

    m_trPos[0][i] = m_xmin + (xrng * ((myfloat) rand() / RAND_MAX));
    m_trPos[1][i] = m_ymin + (yrng * ((myfloat) rand() / RAND_MAX));

    GetOrbitalVelocity(0, m_trPos[0][i], m_trPos[1][i], vel[0], vel[1]);
    m_trPosOld[0][i] = m_trPos[0][i] - m_dt * vel[0];
    m_trPosOld[1][i] = m_trPos[1][i] - m_dt * vel[1];

    m_trAcc[0][i] = 0;
    m_trAcc[1][i] = 0;

    m_trMinDist[i] = 9e19;
    m_trMaxDist[i] = -9e19;
    m_trExc[i] = 1;
  }
}

//----------------------------------------------------------------------------------------

myfloat PlanetarySystem::GetTime() const
{
  return m_time;
}

//----------------------------------------------------------------------------------------

myfloat PlanetarySystem::GetFPS() const
{
  return m_fps;
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::Move()
{
  MovePlanets();
  MoveVerlet();

  gettimeofday(&m_tv2, NULL);

  m_time += m_dt;
  ++m_frames;
  myfloat t1 = (myfloat) (m_tv1.tv_sec * 1000 + (m_tv1.tv_usec / 1000)),
          t2 = (myfloat) (m_tv2.tv_sec * 1000 + (m_tv2.tv_usec / 1000)),
          dt = t2 - t1;

  if (dt > 1000)
  {
    gettimeofday(&m_tv1, NULL);
    m_fps = (myfloat) m_frames / (dt * 0.001);
    m_frames = 0;
  }
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::MovePlanets()
{
  // Calculate accellerations
  myfloat r[2], buf[2];
  myfloat * acc[2], *pos[2], *pos_old[2];

  acc[0] = &(m_plAcc[0][0]);
  acc[1] = &(m_plAcc[1][0]);
  pos[0] = &(m_plPos[0][0]);
  pos[1] = &(m_plPos[1][0]);
  pos_old[0] = &(m_plPosOld[0][0]);
  pos_old[1] = &(m_plPosOld[1][0]);

  // reset accellerations
  memset(acc[0], 0, sizeof (myfloat) * m_nNumPlanets);
  memset(acc[1], 0, sizeof (myfloat) * m_nNumPlanets);

  // calculate accellerations
  for (int i = 0; i < (int) m_nNumPlanets; ++i)
  {
    for (int j = i + 1; j < (int) m_nNumPlanets; ++j)
    {
      r[0] = pos[0][i] - pos[0][j];
      r[1] = pos[1][i] - pos[1][j];
      myfloat dist = sqrt(r[0] * r[0] + r[1] * r[1]);

      myfloat frc_x = -GAMMA_ * m_plMass[i] * m_plMass[j] / (dist * dist * dist) * r[0],
              frc_y = -GAMMA_ * m_plMass[i] * m_plMass[j] / (dist * dist * dist) * r[1];

      acc[0][i] += frc_x / m_plMass[i];
      acc[1][i] += frc_y / m_plMass[i];
      acc[0][j] -= frc_x / m_plMass[j];
      acc[1][j] -= frc_y / m_plMass[j];

    }
  }

  // integrate
  for (int i = 0; i < (int) m_nNumPlanets; ++i)
  {
    buf[0] = pos[0][i];
    buf[1] = pos[1][i];
    pos[0][i] = 2 * pos[0][i] - pos_old[0][i] + acc[0][i] * m_dt * m_dt;
    pos[1][i] = 2 * pos[1][i] - pos_old[1][i] + acc[1][i] * m_dt * m_dt;
    pos_old[0][i] = buf[0];
    pos_old[1][i] = buf[1];
  }
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::MoveVerlet()
{
  // Calculate accellerations
  myfloat r[2], buf[2];
  myfloat * acc[2], *pos[2], *pos_old[2];

  acc[0] = &(m_trAcc[0][0]);
  acc[1] = &(m_trAcc[1][0]);
  pos[0] = &(m_trPos[0][0]);
  pos[1] = &(m_trPos[1][0]);
  pos_old[0] = &(m_trPosOld[0][0]);
  pos_old[1] = &(m_trPosOld[1][0]);
  int i;

omp_set_num_threads(4);

#pragma omp parallel for default(shared) private(buf, r)
  for (i = 0; i < (int) m_nNumActiveTracer; ++i)
  {
    acc[0][i] = 0;
    acc[1][i] = 0;

    myfloat dist0 = 0; // distance tracer-sun

    for (std::size_t j = 0; j < m_nNumPlanets; ++j)
    {
      r[0] = pos[0][i] - m_plPos[0][j];
      r[1] = pos[1][i] - m_plPos[1][j];
      myfloat dist = sqrt(r[0] * r[0] + r[1] * r[1]),
              grav = m_plMass[j] / (dist * dist * dist);

      if (j == 0)
        dist0 = dist;

      acc[0][i] += grav * r[0];
      acc[1][i] += grav * r[1];

      if (dist < m_plRad[j]*10)
      {
        std::cout << "Tracer collided with " << m_plName[j] << "\n";
        ResetTracer(i);
        break;
      }
    }

    acc[0][i] *= -GAMMA_;
    acc[1][i] *= -GAMMA_;

    // Integrate (Verlet)
    buf[0] = pos[0][i];
    buf[1] = pos[1][i];
    pos[0][i] = 2 * pos[0][i] - pos_old[0][i] + acc[0][i] * m_dt * m_dt;
    pos[1][i] = 2 * pos[1][i] - pos_old[1][i] + acc[1][i] * m_dt * m_dt;
    pos_old[0][i] = buf[0];
    pos_old[1][i] = buf[1];

    // excentricity data
    m_trMinDist[i] = std::min(dist0, m_trMinDist[i]);
    m_trMaxDist[i] = std::max(dist0, m_trMaxDist[i]);
    if (pos[0][i] > m_plPos[0][0] && pos_old[0][i] < m_plPosOld[0][0] && pos[1][i] > m_plPos[1][0])
    {
      m_trExc[i] = m_trMinDist[i] / m_trMaxDist[i];
      m_trMinDist[i] = 9e29;
      m_trMaxDist[i] = 0;

      if (m_trExc[i]<0.75)
      {
        ResetTracer(i);
        continue;
      }
 
    }

    if (pos[0][i] < (m_plPos[0][0] + m_xmin) ||
            pos[0][i] > (m_plPos[0][0] + m_xmax) ||
            pos[1][i] < (m_plPos[1][0] + m_ymin) ||
            pos[1][i] > (m_plPos[1][0] + m_ymax))
      ResetTracer(i);
  }
}

//----------------------------------------------------------------------------------------

int PlanetarySystem::GetPlanetDataArrays(myfloat **x,
                                         myfloat **y,
                                         myfloat **rad,
                                         myfloat **mass,
                                         int **col)
{
  x && (*x = &m_plPos[0][0]);
  y && (*y = &m_plPos[1][0]);
  col && (*col = &m_plColor[0]);
  rad && (*rad = &m_plRad[0]);
  mass && (*mass = &m_plMass[0]);

  return m_nNumPlanets;
}

//----------------------------------------------------------------------------------------

int PlanetarySystem::GetTracerDaraArrays(myfloat **x, myfloat **y, myfloat **vx, myfloat **vy)
{
  x && (*x = &m_trPos[0][0]);
  y && (*y = &m_trPos[1][0]);
  return m_trPos[0].size();
}

//----------------------------------------------------------------------------------------

int PlanetarySystem::GetStatData(myfloat **exc)
{
  exc && (*exc = &m_trExc[0]);
  return m_trExc.size();
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::SetTracerTo(myfloat x, myfloat y)
{
  static std::size_t idx = 0;

  ++idx;

  if (idx >= m_trPos[0].size())
    idx = 0;

  myfloat vel[2];

  GetOrbitalVelocity(m_center, x, y, vel[0], vel[1]);
  m_trPos[0][idx] = x;
  m_trPos[1][idx] = y;
  m_trPosOld[0][idx] = m_trPos[0][idx] - m_dt * vel[0];
  m_trPosOld[1][idx] = m_trPos[1][idx] - m_dt * vel[1];
}

//----------------------------------------------------------------------------------------

void PlanetarySystem::GetOrbitalVelocity(int idx_main, myfloat x, myfloat y, myfloat &vx, myfloat &vy)
{
  // Calculate distance from the planet with index idx_main
  myfloat r[2], dist;
  r[0] = x - m_plPos[0][idx_main];
  r[1] = y - m_plPos[1][idx_main];
  dist = sqrt(r[0] * r[0] + r[1] * r[1]);

  // Based on the distance from the sun calculate the velocity needed to maintain a circular orbit
  myfloat v = sqrt(GAMMA_ * m_plMass[idx_main] / dist);

  // Calculate a suitable vector perpendicular to r for the velocity of the tracer
  vx = (r[1] / dist) * v;
  vy = (-r[0] / dist) * v;
}
