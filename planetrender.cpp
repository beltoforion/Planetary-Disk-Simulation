#include "planetrender.hpp"

#include <cassert>
#include <stdexcept>
#include <sstream>


//----------------------------------------------------------------------------

PlanetaryRenderer::PlanetaryRenderer()
: m_pScreen(NULL)
, m_pPlanets(NULL)
, m_width(0)
, m_height(0)
, m_xcenter(0)
, m_ycenter(0)
, m_fov(0)
, m_scale(0)
, m_center(1)
, m_frameTime(0)
, m_resonance(0)
, m_bShowResonance(true)
, m_bContSource(false)
, m_bShowSolarSystem(false)
{}

//----------------------------------------------------------------------------

void PlanetaryRenderer::InitSDL(int width, int height, PlanetarySystem *pPlanets)
{
  m_pPlanets = pPlanets;

  if (SDL_Init(SDL_INIT_VIDEO) == -1)
    throw std::runtime_error(SDL_GetError());

  atexit(SDL_Quit);
  m_pScreen = SDL_SetVideoMode(width, height, 16, SDL_HWSURFACE);
  if (!m_pScreen)
    throw std::runtime_error(SDL_GetError());

  m_width = width;
  m_height = height;
  m_xcenter = width >> 1;
  m_ycenter = height >> 1;

  m_fov = 2000000000;
  m_scale = m_width / m_fov;
  m_bShowResonance = true;
  m_center = 0;
}

//----------------------------------------------------------------------------

void PlanetaryRenderer::Render()
{
  SDL_FillRect(m_pScreen, NULL, 0);

  DrawPlanets();

  std::size_t num = m_pPlanets->GetNumPlanets();
  for (std::size_t i=0; i<num; ++i)
  {
    if (m_resonance & (1 << i) )
      DrawResonancemarker(i);
  }
  
  if (m_bShowSolarSystem)
    DrawSolarSystem();
  
  DrawTracer();
  DrawStat();
  
  SDL_UpdateRect(m_pScreen, 0, 0, m_width, m_height);
}

//----------------------------------------------------------------------------

void PlanetaryRenderer::DrawStat()
{
  char szBuf[1024];

  sprintf(szBuf, "Number of tracers: %d", (int)m_pPlanets->GetNumActive());
  stringColor(m_pScreen, 20, 20, szBuf, 0xffffffff);

  sprintf(szBuf, "Total time: %2.2f yr", (double)(m_pPlanets->GetTime()/(86400*365.25)));
  stringColor(m_pScreen, 20, 35, szBuf, 0xffffffff);

  sprintf(szBuf, "Steps per second: %d", (int)m_pPlanets->GetFPS());
  stringColor(m_pScreen, 20, 50, szBuf, 0xffffffff);
  
  if (m_pPlanets->GetFPS())
  {
    double total_per_sec = m_pPlanets->GetFPS() * m_pPlanets->GetNumActive();
    sprintf(szBuf, "Time per tracer: %2.2g us", (double)1000000 / total_per_sec);
    stringColor(m_pScreen, 20, 65, szBuf, 0xffffffff);
  }

  sprintf(szBuf, "Frame delay: %d", m_frameTime);
  stringColor(m_pScreen, 20, 80, szBuf, 0xffffffff);
}

//----------------------------------------------------------------------------
void PlanetaryRenderer::PollMessages()
{
  static int x(0), y(0);

  if (m_bContSource)
    m_pPlanets->SetTracerTo(x, y);

  while (SDL_PollEvent(&m_event))
  {
    switch (m_event.type)
    {
      case SDL_QUIT:
        throw std::runtime_error("Quitting the program");
        break;

      case SDL_MOUSEBUTTONDOWN:
        m_pPlanets->SetTracerTo((m_event.button.x - m_xcenter) * m_fov / m_width,
                                (m_event.button.y - m_ycenter) * m_fov / m_height);
        x = (m_event.button.x - m_xcenter) * m_fov / m_width;
        y = (m_event.button.y - m_ycenter) * m_fov / m_height;
        break;

      case SDL_KEYDOWN:
        switch (m_event.key.keysym.sym)
        {
          case SDLK_PLUS:
            m_fov *= 1.2;
            break;
          case SDLK_MINUS:
            m_fov *= 0.8;
            break;
          case SDLK_s:
            m_bShowSolarSystem ^= true;
            break;
          case SDLK_r:
            m_bShowResonance ^= true;
            break;
          case SDLK_c:
            m_bContSource ^= true;
            break;
          case SDLK_0:
            if (m_pPlanets->GetNumPlanets()>=5)
              m_center = 4;
            break;
          case SDLK_1:
            if (m_pPlanets->GetNumPlanets()>=1)
              m_center = 0;
            break;
          case SDLK_2:
            if (m_pPlanets->GetNumPlanets()>=2)
              m_center = 1;
            break;
          case SDLK_3:
            if (m_pPlanets->GetNumPlanets()>=3)
              m_center = 2;
            break;
          case SDLK_4:
            if (m_pPlanets->GetNumPlanets()>=4)
              m_center = 3;
            break;
          case SDLK_F1:
            m_resonance ^= 1<<0;
            break;
          case SDLK_F2:
            m_resonance ^= 1<<1;
            break;
          case SDLK_F3:
            m_resonance ^= 1<<2;
            break;
          case SDLK_F4:
            m_resonance ^= 1<<3;
            break;
          case SDLK_KP_MINUS:
            m_frameTime -= 5000;
            m_frameTime = std::max(0, m_frameTime);
            break;
          case SDLK_KP_PLUS:
            m_frameTime += 5000;
            break;
          default:
            break;
        }
        m_scale = m_width / m_fov;
        break;

      default:
        break;
    } // switch event type
  }
}

//----------------------------------------------------------------------------

int PlanetaryRenderer::GetFrameTime() const
{
  return m_frameTime;
}

//----------------------------------------------------------------------------

void PlanetaryRenderer::DrawPlanets()
{
  assert(m_pPlanets);

  myfloat *xpos, *ypos, *rad;
  int *col;

  lineColor(m_pScreen, m_xcenter, 0, m_xcenter, m_height, 0x4444aaff);
  lineColor(m_pScreen, 0, m_ycenter, m_width, m_ycenter, 0x4444aaff);

  // Draw planets
  int np = m_pPlanets->GetPlanetDataArrays(&xpos, &ypos, &rad, NULL, &col);
  for (std::size_t i = 0; i < (std::size_t)np; ++i)
  {
    // Calculate on screen location
    int x = m_xcenter + (xpos[i] - xpos[m_center]) * m_scale;
    int y = m_ycenter + (ypos[i] - ypos[m_center]) * m_scale;

    if (x < 0 || x > m_width || y < 0 || y > m_height)
      continue;
    
    // all the other planets
    lineColor(m_pScreen, x, y - 15, x, y + 15, col[i]);
    lineColor(m_pScreen, x - 15, y, x + 15, y, col[i]);
    filledCircleColor(m_pScreen, x, y, rad[i] * m_scale, col[i]);
  }
}

//----------------------------------------------------------------------------

void PlanetaryRenderer::DrawSolarSystem()
{
  // Draw a simplified version of the solar system (circular orbits)  
  // This is usefull for estimating distances ond orbits
  const int col[] = { 0xffffffff,    // merkur
                      0xfbe66fff,    // Venus
                      0x007b62ff,    // Earth,
                      0xf2622aff,    // Mars
                      0xff6123ff,    // Jupiter
                      0xffd800ff,    // Saturn
                      0x6590ffff,    // Uranus
                      0x0536b2ff };  // Neptune
  const myfloat dist[] = { 57910000.0, 
                           108210000.0,
                           149600000.0,
                           227920000.0,
                           778570000.0,
                           1433530000.0,
                           2872460000.0, 
                           4495060000.0 };
  const char *names[] = {"Mercury", 
                         "Venus", 
                         "Earth", 
                         "Mars", 
                         "Jupiter",
                         "Saturn",
                         "Uranus", 
                         "Neptune"};                         
  myfloat *xpos, *ypos;
  int np = m_pPlanets->GetPlanetDataArrays(&xpos, &ypos, NULL, NULL, NULL);
  std::size_t sz = sizeof(dist)/sizeof(myfloat);
  for (std::size_t i=0; i<sz; ++i)
  {
    circleColor(m_pScreen, 
                m_xcenter + (xpos[0] - xpos[m_center]) * m_scale, 
                m_ycenter + (ypos[0] - ypos[m_center]) * m_scale, 
                dist[i] * m_scale, 
                col[i]);

    stringColor(m_pScreen, 
                m_xcenter + (xpos[0] - xpos[m_center]) * m_scale, 
                m_ycenter + (ypos[0] - ypos[m_center]) * m_scale + dist[i] * m_scale - 10, 
                names[i],
                col[i]);
  }
}

//----------------------------------------------------------------------------

void PlanetaryRenderer::DrawResonancemarker(std::size_t idx)
{
  assert(m_pPlanets);

  myfloat *xpos, *ypos, *mass;
  int *col, idx_main = 0;

  m_pPlanets->GetPlanetDataArrays(&xpos, &ypos, NULL, &mass, &col);
  myfloat r[2];
  r[0] = xpos[idx_main] - xpos[idx];
  r[1] = ypos[idx_main] - ypos[idx];

  // How far is the planet from the sun?
  myfloat dist = sqrt(r[0] * r[0] + r[1] * r[1]);

  // What is the absolute velocity a planet this far from the sun would need for a circular orbit?
  myfloat vx, vy;
  m_pPlanets->GetOrbitalVelocity(idx_main, xpos[idx_main] + dist, ypos[idx_main], vx, vy);

  // How many time would it take for a planet with this velocity to orbit sun? 
  myfloat period = 2 * M_PI * dist / fabs(vy); // in seconds
  myfloat periods_res[6];
  periods_res[0] = period;
  periods_res[1] = period * 1.0 / 2.0;
  periods_res[2] = period * 1.0 / 3.0;
  periods_res[3] = period * 1.0 / 4.0;
  periods_res[4] = period * 2.0 / 3.0;
  periods_res[5] = period * 2.0 / 5.0;

  const char *names[] = { " 1:1", " 2:1", " 3:1", " 4:1", " 3:2", " 5:2"};
  // Find out what orbits are needed to have these sppeds and draw them in the window
  for (std::size_t i = 0; i < sizeof (periods_res) / sizeof (myfloat); ++i)
  {
    myfloat rad = pow((6.67428e-11 * mass[idx_main] * periods_res[i] * periods_res[i]) / (4 * M_PI * M_PI), 1.0 / 3);

    stringColor(m_pScreen, 
                m_xcenter + (xpos[idx_main] - xpos[m_center]) * m_scale,
                m_ycenter + (ypos[idx_main] - ypos[m_center]) * m_scale - rad / 1000 * m_scale - 10,
                names[i],
                col[idx]);
    
    circleColor(m_pScreen,
                m_xcenter + (xpos[idx_main] - xpos[m_center]) * m_scale,
                m_ycenter + (ypos[idx_main] - ypos[m_center]) * m_scale,
                rad / 1000 * m_scale,
                col[idx]);
  }
}

//----------------------------------------------------------------------------

void PlanetaryRenderer::DrawTracer()
{
  assert(m_pPlanets);
  assert(m_pScreen);

  myfloat *xpos, *ypos,
          *pxpos, *pypos,
          *exc;

  m_pPlanets->GetPlanetDataArrays(&pxpos, &pypos, NULL, NULL, NULL);
  m_pPlanets->GetTracerDaraArrays(&xpos, &ypos, NULL, NULL);
  m_pPlanets->GetStatData(&exc);
  
  for (std::size_t i = 0; i < m_pPlanets->GetNumActive(); ++i)
  {
    int x = m_xcenter + (xpos[i] - pxpos[m_center]) * m_scale;
    int y = m_ycenter + (ypos[i] - pypos[m_center]) * m_scale;

    if (x > 0 && x < m_width && y > 0 && y < m_height)
    {
      //int color = std::min((myfloat)((1-exc[i])*10000.0), (myfloat)255.0);
      pixelColor(m_pScreen, x, y, 0xffffffff);
      //pixelRGBA(m_pScreen, x, y, color, 255-color, 0, 255);
    }
  }
}
