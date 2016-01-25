/* 
 * File:   planetrender.hpp
 * Author: user
 *
 * Created on 1. September 2008, 17:45
 */

#ifndef _PLANETRENDER_HPP
#define	_PLANETRENDER_HPP

//--- Standard includes ------------------------------------------------------
#include <vector>

//--- SDL includes -----------------------------------------------------------
#include <SDL/SDL.h>
#include <SDL/SDL_gfxPrimitives.h>

//--- Simulation -------------------------------------------------------------
#include "simtypes.hpp"
#include "planetarysystem.hpp"


class PlanetaryRenderer
{
public:
  
  PlanetaryRenderer();
  
  void InitSDL(int width, int height, PlanetarySystem *pPlanets);
  void Render();
  void PollMessages();
  int GetFrameTime() const;
  
private:

  void DrawPlanets();
  void DrawSolarSystem();
  void DrawStat();
  void DrawResonancemarker(std::size_t idx);
  void DrawTracer();

  SDL_Surface *m_pScreen;
  SDL_Event m_event;

  PlanetarySystem *m_pPlanets;
  int m_width;
  int m_height;
  int m_xcenter;
  int m_ycenter;
  myfloat m_fov; // the field of view
  myfloat m_scale; // The scaling factor to squeeze the fov into the window
  int m_center;
  int m_frameTime;
  int m_resonance;
  bool m_bShowResonance;
  bool m_bContSource;
  bool m_bShowSolarSystem;
};

#endif	/* _PLANETRENDER_HPP */

