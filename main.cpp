//--- Standard includes ------------------------------------------------------
#include <stdexcept>
#include <iostream>
#include <string>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include <fstream>
#include <sys/time.h>

//--- SDL includes -----------------------------------------------------------
#include <SDL/SDL.h>
#include <SDL/SDL_gfxPrimitives.h>

//--- Simulation -------------------------------------------------------------
#include "planetarysystem.hpp"
#include "planetrender.hpp"
#include "simtypes.hpp"


PlanetaryRenderer s_pRenderer;
PlanetarySystem s_planets;
pthread_mutex_t s_mtxCalc;
volatile bool s_bRunning = true;

//----------------------------------------------------------------------------

void *DoCalculation(void *threadID)
{
  s_planets.InitTracerRandom(30000, 1e9);
  s_planets.DoTracerReset(true);

  while (s_bRunning)
  {
    s_planets.Move();
  }

  
  pthread_exit(NULL);
}

//----------------------------------------------------------------------------

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cout << "usage: " << argv[0] << "config.cfg\n";
    return (EXIT_FAILURE);
  }

  s_planets.InitFromFile(argv[1]);  

  pthread_t thread;
  pthread_attr_t attr;

  // initialize threads
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_create(&thread, &attr, DoCalculation, (void *) & thread);

  try
  {
    s_pRenderer.InitSDL(800, 800, &s_planets);

    for (;;)
    {
      usleep(s_pRenderer.GetFrameTime());
      s_pRenderer.Render();
      s_pRenderer.PollMessages();
    }
  }

  catch(std::exception & exc)
  {
    std::cout << exc.what();
    s_bRunning = false;
  }

  // wait for thread termination
  pthread_join(thread, NULL);
  pthread_attr_destroy(&attr);

  return (EXIT_SUCCESS);
}

