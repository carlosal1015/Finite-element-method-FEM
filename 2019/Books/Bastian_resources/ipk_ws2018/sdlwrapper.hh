#ifndef SDLCANVAS_HH
#define SDLCANVAS_HH

#include <array>
#include <string>
#include <iostream>

#include "SDL2/SDL.h"

/*
 * SDL internals wrapped in a modern C++ class
 */
class SDLCanvas
{
  // counts number of open windows
  static int numCanvas;

  // SDL stuff
  SDL_Window*   window;
  SDL_Renderer* renderer;
  SDL_Texture*  texture;
  SDL_Event     event;

  public:

  // create a window with given name of given width and height
  SDLCanvas(const std::string& name, int winWidth, int winHeight)
  {
    // initialize if we are the first
    if (numCanvas == 0)
    {
      if (SDL_WasInit(SDL_INIT_VIDEO) > 0)
        throw std::runtime_error("SDL already initialized, only one instance allowed");

      if (SDL_Init(SDL_INIT_VIDEO) < 0)
        throw std::runtime_error("Failed to initialize SDL");
    }

    // window == display
    window = SDL_CreateWindow(name.c_str(), SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED, winWidth, winHeight, 0);
    if (not window)
      throw std::runtime_error("Failed to create SDL window");

    // renderer == interface between canvas and display
    renderer = SDL_CreateRenderer(window, -1, 0);
    if (not renderer)
      throw std::runtime_error("Failed to create SDL renderer");

    // texture == canvas
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, winWidth, winHeight);
    if (not texture)
      throw std::runtime_error("Failed to create SDL texture");

    SDL_SetRenderTarget(renderer, texture);

    // one more window
    numCanvas++;
  }

  ~SDLCanvas()
  {
    SDL_DestroyWindow(window);

    // one less window
    numCanvas--;

    // clean up if we are the last
    if (numCanvas == 0)
      SDL_Quit();
  }

  // set color used for subsequent pixels
  void setDrawColor(const std::array<int,3>& color)
  {
    SDL_SetRenderDrawColor(renderer, color[0], color[1], color[2], 0x00);
  }

  // draw a pixel with current color
  void drawPixel(int x, int y)
  {
    SDL_RenderDrawPoint(renderer, x, y);
  }

  // explicitly choose color for pixel
  void drawPixel(int x, int y, std::array<int,3> color)
  {
    SDL_SetRenderDrawColor(renderer, color[0], color[1], color[2], 0x00);
    SDL_RenderDrawPoint(renderer, x, y);
  }

  // fill texture with current color
  void clear()
  {
    SDL_RenderClear(renderer);
  }

  // explicitly choose color for pixel
  void clear(const std::array<int,3>& color)
  {
    setDrawColor(color);
    SDL_RenderClear(renderer);
  }

  // use state of canvas to update displayed image
  void display()
  {
    SDL_SetRenderTarget(renderer, NULL);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
    SDL_SetRenderTarget(renderer, texture);
  }

  // check whether escape was pressed
  bool windowClosed()
  {
    const Uint8* keystates = SDL_GetKeyboardState(NULL);
    if (keystates[SDL_SCANCODE_ESCAPE])
      return true;

    SDL_Event event;
    SDL_PollEvent(&event);
    if(event.type == SDL_QUIT)
      return true;

    return false;
  }
};

// initially there are no windows
int SDLCanvas::numCanvas = 0;

#endif // SDLCANVAS_HH
