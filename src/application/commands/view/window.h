#ifndef CUBIQUITY_WINDOW_H
#define CUBIQUITY_WINDOW_H

// FIXME - We should not need to include logging at this point, but there seems
// to be an issue that including SDL.h can somehow break logging if SDL.h is
// included first. It is probably related to libfmt  (possibly due to header
// only mode). It may be related to <algorithm> because including <algorithm>
// prior to SDL.h is an alternative solution to this problem. Perhaps SDL
// suppresses certain later includes?
// This only seems to occur when cross-compiling to Windows with MinGW. In the
// future I hope to remove the libfmt dependancy (once std::print and
// std::format are better supported) so I won't worry about it too much for now.
#include "base/logging.h"

// Useful tips on setting up and using SDL:
// 	https://nullprogram.com/blog/2023/01/08/
// 	https://stackoverflow.com/questions/64396979/how-
// 		do-i-use-sdl2-in-my-programs-correctly
#include "SDL.h"

enum class KeyState { Up, Down };
enum class MouseButtonState { Up, Down };
enum class WindowType { Software, OpenGL };

class Window
{
public:

	// Constuctor
	Window(WindowType windowType)
		: mWindowType(windowType) {}

	// Main interface
	void show(int width, int height);
	void close();

	// Dimensions
	int width();
	int height();

	// Input state checking
	KeyState keyState(SDL_Scancode scancode);
	MouseButtonState mouseButtonState(Uint32 button);

	// For software rendering
	SDL_Surface* surface();

protected:

	// Event handling
	virtual void onInitialise() {}
	virtual void onUpdate(float deltaTime) {}
	virtual void onShutdown() {}

	virtual void onKeyDown(const SDL_KeyboardEvent& event) {}
	virtual void onKeyUp(const SDL_KeyboardEvent& event) {}
	virtual void onMouseButtonDown(const SDL_MouseButtonEvent& event) {}
	virtual void onMouseButtonUp(const SDL_MouseButtonEvent& event) {}
	virtual void onMouseMotion(const SDL_MouseMotionEvent& event) {}
	virtual void onWindowSizeChanged(int width, int height) {}

private:

	bool mRunning = true;
	SDL_Window* mSDLWindow = nullptr;
	WindowType mWindowType;
};

#endif // CUBIQUITY_WINDOW_H
