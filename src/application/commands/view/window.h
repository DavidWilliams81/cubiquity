#ifndef CUBIQUITY_WINDOW_H
#define CUBIQUITY_WINDOW_H

// Include SDL
#define SDL_MAIN_HANDLED
#include <SDL.h>

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
