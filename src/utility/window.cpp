#include "window.h"

// Include GLAD
#include <glad/glad.h>

// Include standard headers
#include <iostream>

void Window::show(int width, int height)
{
	SDL_Init(SDL_INIT_VIDEO);
	Uint32 flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIDDEN;

	// Configure for OpenGL
	if (mWindowType == WindowType::OpenGL)
	{
		flags |= SDL_WINDOW_OPENGL;
		SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
		SDL_GL_SetAttribute(SDL_GL_ACCELERATED_VISUAL, 1);
		SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
		SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
		SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
		SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	}

	// Create the window
	mSDLWindow = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, flags);

	// More OpenGL configuration
	SDL_GLContext context = nullptr;
	if (mWindowType == WindowType::OpenGL)
	{
		context = SDL_GL_CreateContext(mSDLWindow);
		if (!gladLoadGLLoader((GLADloadproc)SDL_GL_GetProcAddress))
		{
			std::cerr << "Failed to initialize the OpenGL context." << std::endl;
			return;
		}
	}

	// Perform user-provided initialisation.
	onInitialise();

	// Now we are ready to go. Show the window and give it focus (the
	// terminal may have had focus if user was entering something there).
	SDL_ShowWindow(mSDLWindow);
	SDL_RaiseWindow(mSDLWindow);

	// Some params are better initialised on size change (camera, fullscreen render targets, etc).
	// It would be nice if SDL generated a size changed event automatically when the window was created
	// or shown, but perhaps OpenGL wouldn't have been initilised by then anyway. Do it manually.
	onWindowSizeChanged(width, height);

	Uint32 previousUpdateStartTimeMs = SDL_GetTicks();

	while(mRunning)
	{
		SDL_Event event;
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
			case SDL_QUIT:
				close();
				break;
			case SDL_KEYDOWN:
				onKeyDown(event.key);
				break;
			case SDL_KEYUP:
				onKeyUp(event.key);
				break;
			case SDL_MOUSEBUTTONDOWN:
				onMouseButtonDown(event.button);
				break;
			case SDL_MOUSEBUTTONUP:
				onMouseButtonUp(event.button);
				break;
			case SDL_MOUSEMOTION:
				onMouseMotion(event.motion);
				break;
			case SDL_WINDOWEVENT:
				if (event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
				{
					onWindowSizeChanged(event.window.data1, event.window.data2);
				}
				break;
			}
		}

		// Compute time difference between current and last update
		Uint32 updateStartTimeMs = SDL_GetTicks();
		float deltaTimeSeconds = static_cast<float>(updateStartTimeMs - previousUpdateStartTimeMs) / 1000.0f;
		previousUpdateStartTimeMs = updateStartTimeMs;

		// Run user-provided update code
		onUpdate(deltaTimeSeconds);

		// Update the window in the appropriate way.
		switch (mWindowType)
		{
		case WindowType::Software:
			SDL_UpdateWindowSurface(mSDLWindow);
			break;
		case WindowType::OpenGL:
			SDL_GL_SwapWindow(mSDLWindow);
			break;
		}
	}

	// User-provided shutdown handler
	onShutdown();

	if (mWindowType == WindowType::OpenGL)
	{
		SDL_GL_DeleteContext(context);
	}

	SDL_DestroyWindow(mSDLWindow);
	SDL_Quit();
}

void Window::close()
{
	mRunning = false;
}

KeyState Window::keyState(SDL_Scancode scancode)
{
	const Uint8* state = SDL_GetKeyboardState(NULL);
	return state[scancode] ? KeyState::Down : KeyState::Up;
}

MouseButtonState Window::mouseButtonState(Uint32 button)
{
	const Uint32 state = SDL_GetMouseState(NULL, NULL);
	return state & SDL_BUTTON(button) ? MouseButtonState::Down : MouseButtonState::Up;
}

SDL_Surface* Window::surface()
{
	return SDL_GetWindowSurface(mSDLWindow);
}
