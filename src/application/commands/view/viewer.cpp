#include "viewer.h"

#include "rendering.h"
#include "utility.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace Cubiquity;

Viewer::Viewer(const std::string& filename, WindowType windowType)
	: Window(windowType)
{
	std::cout << "Opening volume \'" << filename << "\'... ";
	if (!mVolume.load(filename))
	{
		std::cout << " failed to open volume!" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << " done" << std::endl;

	mMaterials.load(getMaterialsPath(filename));
}

void Viewer::onInitialise()
{
	mVolume.setTrackEdits(true);

	//Box3i bounds = computeBounds(mVolume, [](MaterialId matId) { return matId != 0; });
	auto result = estimateBounds(mVolume);
	MaterialId outsideMaterialId = result.first;
	Box3i bounds = result.second;

	if (outsideMaterialId == 0) // Solid object, point camera at centre and move it back
	{
		Vector3d lower = static_cast<Vector3d>(bounds.lower());
		Vector3d centre = static_cast<Vector3d>(bounds.centre());
		Vector3d upper = static_cast<Vector3d>(bounds.upper());
		double halfDiagonal = length(upper - lower) * 0.5;

		// Centred along x, then back and up a bit
		mCamera.position = Vector3d({ centre.x(), centre.y() - halfDiagonal, centre.z() + halfDiagonal });

		// Look down 45 degrees
		mCamera.pitch = -(Pi / 4.0f);
		mCamera.yaw = 0.0f;
	}
	else // Hollow object, place camera at centre.
	{
		Vector3d centre = static_cast<Vector3d>(bounds.centre());
		centre += Vector3d({ 0.1, 0.1, 0.1 }); // Hack to help not be on a certain boundary which triggers assert in debug mode.
		mCamera.position = Vector3d({ centre.x(), centre.y(), centre.z() });

		// Look straight ahead
		mCamera.pitch = 0.0f;
		mCamera.yaw = 0.0f;
	}
}

void Viewer::onUpdate(float deltaTime)
{	
	// Hold shift to move fast
	float speedMultiplier = 1.0f;
	if (keyState(SDL_SCANCODE_LCTRL) == KeyState::Down) speedMultiplier = 0.1f;
	if (keyState(SDL_SCANCODE_LSHIFT) == KeyState::Down) speedMultiplier = 10.0f;

	// Move forward
	if (keyState(SDL_SCANCODE_W) == KeyState::Down)
	{
		mCamera.position += mCamera.forward() * static_cast<double>(deltaTime * speedMultiplier * CameraMoveSpeed);
		onCameraModified();
	}
	// Move backward
	if (keyState(SDL_SCANCODE_S) == KeyState::Down)
	{
		mCamera.position -= mCamera.forward() * static_cast<double>(deltaTime * speedMultiplier * CameraMoveSpeed);
		onCameraModified();
	}
	// Strafe right
	if (keyState(SDL_SCANCODE_D) == KeyState::Down)
	{
		mCamera.position += mCamera.right() * static_cast<double>(deltaTime * speedMultiplier * CameraMoveSpeed);
		onCameraModified();
	}
	// Strafe left
	if (keyState(SDL_SCANCODE_A) == KeyState::Down)
	{
		mCamera.position -= mCamera.right() * static_cast<double>(deltaTime * speedMultiplier * CameraMoveSpeed);
		onCameraModified();
	}
}

void Viewer::onKeyUp(const SDL_KeyboardEvent& event)
{
	if (event.keysym.sym == SDLK_ESCAPE)
	{
		close();
	}

	if ((event.keysym.mod & KMOD_CTRL) && event.keysym.sym == SDLK_z)
	{
		event.keysym.mod& KMOD_SHIFT ? mVolume.redo() : mVolume.undo();
		onVolumeModified();
	}
}

void Viewer::onMouseMotion(const SDL_MouseMotionEvent& event)
{
	if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		// Compute new orientation. Window origin is top-left (not bottom-
		// left) So we subtract the relative y motion instead of adding it.
		mCamera.yaw += CameraTurnSpeed * event.xrel;
		mCamera.pitch -= CameraTurnSpeed * event.yrel;
		onCameraModified();
	}
}

void Viewer::onMouseButtonDown(const SDL_MouseButtonEvent& event)
{
	if (event.button == SDL_BUTTON(SDL_BUTTON_LEFT))
	//if (mouseButtonState(SDL_BUTTON_LEFT) == MouseButtonState::Down)
	{
		//std::cout << "x = " << event.x << ", y = " << event.y << std::endl;

		Ray3d ray = mCamera.rayFromViewportPos(event.x, event.y, width(), height());

		RayVolumeIntersection intersection = intersectVolume(mVolume, ray);
		if (intersection)
		{
			SphereBrush brush(static_cast<Vector3f>(intersection.position), 30);
			mVolume.fillBrush(brush, 0);

			onVolumeModified();
		}
	}

	if(event.button == SDL_BUTTON(SDL_BUTTON_RIGHT))
	//if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		SDL_SetRelativeMouseMode(SDL_TRUE);
	}
}

void Viewer::onMouseButtonUp(const SDL_MouseButtonEvent& event)
{
	if (event.button == SDL_BUTTON(SDL_BUTTON_RIGHT))
	//if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		SDL_SetRelativeMouseMode(SDL_FALSE);
	}
}
