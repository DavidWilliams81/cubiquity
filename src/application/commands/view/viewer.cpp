#include "viewer.h"

#include "base/logging.h"

#include "cubiquity.h"
#include "visibility.h"
#include "raytracing.h"
#include "utility.h"

#include <filesystem>
#include <fstream>
#include <vector>

using namespace Cubiquity;

Viewer::Viewer(const std::string& filename, WindowType windowType)
	: Window(windowType)
{
	log_info("Opening volume '{}'", filename);
	if (!mVolume.load(filename))
	{
		log_error("Failed to open volume!");
		exit(EXIT_FAILURE);
	}
	log_info("Done");

	// FIXME - This cube doesn't pathtrace properly
	/*for (int z = -8; z < 7; z++)
	{
		for (int y = -8; y < 7; y++)
		{
			for (int x = -8; x < 7; x++)
			{
				mVolume.setVoxel(x, y, z, 250);
			}
		}
	}*/

	Metadata metadata = loadMetadataForVolume(filename);

	// Build an array of colours from the material data for uploading to the GPU.
	std::fill(begin(mColours), end(mColours), Metadata::Warning.diffuse);
	for (int i = 0; i < metadata.materials.size(); i++) {
		mColours[i] = metadata.materials[i].diffuse;
	}
}

void Viewer::onInitialise()
{
	onVolumeModified();

	mVolume.setTrackEdits(true);

	//Box3i bounds = computeBounds(mVolume, [](MaterialId matId) { return matId != 0; });
	Timer timer;

	uint8 outside_material;
	int32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&mVolume, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);

	Vector3d lower({ static_cast<float>(lower_x), static_cast<float>(lower_y), static_cast<float>(lower_z) });
	Vector3d upper({ static_cast<float>(upper_x), static_cast<float>(upper_y), static_cast<float>(upper_z) });
	log_info("Lower bound = ({},{},{})", lower.x(), lower.y(), lower.z());
	log_info("Upper bound = ({},{},{})", upper.x(), upper.y(), upper.z());
	log_info("Bounds estimation took {} seconds", timer.elapsedTimeInSeconds());

	Vector3d centre = (lower + upper) * 0.5;

	if (outside_material == 0) // Solid object, point camera at centre and move it back
	{
		double halfDiagonal = length(upper - lower) * 0.5;

		// Centred along x, then back and up a bit
		mCamera.position = Vector3d({ centre.x(), centre.y() - halfDiagonal, centre.z() + halfDiagonal });

		// Look down 45 degrees
		mCamera.pitch = -(Pi / 4.0f);
		mCamera.yaw = 0.0f;
	}
	else // Hollow object, place camera at centre.
	{
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

	mFrameNumber++;
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
		Ray3f ray = static_cast<Ray3f>(mCamera.rayFromViewportPos(event.x, event.y, width(), height()));
		SubDAGArray subDAGs = findSubDAGs(
			Internals::getNodes(volume()).nodes(), getRootNodeIndex(volume()));
		RayVolumeIntersection intersection = intersectVolume(mVolume, subDAGs, ray, false);
		if (intersection.hit)
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
