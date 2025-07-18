#include "viewer.h"

#include "base/logging.h"
#include "base/serialize.h"

#include "cubiquity.h"
#include "extraction.h"
#include "raytracing.h"
#include "utility.h"

#include <filesystem>
#include <fstream>
#include <vector>

Viewer::Viewer(const std::string& volume_path, WindowType windowType)
	: Window(windowType)
{
	// FIXME - Change this to use the loadVolume()
	// function once we have that returning by value
	log_info("Opening volume '{}'", volume_path);
	if (!mVolume.load(volume_path))
	{
		log_error("Failed to open volume!");
		exit(EXIT_FAILURE);
	}

	mMetadata.load(getMetadataPath(volume_path));
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

	// Build an array of colours from the material data for uploading to the GPU.
	vec3 purple = { 1.0f, 0.0f, 1.0f };
	std::fill(begin(mColours), end(mColours),purple ); // Purple
	for (int i = 0; i < mMetadata.material_count(); i++) {
		mColours[i] = mMetadata.find_material_base_color(i);
	}
}

void Viewer::onInitialise()
{
	onVolumeModified();

	mVolume.setTrackEdits(true);

	Cubiquity::Timer timer;

	ivec3 lower = mMetadata.find_lower_bound();
	ivec3 upper = mMetadata.find_upper_bound();
	dvec3 centre = (lower + upper) * 0.5;

	u8 outside_material = mVolume.voxel(I32_MAX, I32_MAX, I32_MAX); // Corner

	if (outside_material == 0) // Solid object, point camera at centre and move it back
	{
		double halfDiagonal = length(upper - lower) * 0.5;

		// Centred along x, then back and up a bit
		mCamera.position = dvec3({ centre.x, centre.y - halfDiagonal, centre.z + halfDiagonal });

		// Look down 45 degrees
		mCamera.pitch = -(Pi / 4.0f);
		mCamera.yaw = 0.0f;
	}
	else // Hollow object, place camera at centre.
	{
		centre += dvec3({ 0.1, 0.1, 0.1 }); // Hack to help not be on a certain boundary which triggers assert in debug mode.
		mCamera.position = dvec3({ centre.x, centre.y, centre.z });

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
		Cubiquity::SubDAGArray subDAGs = Cubiquity::findSubDAGs(
			Cubiquity::Internals::getNodes(volume()).nodes(), Cubiquity::getRootNodeIndex(volume()));
		Cubiquity::RayVolumeIntersection intersection = intersectVolume(mVolume, subDAGs,
			ray.mOrigin.x, ray.mOrigin.y, ray.mOrigin.z,
			ray.mDir.x, ray.mDir.y, ray.mDir.z,
			false);
		if (intersection.hit)
		{
			Cubiquity::SphereBrush brush(intersection.position.x, intersection.position.y, intersection.position.z, 30);
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
