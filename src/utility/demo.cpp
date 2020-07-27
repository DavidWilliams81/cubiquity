#include "demo.h"

#include "rendering.h"
#include "utility.h"

#include <filesystem>
#include <iostream>
#include <vector>

using namespace Cubiquity;
namespace fs = std::filesystem;

// See https ://stackoverflow.com/a/11142540
void findVolumes(const fs::path& root, std::vector<fs::path>& paths)
{
	std::cout << "\t Searching " << root << std::endl;
	if (!fs::exists(root) || !fs::is_directory(root)) return;

	fs::recursive_directory_iterator it(root);
	fs::recursive_directory_iterator endit;

	while (it != endit)
	{
		if (fs::is_regular_file(*it) && it->path().extension() == ".vol")
		{
			paths.push_back(it->path());
		}
		++it;
	}
}

Demo::Demo(int argc, char** argv, WindowType windowType)
	: Window(windowType)
	, mArgs(argc, argv)
{
}

void Demo::onInitialise()
{
	std::string path;
	if(mArgs.positional().size() >= 1)
	{
		path = mArgs.positional()[0];
	}
	else
	{
		std::cout << "No volume specified - searching filesystem..." << std::endl;
		std::vector<fs::path> paths;
		findVolumes("data", paths);
		findVolumes("../data", paths);

		if (paths.size() > 0)
		{
			std::cout << std::endl << "Found the following volumes:" << std::endl;
			for (uint ct = 0; ct < paths.size(); ct++)
			{
				std::cout << "\t" << ct + 1 << ". " << paths[ct] << std::endl;
			}
			int pathIndex;
			do
			{
				std::cout << std::endl << "Choose a volume (enter a number 1 - " << paths.size() << "): ";
				std::string input;
				std::cin >> input;
				try
				{
					pathIndex = std::stoi(input);
				}
				catch (std::exception& e)
				{
					continue;
				}
			} while (pathIndex < 1 || pathIndex > paths.size());
			path = paths[pathIndex - 1].string();
		}
		else
		{
			std::cout << "No volumes found!" << std::endl;
		}
	}

	if (!path.empty())
	{
		std::cout << "Opening volume \'" << path << "\'... ";
		if (!mVolume.load(path))
		{
			std::cout << " failed to open volume!" << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << " done" << std::endl;

		Box3i bounds = computeBounds(mVolume, [](MaterialId matId) { return matId != 0; });

		// If our volume consists of empty space carved out of a solid volume then the bounds will just be the entire
		// (huge) volume which is not useful for deciding where to put the camera. We hard-code the camera to the
		// origin in this case.
		if (static_cast<uint64_t>(bounds.upper().x()) - static_cast<uint64_t>(bounds.lower().x()) > 1000000000)
		{
			std::cout << "Invalid bounds (probably the quake map?" << std::endl;
			mCamera.position = Vector3f(0.0f, 0.0f, 0.0f);
		}
		else
		{
			Vector3f centre = bounds.centre();
			float halfDiagonal = length((bounds.upper() - bounds.lower())) * 0.5f;

			// Note that the last param uses '-' instead of '+' This is just a hack for our current test data.
			mCamera.position = Vector3f(centre.x() + halfDiagonal, centre.y() + halfDiagonal, centre.z() - halfDiagonal);
		}

		// Hacky values found to work with our current test dta
		mCamera.yaw = Pi + ((3*Pi) / 4.0f);
		mCamera.pitch = -(Pi / 4.7f);
	}
}

void Demo::onUpdate(float deltaTime)
{	
	// Move forward
	if (keyState(SDL_SCANCODE_W) == KeyState::Down)
	{
		mCamera.position += mCamera.forward() * deltaTime * CameraMoveSpeed;
		onCameraModified();
	}
	// Move backward
	if (keyState(SDL_SCANCODE_S) == KeyState::Down)
	{
		mCamera.position -= mCamera.forward() * deltaTime * CameraMoveSpeed;
		onCameraModified();
	}
	// Strafe right
	if (keyState(SDL_SCANCODE_D) == KeyState::Down)
	{
		mCamera.position += mCamera.right() * deltaTime * CameraMoveSpeed;
		onCameraModified();
	}
	// Strafe left
	if (keyState(SDL_SCANCODE_A) == KeyState::Down)
	{
		mCamera.position -= mCamera.right() * deltaTime * CameraMoveSpeed;
		onCameraModified();
	}
}

void Demo::onKeyUp(const SDL_KeyboardEvent& event)
{
	if (event.keysym.sym == SDLK_ESCAPE)
	{
		close();
	}
}

void Demo::onMouseMotion(const SDL_MouseMotionEvent& event)
{
	if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		// Compute new orientation
		mCamera.yaw -= CameraTurnSpeed * event.xrel;
		mCamera.pitch -= CameraTurnSpeed * event.yrel;
		onCameraModified();
	}
}

void Demo::onMouseButtonDown(const SDL_MouseButtonEvent& event)
{
	if (event.button == SDL_BUTTON(SDL_BUTTON_LEFT))
	//if (mouseButtonState(SDL_BUTTON_LEFT) == MouseButtonState::Down)
	{
		//std::cout << "x = " << event.x << ", y = " << event.y << std::endl;

		Ray3d ray = mCamera.rayFromViewportPos(event.x, event.y, width(), height());

		RayVolumeIntersection intersection = ray_parameter(mVolume, ray);
		if (intersection)
		{
			//Vector3d surfaceColour = decodeMaterial(intersection.material);
			mVolume.setVoxel(static_cast<Vector3i>(intersection.position), 0); // FIXME - Should round psition?
			Vector3i centre = static_cast<Vector3i>(intersection.position);
			for (int32 z = centre.z() - 30; z < centre.z() + 30; z++)
			{
				for (int32 y = centre.y() - 30; y < centre.y() + 30; y++)
				{
					for (int32 x = centre.x() - 30; x < centre.x() + 30; x++)
					{
						Vector3d toCentre = intersection.position - Vector3d(x, y, z);
						if (dot(toCentre, toCentre) < 900)
						{
							mVolume.setVoxel(x, y, z, 0);
						}
					}
				}
			}

			onVolumeModified();
		}
	}

	if(event.button == SDL_BUTTON(SDL_BUTTON_RIGHT))
	//if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		SDL_SetRelativeMouseMode(SDL_TRUE);
	}
}

void Demo::onMouseButtonUp(const SDL_MouseButtonEvent& event)
{
	if (event.button == SDL_BUTTON(SDL_BUTTON_RIGHT))
	//if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		SDL_SetRelativeMouseMode(SDL_FALSE);
	}
}
