#include <sys/stat.h>
#include <limits>
#include <iostream>
#include <SDL.h>
#include <random>
#include <chrono>
#include <math.h>
#include <vector>
#include <thread>
#include <random>
#include <Windows.h>
#include <filesystem>
#include <cstring>
#include <string>
#include <fstream>

#include "SDLInputManager.h"
#include "Common.hpp"
#include "Object.hpp"
#include "Scene.hpp"
#include "imgui_impl_sdlrenderer.h"
#include "imgui_impl_sdl2.h"
#include "tinyfiledialogs.h"


#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define THREADS 16

float SCREEN_SCALE = .5;
int FOV = 55;
int MAXBOUNCES = 2;
int TARGETFRAMES = 4096;
int ACCUMULATIONFRAMES = 1;
bool SIMPLEDRAW = true;

#undef main

using std::cout; using std::endl;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;


int frameCount = 0;
float progressiveResolutionScaler;
bool setFrame = false;

SDL_Surface* renderSurface;
SDL_Texture* renderTexture;
SDL_Renderer* renderer;
Object* selectedObject = NULL;

float3 SunDirection = float3(1, -1, -1);
Color SkyColor = Color(.2, .35, 1.0f) * 10.0f;
Color HorizonColor = Color(1.0, 0.9f, 0.5f) * 5.0f;
Color GroundColor = Color(.08f, .06f, .03f);
Color SunColor = Color(500, 500, 500);
Color colorBuffer[SCREEN_HEIGHT * SCREEN_WIDTH];
std::vector<Object*> ObjectsToRender;

void SetScreenPixel(const int& x, const int& y, const Color& color) {
	Uint32* targetPixel = (Uint32*)((Uint8*)renderSurface->pixels + (SCREEN_HEIGHT - 1 - y) * renderSurface->pitch + x * renderSurface->format->BytesPerPixel);
	if (!setFrame) {
		float weight = 1.0 / ACCUMULATIONFRAMES;
		colorBuffer[x + y * SCREEN_WIDTH] = colorBuffer[x + y * SCREEN_WIDTH] * (1 - weight) + color * weight;
	}
	else {
		colorBuffer[x + y * SCREEN_WIDTH] = color;
	}
	//Color finalColor = Color(x / (float)SCREEN_WIDTH, y / (float)SCREEN_HEIGHT,0);
	Color finalColor = colorBuffer[x + y * SCREEN_WIDTH];
	finalColor = finalColor / (Color(1, 1, 1) + finalColor);
	*targetPixel = finalColor;
}
Color GetEnvironmentColor(const float3& rayDirection) {
	float upd = float3::Dot(rayDirection, WORLDUP);
	Color Sun = float3::Dot(rayDirection, SunDirection * -1) > 0.99 ? SunColor : Color(0, 0, 0);
	if (upd > 0) {
		Color t = Color::Lerp(HorizonColor, SkyColor, powf(upd, 0.1f));
		t = Color::Lerp(t, SkyColor * 0.1f, upd);
		return t + Sun;
	}
	else {
		upd = fabsf(upd);
		return Color::Lerp(HorizonColor, GroundColor, powf(upd, .05f)) + Sun;
	}
}
float3 GetRandomDirection() {
	float3 sr;
	do {
		sr.x = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr.y = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr.z = ((float)rand() / RAND_MAX - 0.5f) * 2;
	} while (sr.x * sr.y * sr.z > 1);
	return sr.Normalized();
}
float3 GetRandomNormalOrientedHemisphere(const float3& normal) {
	float3 sr = GetRandomDirection();
	if (float3::Dot(sr, normal) < 0) {
		sr = sr * -1;
	}
	return sr;
}
float3 GetRayDirection(const Transform& cameraTransform, const int& pixelX, const int& pixelY) {
	const float clipDistance = .01f;

	float nX = (pixelX / (float)SCREEN_WIDTH) * 2 - 1;
	float nY = (pixelY / (float)SCREEN_HEIGHT) * 2 - 1;
	float aspecRatio = (float)SCREEN_WIDTH / (float)SCREEN_HEIGHT;
	float hFov = FOV * M_PI / 180.0f;
	float3 forwardDir = cameraTransform.forward * clipDistance;
	float rd = (clipDistance * tanf(hFov / 2.0f)) * aspecRatio;
	float ld = (clipDistance * tanf(hFov / 2.0f));
	float3 u = cameraTransform.right * rd * nX;
	float3 v = cameraTransform.up * ld * nY;

	float3 pos = (u + v + forwardDir).Normalized();

	return pos;
}
RayHitObject GetClosestObject(const float3& rayOrigin, const float3& rayDirection) {
	int objectCount = ObjectsToRender.size();
	RayHitObject returnResults;
	float shortestDistance = std::numeric_limits<float>().infinity();
	for (size_t i = 0; i < objectCount; i++)
	{
		Rayhit hitResults = ObjectsToRender[i]->Raytrace(rayOrigin, rayDirection);
		if (hitResults.valid) {

			if (hitResults.distance < shortestDistance) {
				returnResults.objectReference = ObjectsToRender[i];
				shortestDistance = hitResults.distance;
				returnResults.rayHit = hitResults;
			}
		}
	}
	return returnResults;
}
Color RaytraceScene(const float3& rayOrigin, const float3& rayDirection) {
	auto hit = GetClosestObject(rayOrigin, rayDirection);
	if (!hit.rayHit.valid) {
		return GetEnvironmentColor(rayDirection);
	}

	if (SIMPLEDRAW) {
		Color reflectedColor = GetEnvironmentColor(rayDirection.Reflect(hit.rayHit.normal));
		float k = hit.objectReference->material.SpecularAmount;
		float s = hit.objectReference->material.Smoothness;

		float fresnal = 0;
		if (hit.objectReference == selectedObject) {
			fresnal = 1 - float3::Dot(-hit.rayHit.normal, rayDirection);
			fresnal = max(fresnal, 0.0f);
			fresnal = smoothstep(0.0f, 0.5f, fresnal);
		}

		return Color::Lerp(hit.objectReference->material.BaseColor * (1-k) + reflectedColor * k * s + hit.objectReference->material.EmissiveColor, Color(3,3,0), fresnal);
	}

	Color incomingLight = hit.objectReference->material.EmissiveColor;
	Color hitColor = hit.objectReference->material.BaseColor;
	float3 sray = rayDirection;
	bool specularProb = hit.objectReference->material.SpecularAmount >= ((float)rand() / RAND_MAX);
	const float lightEnergyDissipation = 0.8f;
	for (int i = 0; i < MAXBOUNCES; i++)
	{
		if (i != 0) {
			hitColor *= lightEnergyDissipation;
		}
		float3 reflectedRay = sray.Reflect(hit.rayHit.normal);
		//sray = (hit.rayHit.normal + GetRandomDirection()).Normalized();
		sray = GetRandomNormalOrientedHemisphere(hit.rayHit.normal);
		sray = float3::Lerp(sray, reflectedRay, hit.objectReference->material.Smoothness * specularProb);
		sray = sray.Normalized();
		hit = GetClosestObject(hit.rayHit.point + hit.rayHit.normal * .00001f, sray);
		if (!hit.rayHit.valid) {
			incomingLight += GetEnvironmentColor(sray) * hitColor;
			break;
		}
		specularProb = hit.objectReference->material.SpecularAmount >= ((float)rand() / RAND_MAX);
		incomingLight += hit.objectReference->material.EmissiveColor * hitColor;
		hitColor *= Color::Lerp(hit.objectReference->material.BaseColor, hit.objectReference->material.SpecularColor, specularProb);
	}

	//float3 lightVector = float3(0, 1, 0).Normalized();
	//Color lightColor = SunColor;
	//float3 viewDirection = rayDirection * -1;
	//float3 halfVector = (viewDirection + lightVector).Normalized();
	//float baseReflectance = .25;
	//auto sray = GetRandomNormalOrientedHemisphere(hit.normal);
	//auto reflectionRay = rayDirection.Reflect(hit.normal);
	//reflectionRay = float3::Lerp(reflectionRay, sray, roughness * roughness).Normalized();
	//auto reflectionObject = GetClosestObject(hit.point + hit.normal * 0.05f, reflectionRay);
	//if (reflectionObject.valid) {
	//	reflectionColor = Color(0, 0, 0);
	//}
	//else {
	//	reflectionColor = get_environment_color(reflectionRay);
	//}

	//cook torrence specular lighting: fcoock-torrence = DGF / (4 * <V, N> * <L, N>)
	//float cookTorrenceNumerator = PBRModel::Fresnal(baseReflectance, viewDirection, halfVector) * PBRModel::NormaliDistribution(roughness, halfVector, hit.normal) * PBRModel::GeometricAttenuation(roughness, hit.normal, lightVector) * PBRModel::GeometricAttenuation(roughness, hit.normal, viewDirection);
	//float cookTorrenceDenominator = max(float3::Dot(viewDirection, hit.normal, true) * float3::Dot(lightVector, hit.normal, true) * 4.0f, .000001f);
	//float3 cookTorrence = float3(1,1,1) * cookTorrenceNumerator / cookTorrenceDenominator;
	//float ks = flerpf(PBRModel::Fresnal(baseReflectance, viewDirection, halfVector), 1, metalness);
	//float kd = 1 - ks;
	//float3 brdf = float3(diffusedColor.r, diffusedColor.g, diffusedColor.b) * kd + cookTorrence * ks;
	//Color outgoing = emmission + brdf + lightColor * float3::Dot(lightVector, hit.normal, true);

	return incomingLight;
}

#pragma region ThreadUniqueRandFunctions
std::uniform_int_distribution<int>* dice_distribution = new std::uniform_int_distribution<int>(0, INT_MAX);
auto randFunctions = new std::_Binder<std::remove_cv<std::_Unforced>::type, std::uniform_int_distribution<int>&, std::mt19937&>*[THREADS];
bool suspendAllThreads = false;
auto renderWorkers = new std::thread * [THREADS];
auto threadGroupStatus = new bool[THREADS];
#pragma endregion

void renderArea(
	unsigned int index,
	unsigned int minX, unsigned int maxX, 
	unsigned int minY, unsigned int maxY, 
	const Transform* camera) {

	//min inclusive
	//max exclusive
	while (!suspendAllThreads) {
		if (threadGroupStatus[index] == false) {
			int steps = ceil(1 / ((float)SCREEN_SCALE * progressiveResolutionScaler));
			Color color;
			for (size_t i = minX; i < maxX; i += steps)
			{
				for (size_t j = minY; j < maxY; j += steps)
				{
					float3 rayDirection = GetRayDirection(*camera, i, j);
					color = RaytraceScene((*camera).position, rayDirection);

					for (size_t i1 = 0; i1 < steps && i + i1 < maxX; i1++)
					{
						for (size_t j1 = 0; j1 < steps && j + j1 < maxY; j1++)
						{
							SetScreenPixel(i + i1, j + j1, color);
						}
					}
				}
			}
			threadGroupStatus[index] = true;
		}
		else {
			Sleep(1);
		}
	}
};

int main()
{

#pragma region App Init
	srand(0);
	SunDirection = SunDirection.Normalized();
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window* window = SDL_CreateWindow("Raytracer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	renderSurface = SDL_GetWindowSurface(window);

	SDL_FillRect(renderSurface, 0, Color(0, 0, 0, 0));
	progressiveResolutionScaler = 1;

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	io.ConfigDragClickToInputText = true;
	ImGui::StyleColorsDark();
	ImGui_ImplSDL2_InitForSDLRenderer(window, renderer);
	ImGui_ImplSDLRenderer_Init(renderer);
	
	std::string sceneFolder = "./Scenes";
	struct stat sb;
	if (stat(sceneFolder.c_str(), &sb) != 0) {
		CreateDirectoryA(sceneFolder.c_str(), NULL);
	}

#pragma endregion

#pragma region SceneSetup

	Scene scene1("./Scenes/Scene1.json");
	scene1.Load();
	ObjectsToRender = scene1.GetObjects();

	Transform camera;
	camera.position = float3(0, 0, 0);
	camera.RotateAboutAxis(0, WORDLRIGHT);

	/*Object* d;
	int s = 8;
	int zOffset = 10;
	float half = s / 2;
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
		{
			d = new Sphere(0.2, float3((i - half) * 0.6, -1, (j + zOffset) * 0.3));
			d->material = Material(Color(
				(rand() % 256) / (255.0f),
				(rand() % 256) / (255.0f),
				(rand() % 256) / (255.0f)
			), Color(0, 0, 0), 1, 0);
			scene1.AddObject(d);
		}
	}
	d = new Sphere(1, float3(0, 0, 5));
	d->material = Material(Color(1, 1, 1), Color(0, 0, 0), 0.5, 0);
	scene1.AddObject(d);
	d = new Sphere(2, float3(4, 4, 8));
	d->material = Material(Color(1, 1, 1), Color(50, 50, 50), 0.5, 0);
	scene1.AddObject(d);
	d = new Sphere(1000, float3(0, -1001.2f, 5));
	d->material = Material(Color(1, 1, 1), Color(0, 0, 0), 0.5, 0);
	scene1.AddObject(d);
	scene1.Save();*/

#pragma endregion

#pragma region ThreadSetup
	int div = ceil(SCREEN_WIDTH / (THREADS)) + 1;
	for (size_t i = 0; i < THREADS; i++)
	{
		threadGroupStatus[i] = false;
		auto mt1 = std::mt19937();
		mt1.seed(i);
		auto d = std::bind(*dice_distribution, mt1);

		randFunctions[i] = new std::_Binder<std::remove_cv<std::_Unforced>::type, std::uniform_int_distribution<int>&, std::mt19937&>(d);
		int initialX = div * i;
		int nextX = min(initialX + div, SCREEN_WIDTH);
		renderWorkers[i] = new std::thread(renderArea, i, initialX, nextX, 0, SCREEN_HEIGHT, &camera);
	}
#pragma endregion

	long long frametime = 0;
	double totalframetime = 0;
	float fps = 0;
	bool pause = false;
	bool doSetFrame = false;
	float delta = 0;
	int mouseX = 0;
	int mouseY = 0;
	float mouseSpeed = .08;
	float moveSpeed = 1;
	SDL_Event e;
	auto inputManager = SDLInputManager(&e);

	char savePath[512];
	strcpy_s(savePath, scene1.GetFilePath().c_str());
	char loadPath[512];
	strcpy_s(loadPath, scene1.GetFilePath().c_str());

	std::chrono::steady_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	while (true) {
		SDL_PumpEvents();
		inputManager.SetState();
		ImGui_ImplSDL2_ProcessEvent(&e);

		if (e.window.event == SDL_WINDOWEVENT_CLOSE) {
			break;
		}

		//check if all thread groups are done working on the render
		bool threadGate = false;
		for (size_t i = 0; i < THREADS; i++)
		{
			if (threadGroupStatus[i] == false) {
				threadGate = true;
				break;
			}
		}
		if (threadGate) {
			continue;
		}
		// thread safe after this point


		if (inputManager.OnKeyDown(SDL_SCANCODE_P)) {
			pause = !pause;
		}
		inputManager.GetMouseInput(mouseX, mouseY);
		if (inputManager.OnMouse(SDL_BUTTON_RIGHT)) {
			doSetFrame = true;
			camera.RotateAboutAxis(mouseX * mouseSpeed * 0.03, WORLDUP);
			camera.RotateAboutAxis(mouseY * mouseSpeed * 0.03, camera.right);
		}

		ImGui_ImplSDLRenderer_NewFrame();
		ImGui_ImplSDL2_NewFrame();
		ImGui::NewFrame();
		{
			static float f = 0.0f;
			static int counter = 0;
			ImGui::Begin("Inspector");  
			if (ImGui::BeginMenu("Scene File"))
			{
				if (ImGui::MenuItem("Open..", "Ctrl+O")) { 
					//TCHAR wpath[MAX_PATH];
					//char path[MAX_PATH];
					//size_t nNumCharConverted;
					//GetCurrentDirectory(MAX_PATH, wpath);
					//wcstombs_s(&nNumCharConverted, path, MAX_PATH, wpath, MAX_PATH);
					const char* ext = "*.json";
					const char* const* extension = &ext;
					char* filePath = tinyfd_openFileDialog("Open Scene", loadPath, 1, extension, NULL, 0);
					if (filePath != NULL) {
						strcpy_s(loadPath, filePath);
						scene1.Unload();
						scene1 = Scene(std::string(loadPath));
						scene1.Load();
						ObjectsToRender = scene1.GetObjects();
						doSetFrame = true;
						selectedObject = NULL;
					}
				}
				if (ImGui::MenuItem("Save", "Ctrl+S")) { 
					const char* ext = "*.json";
					const char* const* extension = &ext;
					char* filePath = tinyfd_saveFileDialog("Open Scene", savePath, 1, extension, NULL);
					if (filePath != NULL) {
						strcpy_s(savePath, filePath);
						scene1.SaveAs(savePath);
					}
				}
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("Create"))
			{
				if (ImGui::MenuItem("Cube")) {
					selectedObject = new Box(float3(1, 1, 1));
					selectedObject->transform.position = camera.position + camera.forward * 5;
					ObjectsToRender.push_back(selectedObject);
					scene1.AddObject(selectedObject);
				}
				if (ImGui::MenuItem("Sphere")) {
					selectedObject = new Sphere(.5f, camera.position + camera.forward * 5);
					ObjectsToRender.push_back(selectedObject);
					scene1.AddObject(selectedObject);
				}
				ImGui::EndMenu();
			}
			if (selectedObject != NULL) {
				if (ImGui::CollapsingHeader("Object Properties")) {
					doSetFrame = true;
					selectedObject->OnGUI();
				}
			}
			if (ImGui::CollapsingHeader("Settings")) {
				if (ImGui::Button("Switch Render Mode")) {
					SIMPLEDRAW = !SIMPLEDRAW;
					doSetFrame = true;
				}

				ImGui::InputInt("Max Frames", &TARGETFRAMES);
				int beforeValue = FOV;
				ImGui::SliderInt("FOV", &beforeValue, 15, 103);
				if (beforeValue != FOV) {
					FOV = beforeValue;
					doSetFrame = true;
				}
				beforeValue = MAXBOUNCES;
				ImGui::InputInt("Light Bounces", &beforeValue);
				if (beforeValue != MAXBOUNCES) {
					MAXBOUNCES = beforeValue;
					MAXBOUNCES = max(MAXBOUNCES, 0);
					doSetFrame = true;
				}

				ImGui::SliderFloat("Render Scale", &SCREEN_SCALE, 0.25f, 1.0f);
				if (SIMPLEDRAW) {
					SCREEN_SCALE = SDL_clamp(SCREEN_SCALE, 0.25f, 0.5f);
				}
			}

			ImGui::Text("Application average %.3f \nms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
			ImGui::End();
		}

		if (!io.WantCaptureMouse && !io.WantCaptureKeyboard) {
			if (selectedObject != NULL) {
				if (inputManager.OnKeyDown(SDL_SCANCODE_DELETE)) {
					scene1.RemoveObject(selectedObject);
					delete_from_vector(ObjectsToRender, selectedObject);
					delete selectedObject;
					selectedObject = NULL;
					doSetFrame = true;
				}
			}
			float speed = moveSpeed * delta;
			float3 previousPos = camera.position;
			if (inputManager.OnKey(SDL_SCANCODE_LSHIFT)) {
				speed = 2 * delta;
			}
			if (inputManager.OnKey(SDL_SCANCODE_W)) {
				camera.position = camera.position + camera.forward * speed;
			}
			if (inputManager.OnKey(SDL_SCANCODE_D)) {
				camera.position = camera.position + camera.right * speed;
			}
			if (inputManager.OnKey(SDL_SCANCODE_A)) {
				camera.position = camera.position - camera.right * speed;
			}
			if (inputManager.OnKey(SDL_SCANCODE_S)) {
				camera.position = camera.position - camera.forward * speed;
			}
			if (inputManager.OnKey(SDL_SCANCODE_E)) {
				camera.position = camera.position + camera.up * speed;
			}
			if (inputManager.OnKey(SDL_SCANCODE_Q)) {
				camera.position = camera.position - camera.up * speed;
			}
			if (camera.position != previousPos) {
				doSetFrame = true;
			}
			if (inputManager.OnMouseDown(SDL_BUTTON_LEFT)) {
				if (selectedObject != NULL) {
					selectedObject = NULL;
				}
				else {
					int x, y;
					inputManager.GetMouseScreenPosition(x, y);
					y = SCREEN_HEIGHT - y;

					auto hit = GetClosestObject(camera.position, GetRayDirection(camera, x, y));

					if (hit.rayHit.valid) {
						selectedObject = hit.objectReference;
					}
					else {
						selectedObject = NULL;
					}
				}
			}

		}



		renderTexture = SDL_CreateTextureFromSurface(renderer, renderSurface);
		ImGui::Render();
		SDL_RenderSetScale(renderer, io.DisplayFramebufferScale.x, io.DisplayFramebufferScale.y);
		SDL_RenderClear(renderer);
		SDL_RenderCopy(renderer, renderTexture, NULL, NULL);
		ImGui_ImplSDLRenderer_RenderDrawData(ImGui::GetDrawData());
		SDL_RenderPresent(renderer);
		SDL_DestroyTexture(renderTexture);

		auto t2 = std::chrono::high_resolution_clock::now();
		frametime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		delta = frametime / 1000.0f;
		totalframetime += SIMPLEDRAW ? 0 : delta;
		fps = 1000.0f / frametime;
		std::string titleContent = "";
		titleContent += "fps: " + std::to_string(fps) + " | ";
		titleContent += "total time (seconds): " + std::to_string(totalframetime) + " | ";
		titleContent += "ACCUMULATIONFRAMES: " + std::to_string(ACCUMULATIONFRAMES) + " | ";
		SDL_SetWindowTitle(window, titleContent.c_str());
		inputManager.ResetState();

		t1 = std::chrono::high_resolution_clock::now();
		
		if (pause || ACCUMULATIONFRAMES == TARGETFRAMES) {
			continue;
		}

		if (doSetFrame) {
			setFrame = true;
			ACCUMULATIONFRAMES = 1;
			totalframetime = 0;
			progressiveResolutionScaler = 1.0f / 4.0f;
			doSetFrame = false;
		}
		else {
			setFrame = false;
			if (progressiveResolutionScaler != 1) {
				setFrame = true;
			}
			progressiveResolutionScaler = 1;
			ACCUMULATIONFRAMES += SIMPLEDRAW ? 0 : 1;
		}
		//resume thread render
		for (size_t i = 0; i < THREADS; i++)
		{
			threadGroupStatus[i] = false;
		}
	}

	suspendAllThreads = true;
	for (size_t i = 0; i < THREADS; i++)
	{
		renderWorkers[i]->join();
		delete renderWorkers[i];
		delete randFunctions[i];
	}
	delete[] renderWorkers;
	delete[] threadGroupStatus;
	delete[] randFunctions;

	ImGui_ImplSDLRenderer_Shutdown();
	ImGui_ImplSDL2_Shutdown();
	ImGui::DestroyContext();
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}