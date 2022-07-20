#include <limits>
#include <iostream>
#include <SDL.h>
#include <random>
#include <chrono>
#include <math.h>
#include <vector>
#include <thread>

#undef main

using std::cout; using std::endl;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

#define WORDLRIGHT float3(1,0,0)
#define WORLDUP float3(0,1,0)
#define WORLDFORWARD float3(0,0,1)
#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
SDL_Surface* MAINSCREEN;

struct float3 {
	float x, y, z;
	float3(float x, float y, float z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	float Magnitude() const {
		return sqrtf(x * x + y * y + z * z);
	}
	float MagnitudeSqrd() const {
		return (x * x + y * y + z * z);
	}
	static float Dot(const float3& lhs, const float3& rhs) {
		return (lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
	}
	static float3 Cross(const float3& lhs, const float3& rhs) {
		return float3(lhs.y * rhs.z - rhs.y * lhs.z, rhs.x * lhs.z - lhs.x * rhs.z, lhs.x * rhs.y - rhs.x * lhs.y);
	}
	float3() {
		x = 0;
		y = 0;
		z = 0;
	}
	float3 operator-(const float3& rhs) const {
		return float3(x - rhs.x, y - rhs.y, z - rhs.z);
	}
	float3 operator+(const float3& rhs) const {
		return float3(x + rhs.x, y + rhs.y, z + rhs.z);
	}
	float3 operator*(const float& rhs) const {
		return float3(x * rhs, y * rhs, z * rhs);
	}
	float3 operator/(const float& rhs) const {
		return float3(x / rhs, y / rhs, z / rhs);
	}
	float3 Normalized() {
		float length = Magnitude();
		return float3(x / length, y / length, z / length);
	}
	float3 Reflect(const float3& normal) const {
		return *this - normal * (2 * (float3::Dot(*this, normal)));
	}
};
struct Color {
	//1 = 255
	//0 = 0
	float r;
	float g;
	float b;
	float a;

private:
	Uint32 fromRGBA() const {
		auto sr = (int)(r * 255);
		auto sg = (int)(g * 255);
		auto sb = (int)(b * 255);
		auto sa = (int)(a * 255);

		if (sr > 255) sr = 255;
		if (sg > 255) sg = 255;
		if (sb > 255) sb = 255;
		if (sa > 255) sa = 255;

		auto sr8 = (Uint8)sr;
		auto sg8 = (Uint8)sg;
		auto sb8 = (Uint8)sb;
		auto sa8 = (Uint8)sa;

		return (Uint32)(sa8 << 24 | sr8 << 16 | sg8 << 8 | sb8);
	}
public:
	operator Uint32() const { return fromRGBA(); }
	Color operator *(const float& rhs) const {
		return Color(r * rhs, g * rhs, b * rhs, a * rhs);
	}
	Color(float r, float g, float b, float a = 0) {
		if (r < 0) r = 0;
		if (g < 0) g = 0;
		if (b < 0) b = 0;
		if (a < 0) a = 0;
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
	}
	Color(Uint32 val) {
		a = (float)((val >> 24) & 0xff);
		r = (float)((val >> 16) & 0xff);
		g = (float)((val >> 8) & 0xff);
		b = (float)((val & 0xff));
	}
	Color(float3 vector) : Color(vector.x, vector.y, vector.z) {

	}
	Color() {
		r = 1;
		g = 1;
		b = 1;
		a = 1;
	}

	static Color Lerp(const Color& a, const Color& b, float time) {
		time > 1 ? 1 : time;
		time < 0 ? 0 : time;
		return Color(a.r * (1 - time) + b.r * time, a.g * (1 - time) + b.g * time, a.b * (1 - time) + b.b * time, a.a * (1 - time) + b.a * time);
	}
};
struct Transform {
	float3 right = float3(1, 0, 0);
	float3 up = float3(0, 1, 0);
	float3 forward = float3(0, 0, 1);
	float3 position = float3(0, 0, 0);
	float3 scale = float3(1, 1, 1);
	void RotateAboutAxis(float angle, float3 axis) {
		forward = forward * cosf(angle) + float3::Cross(axis, forward) * sinf(angle) + axis * float3::Dot(axis, forward) * (1 - cosf(angle));
		up = up * cosf(angle) + float3::Cross(axis, up) * sinf(angle) + axis * float3::Dot(axis, up) * (1 - cosf(angle));
		right = right * cosf(angle) + float3::Cross(axis, right) * sinf(angle) + axis * float3::Dot(axis, right) * (1 - cosf(angle));
	}
};
struct Rayhit {
	bool valid = false;
	Color color;
	float3 normal;
	float3 point;
	float distance;
};

float3 SunDirection = float3(0, 1, .5);
Color SkyColor = Color(.5, 1, 2);
Color GroundColor = Color(.05f, .06f, .1f);
Color SunColor = Color(1.0f, 1.0f, 1.0f);


void set_pixel(const int& x, const int& y, const Color& color, const SDL_Surface* surface) {

	Uint32* targetPixel = (Uint32*)((Uint8*)surface->pixels + (SCREEN_HEIGHT - 1 - y) * surface->pitch + x * surface->format->BytesPerPixel);

	*targetPixel = color;
}

Color get_environment_color(const float3& rayDirection) {
	auto halfColor = Color::Lerp(GroundColor, SkyColor, 0.5f);
	float upd = float3::Dot(rayDirection, WORLDUP);
	if (upd > 0) {
		return Color::Lerp(halfColor, SkyColor, powf(upd, 0.25f));
	}
	else {
		upd = fabsf(upd);
		return Color::Lerp(halfColor, GroundColor, powf(upd, 0.25f));
	}
}

class Object {
private:
	int myIndex;
	static Rayhit GetClosestObject(const float3& rayOrigin, const float3& rayDirection) {
		int objectCount = allObjects.size();
		Rayhit returnResults;
		float shortestDistance = std::numeric_limits<float>().infinity();
		for (size_t i = 0; i < objectCount; i++)
		{
			Rayhit hitResults = allObjects[i]->Raytrace(rayOrigin, rayDirection);
			if (hitResults.valid) {

				if (hitResults.distance < shortestDistance) {
					shortestDistance = hitResults.distance;
					returnResults = hitResults;
				}
			}
		}
		return returnResults;
	}
	static float3 UP;
	static float3 RandomNormalOrientedHemisphere(const float3& normal, const float3& rayDir) {

		//return rayDir - normal * (2 * (float3::Dot(rayDir, normal)));

		float3 sr;
		sr.x = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr.y = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr.y = sr.y < 0 ? 0 : sr.y;
		sr.z = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr = sr.Normalized();

		//sr creates a random unit vector in a hemisphere oriented up

		//orient sr towards normal




		//auto rotAxis = float3::Cross(normal, UP);
		//float dd = (1 - ((float3::Dot(UP, normal) + 1) / 2)) * M_PI;
		//return sr * cosf(dd) + ( rotAxis * (1 - cosf(dd)) * float3::Dot(sr, rotAxis) + (float3::Cross(sr, rotAxis) * sinf(dd)));

		return sr;
	}
protected:
	virtual Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const{
		return Rayhit();
	}

public:

	/// <summary>
	/// Maximum specular reflection bounces.
	/// </summary>
	static int MaxRaytraceBounces;
	/// <summary>
	/// Number of rays to reflect after ray hit
	/// </summary>
	static int SamplesPerRay;

	static std::vector<Object*> allObjects;

	/// <summary>
	/// 0 = Color,
	/// 1 = Normals,
	/// 2 = Depth
	/// </summary>
	static char RaytraceRenderMode;

	/// <summary>
	/// Ray traces from a world space position and ray direction. Calculates the resulting color.
	/// </summary>
	/// <param name="rayOrigin"></param>
	/// <param name="rayDirection"></param>
	/// <returns></returns>
	static Color RaytraceScene(const float3& rayOrigin, const float3& rayDirection) {
		//implement phong shading model
		//get closest object from the origin in the directoin of the ray
		//get hit objct's color at hit point
		//scatter ray around at hit point based on hit objects normal <- diffused scattered light ray
		//if the scattered ray hits an object, it is occluded, if it does not hit, it takes sky light
		//multiply with objects base color
		auto hit = GetClosestObject(rayOrigin, rayDirection);
		if (hit.valid) {

			Color diffusedColor;
			Color reflectionColor;

			for (size_t i = 0; i < SamplesPerRay; i++)
			{
				auto sray = RandomNormalOrientedHemisphere(hit.normal, rayDirection);
				auto shit = GetClosestObject(hit.point + hit.normal * 0.05f, sray);
				if (shit.valid) {
					diffusedColor = Color(0, 0, 0);
				}
				else {
					diffusedColor = get_environment_color(sray);
				}
			}

			diffusedColor = diffusedColor * (1.0f / SamplesPerRay);
			diffusedColor = Color(hit.color.r * diffusedColor.r, hit.color.g * diffusedColor.g, hit.color.b * diffusedColor.b);


			auto reflectionRay = rayDirection.Reflect(hit.normal);
			auto reflectionObject = GetClosestObject(hit.point + hit.normal * 0.05f, reflectionRay);
			if (reflectionObject.valid) {
				reflectionColor = Color(0,0,0);
			}
			else {
				reflectionColor = get_environment_color(reflectionRay);
			}


			return reflectionColor;
		}

		return get_environment_color(rayDirection);
	}


	/// <summary>
	/// 1 = glossy, 0 = rough
	/// </summary>
	float roughness = 0;
	Transform transform;
	/// <summary>
	/// Base color
	/// </summary>
	Color color;


	Object() {
		allObjects.push_back(this);
		myIndex = allObjects.size() - 1;
	}
	virtual ~Object() {
		allObjects.erase(allObjects.begin() + myIndex);
	}
};

float3 Object::UP = float3(0, 1, 0);
int Object::MaxRaytraceBounces = 1;
int Object::SamplesPerRay = 1;
char Object::RaytraceRenderMode = 0;
std::vector<Object*> Object::allObjects = std::vector<Object*>();

class Sphere : public Object {
	float radius;

public:
	const float& GetRadius() {
		return radius;
	}
	void SetRadius(const float& rad) {
		this->radius = rad;
		float t = (rad * rad) / 3.0f;
		t = sqrtf(t);
		transform.scale = float3(t, t, t);
	}
	Sphere(float radius, float3 position) {
		transform.position = position;
		SetRadius(radius);
	}
private:
	bool line_sphere_intersection(
		const float3& rayPosition,
		const float3& rayDirection,
		const float3& spherePosition,
		const float& sphereRadius,
		Color& outColor,
		float3& outNormal,
		float3& outPosition,
		float& outDistance) const 
	{
		outDistance = (std::numeric_limits<float>().infinity());
		float3 directionToSphere = spherePosition - rayPosition;
		float distanceToSphere = directionToSphere.Magnitude();

		float tc = (float3::Dot(directionToSphere, rayDirection));
		tc = abs(tc);

		float3 intersectCheck = rayDirection * tc + rayPosition;
		float squaredRadius = sphereRadius * sphereRadius;

		float3 dir2sphereIntersectDir = intersectCheck - spherePosition;
		float d2 = dir2sphereIntersectDir.MagnitudeSqrd();

		if (d2 > squaredRadius) {
			return false;
		}

		float t1c = sqrtf(squaredRadius - d2);

		float t1 = tc - t1c;
		outDistance = t1;

		outPosition = rayPosition + (rayDirection * t1);
		outNormal = (outPosition - spherePosition).Normalized();

		outColor = color;
		return true;
	}
public:
	Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const override {

		Color color;
		float3 normal, hitPoint;
		float distance;
		Rayhit hitResults;

		if (line_sphere_intersection(rayOrigin, rayDir, transform.position, radius, color, normal, hitPoint, distance)) {
			hitResults.color = color;
			hitResults.normal = normal;
			hitResults.point = hitPoint;
			hitResults.distance = distance;
			hitResults.valid = true;
		}
		return hitResults;
	}
};

float3 GetRayDirection(const Transform& cameraTransform, const int& pixelX, const int& pixelY) {
	const float clipDistance = .01f;

	float nX = (pixelX / (float)SCREEN_WIDTH) * 2 - 1;
	float nY = (pixelY / (float)SCREEN_HEIGHT) * 2 - 1;

	//float hFov = (float)horizontalFovDegrees * M_PI / 180.0f;
	float aspecRatio = (float)SCREEN_WIDTH / (float)SCREEN_HEIGHT;

	float hFov = 70 * M_PI / 180.0f;
	float vFov = 43 * M_PI / 180.0f;


	float3 forwardDir = cameraTransform.forward * clipDistance;
	float rd = (clipDistance * tanf(hFov / 2.0f));
	float ld = (clipDistance * tanf(vFov / 2.0f));
	float3 u = cameraTransform.right * rd;
	float3 v = cameraTransform.up * ld;

	float3 pos = ((u * nX + v * nY) + forwardDir).Normalized();

	return pos;
}

void renderQuadrant (int quadrant, int raytraceSamples, const Transform *camera){

	//min inclusive
	//max exclusive
	int xRangeMin = 0, xRangeMax = 0;
	int yRangeMin = 0, yRangeMax = 0;
	switch (quadrant)
	{
	case(1):
		xRangeMin = (SCREEN_WIDTH / 2);
		xRangeMax = SCREEN_WIDTH;

		yRangeMin = (SCREEN_HEIGHT / 2);
		yRangeMax = SCREEN_HEIGHT;
		break;
	case(2):
		xRangeMin = 0;
		xRangeMax = (SCREEN_WIDTH / 2);

		yRangeMin = (SCREEN_HEIGHT / 2);
		yRangeMax = SCREEN_HEIGHT;
		break;
	case(3):
		xRangeMin = 0;
		xRangeMax = (SCREEN_WIDTH / 2);

		yRangeMin = 0;
		yRangeMax = (SCREEN_HEIGHT / 2);
		break;
	case(4):
		xRangeMin = (SCREEN_WIDTH / 2);
		xRangeMax = SCREEN_WIDTH;

		yRangeMin = 0;
		yRangeMax = (SCREEN_HEIGHT / 2);
		break;
	default:
		return;
	}
	
	
	SDL_Surface* renderTarget = MAINSCREEN;
	void* pixels = renderTarget->pixels;
	int pitch = renderTarget->pitch;
	Uint8 bytesPerPixel = renderTarget->format->BytesPerPixel;

	for (size_t i = 0; i < raytraceSamples; i++)
	{
		int randX = rand() % (xRangeMax - 1) + xRangeMin;
		int randY = rand() % (yRangeMax - 1) + yRangeMin;
		float3 rayDirection = GetRayDirection(*camera, randX, randY);
		auto color = Object::RaytraceScene((*camera).position, rayDirection);
		Uint32* targetPixel = (Uint32*)((Uint8*)pixels + (SCREEN_HEIGHT - 1 - randY) * pitch + randX * bytesPerPixel);
		*targetPixel = color;
	}
};


int main()
{
	srand(0);
	SunDirection = SunDirection.Normalized();
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window* window = SDL_CreateWindow("Raytracer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_VULKAN);
	SDL_Surface* screen = SDL_GetWindowSurface(window);
	MAINSCREEN = screen;
	SDL_FillRect(screen, 0, Color(0, 0, 0, 0));
	auto cursor = SDL_GetDefaultCursor();
	SDL_SetCursor(cursor);
	float frametime = 0;
	float fps = 0;
	SDL_SetRelativeMouseMode(SDL_TRUE);

	Transform camera;
	camera.position = float3(0, 0, 0);
	camera.RotateAboutAxis(0, WORDLRIGHT);

	Color closestColor;
	Color color;
	float3 normal;
	float3 position;
	float renderDistance;

	float time = 0;
	float delta = 0;

	const int raytraceSamples = 4096 * 3;
	int mouseX = 0;
	int mouseY = 0;
	float mouseSpeed = 0.0005;

	Object* d;
	int s = 8;
	int zOffset = 10;
	float half = s / 2;
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
		{
			d = new Sphere(0.2, float3((i - half) * 0.6, -1, (j + zOffset) * 0.3));
			d->color = Color(0.5f, 0.5f, 0.5f);
		}
	}
	d = new Sphere(1, float3(0, 0, 5));
	d->color = Color(.5f, .5f, .5f);
	printf("Finished creating objects");


	bool pause = false;
	while (true) {
		auto t1 = std::chrono::high_resolution_clock::now();

		SDL_PumpEvents();
		const Uint8* state = SDL_GetKeyboardState(NULL);
		if (state[SDL_SCANCODE_ESCAPE]) {
			break;
		}

		//game loop
		if (!pause) {

			SDL_GetRelativeMouseState(&mouseX, &mouseY);
			camera.RotateAboutAxis(mouseX * mouseSpeed, WORLDUP);
			camera.RotateAboutAxis(mouseY * mouseSpeed, camera.right);

			float speed = .5;
			if (state[SDL_SCANCODE_LSHIFT]) {
				speed = 1;
			}
			if (state[SDL_SCANCODE_W]) {
				camera.position = camera.position + camera.forward * delta * speed;
			}
			if (state[SDL_SCANCODE_D]) {
				camera.position = camera.position + camera.right * delta * speed;
			}
			if (state[SDL_SCANCODE_A]) {
				camera.position = camera.position - camera.right * delta * speed;
			}
			if (state[SDL_SCANCODE_S]) {
				camera.position = camera.position - camera.forward * delta * speed;
			}
			if (state[SDL_SCANCODE_E]) {
				camera.position = camera.position + camera.up * delta * speed;
			}
			if (state[SDL_SCANCODE_Q]) {
				camera.position = camera.position - camera.up * delta * speed;
			}


			//for (size_t i = 0; i < SCREEN_WIDTH; i++)
			//{
			//	for (size_t j = 0; j < SCREEN_HEIGHT; j++)
			//	{
			//		int randX = i;
			//		int randY = j;
			//		float3 rayDirection = GetRayDirection(camera, randX, randY);
			//		auto color = Object::RaytraceScene(camera.position, rayDirection);
			//		set_pixel(randX, randY, color);
			//	}
			//}
			//pause = true;

			//create 4 threads
			//each thread will work on different quadrants of the screen
			std::thread q1(renderQuadrant, 1, raytraceSamples, &camera);
			//std::thread q2(renderQuadrant, 2, raytraceSamples, &camera);
			//std::thread q3(renderQuadrant, 3, raytraceSamples, &camera);
			//std::thread q4(renderQuadrant, 4, raytraceSamples, &camera);

			q1.join();
			//q2.join();
			//q3.join();
			//q4.join();

			//for (size_t i = 0; i < raytraceSamples; i++)
			//{
			//	int randX = rand() % SCREEN_WIDTH;
			//	int randY = rand() % SCREEN_HEIGHT;
			//	float3 rayDirection = GetRayDirection(camera, randX, randY);
			//	auto color = Object::RaytraceScene(camera.position, rayDirection);			
			//	//set_pixel(randX, randY, color, MAINSCREEN);
			//
			//	Uint32* targetPixel = (Uint32*)((Uint8*)MAINSCREEN->pixels + (SCREEN_HEIGHT - 1 - randY) * MAINSCREEN->pitch + randX * MAINSCREEN->format->BytesPerPixel);
			//	*targetPixel = color;
			//}

			SDL_UpdateWindowSurface(window);
		}

		auto t2 = std::chrono::high_resolution_clock::now();
		frametime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		delta = frametime / 1000.0f;
		time += delta;
		printf("Time(s): %.2f\n", time);

		fps = (frametime > 0) ? 1000.0f / frametime : 10000;
		printf("frametime(ms): %.0f\n", frametime);
		printf("fps: %.3f\n", fps);

	}



	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}