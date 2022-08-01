#include <limits>
#include <iostream>
#include <SDL.h>
#include <random>
#include <chrono>
#include <math.h>
#include <vector>
#include <thread>
#include <random>
#include <functional>
#include <Windows.h>

#undef main

//#define NOISYRENDER
#define MULTITHREAD

using std::cout; using std::endl;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

const int THREADS = 64;
#define WORDLRIGHT float3(1,0,0)
#define WORLDUP float3(0,1,0)
#define WORLDFORWARD float3(0,0,1)
#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define RENDER_SCALE 1

float flerpf(float a, float b, float t) {
	return a * (1 - t) + b * t;
}
struct float3 {
	float x, y, z;

	inline float operator[](const int& idx) const {
		switch (idx) {
		case(0):
			return x;
		case(1):
			return y;
		case(2):
			return z;
		}
	}
	inline float operator+(const int& idx) const {
		switch (idx) {
		case(0):
			return x;
		case(1):
			return y;
		case(2):
			return z;
		}
	}
	inline float operator()(const int& idx) const {
		switch (idx) {
		case(0):
			return x;
		case(1):
			return y;
		case(2):
			return z;
		}
	}

	inline float3 xz() {
		return float3(x, 0, z);
	}
	inline float3 xy() {
		return float3(x, y, 0);
	}
	inline float3 yz() {
		return float3(0, y, z);
	}


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
	static float Dot(const float3& lhs, const float3& rhs, bool noNeg = false) {

		auto d = (lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
		if (noNeg) {
			return d < 0 ? 0 : d;
		}
		else {
			return d;
		}

	}
	const float3& abs() const {
		return float3(fabsf(x), fabsf(y), fabsf(z));
	}
	static float3 Cross(const float3& lhs, const float3& rhs) {
		return float3(lhs.y * rhs.z - rhs.y * lhs.z, rhs.x * lhs.z - lhs.x * rhs.z, lhs.x * rhs.y - rhs.x * lhs.y);
	}
	static float3 Lerp(const float3& a, const float3& b, float t) {
		return float3(flerpf(a(0), b(0), t), flerpf(a(1), b(1), t), flerpf(a(2), b(2), t));
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


	static float3 ComponentMultiply(float3 lhs, float3 rhs) {
		return float3(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
	}
	static float3 ComponentDivide(float3 lhs, float3 rhs) {
		return float3(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z);
	}

	float3& operator-=(const float3& rhs) {
		*this = (*this) - rhs;
		return *this;
	}
	float3& operator+=(const float3& rhs) {
		*this = (*this) + rhs;
		return *this;
	}
	float3& operator*=(const float& rhs) {
		*this = (*this) * rhs;
		return *this;
	}
	float3& operator/=(const float& rhs) {
		*this = (*this) / rhs;
		return *this;
	}

	float3 operator*(const float3& o) {
		return float3(x * o(0), y * o(1), z * o(2));
	}


	const float3& Normalized() const {
		float length = Magnitude();
		return float3(x / length, y / length, z / length);
	}
	const float3& Reflect(const float3& normal) const {
		return *this - normal * (2 * (float3::Dot(*this, normal)));
	}
};
struct Color {
	//1 = 255
	//0 = 0
	float r = 0;
	float g = 0;
	float b = 0;
	float a = 0;

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
	Color operator *(const Color& rhs) const {
		return Color(r * rhs.r, g * rhs.g, b * rhs.b, a * rhs.a);
	}
	Color operator + (const Color& o) const {
		return Color(r + o.r, g + o.g, b + o.b, a + o.a);
	}
	Color operator / (const float& o) const {
		return Color(r / o, g / o, b / o, a / o);
	}
	Color operator - (const Color& o) const {
		return Color(r - o.r, g - o.g, b - o.b, a - o.a);
	}

	Color& operator += (const Color& o) {

		*this = *this + o;

		return *this;
	}
	Color& operator *=(const float& o) {
		*this = *this * o;

		return *this;
	}
	Color& operator -=(const Color& o) {
		*this = *this - o;

		return *this;
	}
	Color& operator /=(const float& o) {
		*this = *this / o;

		return *this;
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
	void* object;
};

float3 SunDirection = float3(.5, 1, .5);
Color SkyColor = Color(1, 1.2, 1.5);
Color GroundColor = Color(0, 0, 0);
Color SunColor = Color(6, 6, 6);
Color** RenderBuffer;


void set_pixel(const int& x, const int& y, const Color& color, const SDL_Surface* surface) {

	if (x >= SCREEN_WIDTH || y >= SCREEN_HEIGHT || x < 0 || y < 0) {
		return;
	}

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
					returnResults.object = allObjects[i];
				}
			}
		}
		return returnResults;
	}
	static float3 RandomNormalOrientedHemisphere(const float3& normal) {

		//return rayDir - normal * (2 * (float3::Dot(rayDir, normal)));
		float3 sr;
		do {
			sr.x = ((float)rand() / RAND_MAX - 0.5f) * 2;
			sr.y = ((float)rand() / RAND_MAX - 0.5f) * 2;
			sr.z = ((float)rand() / RAND_MAX - 0.5f) * 2;
			sr = sr.Normalized();
		} while (float3::Dot(normal, sr) < 0);


		//sr creates a random unit vector in a hemisphere oriented up

		//orient sr towards normal

		return sr;
	}
protected:
	virtual Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const{
		return Rayhit();
	}

public:

	float smoothness = .5f;
	float metalness = 0;

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

			float smoothness = ((Object*)hit.object)->smoothness;
			float metalness = ((Object*)hit.object)->metalness;

			Color diffusedColor = Color(0,0,0);
			Color reflectionColor;
			Color specularColor = Color(0,0,0);
			Color lightColor = SunColor;
			float3 lighDirection = SunDirection;

			for (size_t i = 0; i < SamplesPerRay; i++)
			{
				//auto sray = RandomNormalOrientedHemisphere(hit.normal);
				auto sray = RandomNormalOrientedHemisphere(hit.normal);
				auto shit = GetClosestObject(hit.point + hit.normal * 0.05f, sray);
				if (shit.valid) {
				
				}
				else {
					auto skyCol = get_environment_color(sray);
					diffusedColor.r += skyCol.r;
					diffusedColor.g += skyCol.g;
					diffusedColor.b += skyCol.b;
					diffusedColor.a += skyCol.a;
				}
			}
			diffusedColor = diffusedColor * (1.0f / SamplesPerRay);
			diffusedColor = Color(hit.color.r * diffusedColor.r, hit.color.g * diffusedColor.g, hit.color.b * diffusedColor.b);
			
			
			//diffusedColor = hit.color;
			float lightDot = float3::Dot(lighDirection, hit.normal, true);

			//auto sray = RandomNormalOrientedHemisphere(hit.normal);
			//auto reflectionRay = rayDirection.Reflect(hit.normal);
			//
			//reflectionRay = reflectionRay + (sray * powf((1 - smoothness), 1));
			//reflectionRay = reflectionRay.Normalized();

			////not pbr, but close enough

			//auto reflectionObject = GetClosestObject(hit.point + hit.normal * 0.05f, reflectionRay);
			//if (reflectionObject.valid) {
			//	reflectionColor = hit.color * (1 - ((Object*)hit.object)->metalness);
			//}
			//else {
			//	reflectionColor = get_environment_color(reflectionRay);
			//	float specDot = float3::Dot(reflectionRay, lighDirection);
			//	specDot = specDot < 0 ? 0 : specDot;
			//	float temp = (specDot - smoothness) * (1 / (1 - smoothness));
			//	temp = temp < 0 ? 0 : temp;

			//	//specular size approximation, custom model
			//	specularColor = lightColor * powf(temp, 2) * (smoothness * smoothness * smoothness);
			//}

			////80% of diffused color is present and 20% of reflections are present
			////as metalness goes up, diffuse approaches 0 and reflections dominate.
			//Color diffuseReflection = (diffusedColor * (0.8f * (1 - metalness)) + (reflectionColor * (0.2f * (1 - metalness) + 1.0f * metalness)));

			//add specular reflection
			return diffusedColor;// +specularColor + (diffusedColor * lightColor) * float3::Dot(lighDirection, hit.normal) * (1 - metalness);
		}

		return get_environment_color(rayDirection);
	}


	/// <summary>
	/// 1 = glossy, 0 = rough
	/// </summary>


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

int Object::MaxRaytraceBounces = 1;
int Object::SamplesPerRay = 128;
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
			hitResults.object = (void*)this;
		}
		return hitResults;
	}
};
class Plane : public Object {
public:
	Plane(float3 size, float3 normal, float3 position) : Object(){
		transform.scale = size;
		transform.forward = float3(normal.x, normal.y, normal.z);
		transform.position = float3(position.x, position.y, position.z);
	}
public:
	Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const override {

		Color color;
		float3 normal, hitPoint;
		float distance;
		Rayhit hitResults;
		
		//float denom = float3::Dot(rayDir, transform.forward);

		float num = float3::Dot(transform.forward, rayOrigin - transform.position);
		float denom = float3::Dot(transform.forward, rayDir);
		float t = (num / denom);

		if (t < 0 && denom < 0) {
			
			float3 tt = (rayDir * t);
			float3 intersection = tt + rayOrigin;

			hitResults.distance = tt.Magnitude();

			hitResults.valid = true;
			hitResults.color = this->color;//// *(hitResults.distance / 100.0f);// / 1.0f);
			hitResults.point = intersection;
			hitResults.object = (void*)this;
			hitResults.normal = transform.forward.Normalized();
			return hitResults;
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


std::uniform_int_distribution<int>* dice_distribution = new std::uniform_int_distribution<int>(0, INT_MAX);
//std::mt19937* rnd1 = new std::mt19937();
//std::mt19937* rnd2 = new std::mt19937();
//std::mt19937* rnd3 = new std::mt19937();
//std::mt19937* rnd4 = new std::mt19937();
//auto rand1 = std::bind(*dice_distribution, *rnd1);
//auto rand2 = std::bind(*dice_distribution, *rnd2);
//auto rand3 = std::bind(*dice_distribution, *rnd3);
//auto rand4 = std::bind(*dice_distribution, *rnd4);
//bool finished1 = false;
//bool finished2 = false;
//bool finished3 = false;
//bool finished4 = false;

bool suspendAllThreads = false;

auto renderWorkers = new std::thread * [THREADS];
auto waitWorkerFlags = new bool[THREADS];
auto randFunctions = new std::_Binder<std::remove_cv<std::_Unforced>::type, std::uniform_int_distribution<int>&, std::mt19937&>*[THREADS];


void renderArea(
	unsigned int index,
	unsigned int minX, unsigned int maxX, 
	unsigned int minY, unsigned int maxY, 
	unsigned int raytraceSamples, const Transform* camera, const SDL_Surface* renderTexture) {

	//min inclusive
	//max exclusive
	//auto rand = randn;
	while (!suspendAllThreads) {
		if (waitWorkerFlags[index] == false) {

#ifdef NOISYRENDER

			for (size_t i = 0; i < raytraceSamples; i++)
			{
				int randX = floor(minX + ((double)(*(randFunctions[index]))() / (double)INT_MAX) * (maxX - minX));
				int randY = floor(minY + ((double)(*(randFunctions[index]))() / (double)INT_MAX) * (maxY - minY));
				float3 rayDirection = GetRayDirection(*camera, randX, randY);
				auto color = Object::RaytraceScene((*camera).position, rayDirection);
				set_pixel(randX, randY, color, renderTexture);
			}

#else

			for (size_t i = minX; i < maxX; i++)
			{
				for (size_t j = minY; j < maxY; j++)
				{
					if (i % (RENDER_SCALE) == 0 && j % (RENDER_SCALE) == 0) {
						float3 rayDirection = GetRayDirection(*camera, i, j);
						auto color = Object::RaytraceScene((*camera).position, rayDirection);
						RenderBuffer[i / RENDER_SCALE][j / RENDER_SCALE] = color;
					}
					set_pixel(i, j, RenderBuffer[i / RENDER_SCALE][j / RENDER_SCALE], renderTexture);

				}
			}

#endif // NOISYRENDER

			waitWorkerFlags[index] = true;
		}
		else {
			Sleep(1);
		}
	}
};

void InitBuffer() {

	if (RenderBuffer != NULL) {
		for (size_t i = 0; i < SCREEN_WIDTH / RENDER_SCALE; i++)
		{
			for (size_t j = 0; j < SCREEN_HEIGHT / RENDER_SCALE; j++)
			{
				delete[] RenderBuffer[i];
			}
		}
		delete[] RenderBuffer;
	}

	RenderBuffer = new Color * [SCREEN_WIDTH / RENDER_SCALE];

	for (size_t i = 0; i < SCREEN_WIDTH / RENDER_SCALE; i++)
	{
		RenderBuffer[i] = new Color[SCREEN_HEIGHT / RENDER_SCALE];
		for (size_t j = 0; j < SCREEN_HEIGHT / RENDER_SCALE; j++)
		{
			RenderBuffer[i][j] = Color(0, 0, 0);
		}
	}
}

int main()
{
	srand(0);
	SunDirection = SunDirection.Normalized();
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window* window = SDL_CreateWindow("Raytracer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_VULKAN);
	SDL_Surface* screen = SDL_GetWindowSurface(window);

	InitBuffer();

	SDL_FillRect(screen, 0, Color(0, 0, 0, 0));
	auto cursor = SDL_GetDefaultCursor();
	SDL_SetCursor(cursor);
	long long frametime = 0;
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


	float delta = 0;

	const unsigned int raytraceSamples = 4096 * 12;
	int mouseX = 0;
	int mouseY = 0;
	float mouseSpeed = .1;

	Object* d;
	int s = 8;
	int zOffset = 10;
	float half = s / 2;
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < s; j++)
		{
			d = new Sphere(0.25, float3((i - half) * 0.5 + .25f, (j - half) * 0.5 + .25, 5));
			d->smoothness = (float)i / (float)s;
			d->metalness = (float)j / (float)s;
			d->color = Color(1, 0, 0);
		}
	}
	
	
	int boxSize = 5;
	
	//floor
	d = new Plane(float3(1,1,1), float3(0, boxSize,0), float3(0, -boxSize, 5));
	d->color = Color(0, 0, 1);

	////ceiling
	//d = new Plane(float3(1, 1, 1), float3(0, -boxSize, 0), float3(0, boxSize, 5));
	//d->color = Color(1, 1, 1);
	//
	////left wall
	//d = new Plane(float3(1, 1, 1), float3(boxSize, 0, 0), float3(-boxSize, 0, 5));
	//d->color = Color(1, 0, 0);
	//
	////right wall
	//d = new Plane(float3(1, 1, 1), float3(-boxSize, 0, 0), float3(boxSize, 0, 5));
	//d->color = Color(0, 1, 0);
	//
	////front wall
	//d = new Plane(float3(1, 1, 1), float3(0, 0, -boxSize), float3(0, 0, boxSize + 5));
	//d->color = Color(1, 1, 1);
	//
	////backwall
	//d = new Plane(float3(1, 1, 1), float3(0, 0, boxSize), float3(0, 0, -boxSize ));
	//d->color = Color(1, 1, 1);


	printf("Finished creating objects");

	int ss = raytraceSamples * 0.25f;
	float time = 0;


	//initialize threads
#ifdef MULTITHREAD
	int div = round(SCREEN_WIDTH / (THREADS));
	for (size_t i = 0; i < THREADS; i++)
	{
		waitWorkerFlags[i] = false;
		auto mt1 = std::mt19937();
		mt1.seed(i);
		auto d = std::bind(*dice_distribution, mt1);

		randFunctions[i] = new std::_Binder<std::remove_cv<std::_Unforced>::type, std::uniform_int_distribution<int>&, std::mt19937&>(d);
		int initialX = div * i;
		renderWorkers[i] = new std::thread(renderArea, i, initialX, initialX + div, 0, SCREEN_HEIGHT, raytraceSamples / THREADS, &camera, screen);
	}
#endif // MULTITHREAD
	
	
	std::chrono::steady_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	bool pause = false;
	while (true) {

#ifndef MULTITHREAD
		t1 = std::chrono::high_resolution_clock::now();
#endif

		SDL_PumpEvents();
		const Uint8* state = SDL_GetKeyboardState(NULL);
		if (state[SDL_SCANCODE_ESCAPE]) {
			break;
		}

		//game loop
		if (!pause) {


#ifndef  MULTITHREAD
			SDL_GetRelativeMouseState(&mouseX, &mouseY);
			camera.RotateAboutAxis(mouseX * mouseSpeed * delta, WORLDUP);
			camera.RotateAboutAxis(mouseY * mouseSpeed * delta, camera.right);

			float speed = delta;
			if (state[SDL_SCANCODE_LSHIFT]) {
				speed = 2 * delta;
			}
			if (state[SDL_SCANCODE_W]) {
				camera.position = camera.position + camera.forward * speed;
			}
			if (state[SDL_SCANCODE_D]) {
				camera.position = camera.position + camera.right * speed;
			}
			if (state[SDL_SCANCODE_A]) {
				camera.position = camera.position - camera.right * speed;
			}
			if (state[SDL_SCANCODE_S]) {
				camera.position = camera.position - camera.forward * speed;
			}
			if (state[SDL_SCANCODE_E]) {
				camera.position = camera.position + camera.up * speed;
			}
			if (state[SDL_SCANCODE_Q]) {
				camera.position = camera.position - camera.up * speed;
			}
#endif // ! MULTITHREAD


#ifdef  MULTITHREAD

			bool threadFlag = false;
			for (size_t i = 0; i < THREADS; i++)
			{
				if (waitWorkerFlags[i] == false) {
					threadFlag = true;
					break;
				}
			}

			if (!threadFlag) {
				SDL_GetRelativeMouseState(&mouseX, &mouseY);
				camera.RotateAboutAxis(mouseX* mouseSpeed* delta, WORLDUP);
				camera.RotateAboutAxis(mouseY* mouseSpeed* delta, camera.right);

				float speed = delta * 5;
				if (state[SDL_SCANCODE_LSHIFT]) {
					speed *= 2;
				}
				float3 moveDir = float3(0,0,0);
				if (state[SDL_SCANCODE_W]) {
					moveDir = camera.forward;
				}
				if (state[SDL_SCANCODE_D]) {
					moveDir += camera.right;
				}
				if (state[SDL_SCANCODE_A]) {
					moveDir -= camera.right;
				}
				if (state[SDL_SCANCODE_S]) {
					moveDir -= camera.forward;
				}
				if (state[SDL_SCANCODE_E]) {
					moveDir += camera.up;
				}
				if (state[SDL_SCANCODE_Q]) {
					moveDir -= camera.up;
				}
				if(moveDir.MagnitudeSqrd() != 0)
					camera.position = camera.position + moveDir.Normalized() * speed;

				time += delta;
				SDL_UpdateWindowSurface(window);
				auto t2 = std::chrono::high_resolution_clock::now();
				frametime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
				delta = frametime / 1000.0f;
				fps = 1000.0f / frametime;
				cout << "frametime: " << frametime << "\n";
				cout << "delta: " << delta << "\n";
				cout << "fps: " << fps << "\n";;
				t1 = std::chrono::high_resolution_clock::now();


				for (size_t i = 0; i < THREADS; i++)
				{
					waitWorkerFlags[i] = false;
				}
			}

			//if (finished1 && finished2 && finished3 && finished4) {
			//	SDL_UpdateWindowSurface(window);
			//	auto t2 = std::chrono::high_resolution_clock::now();
			//	frametime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
			//	delta = frametime / 1000.0f;
			//	fps = 1000.0f / frametime;
			//	cout << "frametime: " << frametime << "\n";
			//	cout << "delta: " << delta << "\n";
			//	cout << "fps: " << fps << "\n";;
			//	t1 = std::chrono::high_resolution_clock::now();
			//	finished1 = false;
			//	finished2 = false;
			//	finished3 = false;
			//	finished4 = false;
			//}
			//else {
			//	delta = 0;
			//}

#else
#ifdef NOISYRENDER

			for (size_t i = 0; i < raytraceSamples; i++)
			{
				int randX = rand() % SCREEN_WIDTH;
				int randY = rand() % SCREEN_HEIGHT;
				float3 rayDirection = GetRayDirection(camera, randX, randY);
				auto color = Object::RaytraceScene(camera.position, rayDirection);
				//set_pixel(randX, randY, color, MAINSCREEN);

				Uint32* targetPixel = (Uint32*)((Uint8*)screen->pixels + (SCREEN_HEIGHT - 1 - randY) * screen->pitch + randX * screen->format->BytesPerPixel);
				*targetPixel = color;
			}

#else
			for (size_t i = 0; i < SCREEN_WIDTH; i++)
			{
				for (size_t j = 0; j < SCREEN_HEIGHT; j++)
				{
					float3 rayDirection = GetRayDirection(camera, i, j);
					auto color = Object::RaytraceScene((camera).position, rayDirection);
					set_pixel(i, j, color, screen);
				}
			}

#endif // NOISYRENDER ON


			SDL_UpdateWindowSurface(window);

			auto t2 = std::chrono::high_resolution_clock::now();
			frametime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
			delta = frametime / 1000.0f;
			fps = 1000.0f / frametime;
			cout << "frametime: " << frametime << "\n";
			cout << "delta: " << delta << "\n";
			cout << "fps: " << fps << "\n";;
#endif //  MULTITHREAD
		}
	}


#ifdef MULTITHREAD
	suspendAllThreads = true;

	for (size_t i = 0; i < THREADS; i++)
	{
		renderWorkers[i]->join();
		delete renderWorkers[i];

		delete randFunctions[i];
	}
	delete[] renderWorkers;
	delete[] waitWorkerFlags;
	delete[] randFunctions;

	//q1->join();
	//q2->join();
	//q3->join();
	//q4->join();
	//delete q1;
	//delete q2;
	//delete q3;
	//delete q4;
#endif


	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}