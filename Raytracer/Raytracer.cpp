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


using std::cout; using std::endl;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;

const int THREADS = 16;
#define WORDLRIGHT float3(1,0,0)
#define WORLDUP float3(0,1,0)
#define WORLDFORWARD float3(0,0,1)
#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define SCREEN_SCALE .5
#define FOV 55
#define MAXBOUNCES 2
int ACCUMULATIONFRAMES = 8;

int frameCount = 0;

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
	static float3 Cross(const float3& lhs, const float3& rhs) {
		return float3(lhs.y * rhs.z - rhs.y * lhs.z, rhs.x * lhs.z - lhs.x * rhs.z, lhs.x * rhs.y - rhs.x * lhs.y);
	}
	static float3 Lerp(const float3& a, const float3& b, float t){
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
	bool operator != (const float3& rhs) {
		return x != rhs.x || y != rhs.y || z != rhs.z;
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
	Color operator /(const Color& rhs) const {
		return Color(r / rhs.r, g / rhs.g, b / rhs.b, a / rhs.a);
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
struct Material {
public:
	float Smoothness;
	float Metalness;
	Color BaseColor;
	Color EmissiveColor;
	Material(Color BaseColor, Color Emission, float smoothness, float metalness) {
		this->BaseColor = BaseColor;
		EmissiveColor = Emission;
		Smoothness = smoothness;
		Metalness = metalness;
	}
	Material() {
		Smoothness = 0.5f;
		Metalness = 0.0f;
		BaseColor = Color(1, 1, 1);
		EmissiveColor = Color(0, 0, 0);
	}
};
struct Rayhit {
	bool valid = false;
	Material material;
	float3 normal;
	float3 point;
	float distance;
};

SDL_Surface* renderTarget;
float3 SunDirection = float3(1, -1, -1);
Color SkyColor = Color(.2, .35, 1.0f) * 10.0f;
Color HorizonColor = Color(1.0, 0.9f, 0.5f) * 5.0f;
Color GroundColor = Color(.08f, .06f, .03f);
Color SunColor = Color(500, 500, 500);
Color colorBuffer[SCREEN_HEIGHT * SCREEN_WIDTH];
bool setFrame = 0;

void SetScreenPixel(const int& x, const int& y, const Color& color) {
	Uint32* targetPixel = (Uint32*)((Uint8*)renderTarget->pixels + (SCREEN_HEIGHT - 1 - y) * renderTarget->pitch + x * renderTarget->format->BytesPerPixel);
	float weight = 1.0 / (ACCUMULATIONFRAMES);
	if (setFrame > 0) {
		colorBuffer[x + y * SCREEN_WIDTH] = colorBuffer[x + y * SCREEN_WIDTH] * (1 - weight) + color * weight;
	}
	else {
		ACCUMULATIONFRAMES = 8;
		colorBuffer[x + y * SCREEN_WIDTH] = color;
	}
	if (frameCount > ACCUMULATIONFRAMES) {
		ACCUMULATIONFRAMES++;
	}
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

struct PBRModel {
public:
	//PBR Reflectance model copied from https://www.youtube.com/watch?v=RRE-F57fbXw&t=425s&ab_channel=VictorGordan
	static auto NormaliDistribution(float roughness, const float3& halfVector, const float3& normal) {
		//GGX/Trowbridge-Reitz
		float3 h = halfVector;
		float rr = roughness * roughness;
		float aa = rr * rr;
		float NdotH = float3::Dot(normal, h, true);
		float denom = max(std::_Pi * powf((NdotH * NdotH) * (aa - 1) + 1, 2), FLT_EPSILON);
		return aa / denom;
	};
	static auto GeometricAttenuation(float roughness, const float3& normal, const float3& xDirection) {
		//Schlick-Beckmann model
		float a = roughness * roughness;
		float k = a / 2;
		float denom = max((float3::Dot(normal, xDirection, true) * (1 - k) + k), FLT_EPSILON);
		auto f = float3::Dot(normal, xDirection, true) / denom;
		return f;
	};
	static float Fresnal(const float& baseReflectivity, const float3& viewDirection, const float3& H) {
		//Schlicks fresnal approximation
		float VdotH = float3::Dot(viewDirection, H, true);
		return (baseReflectivity + ((1 - baseReflectivity) * powf(1 - VdotH, 5)));
	};
};

class Object {
private:
	int myIndex;
	const float3 WorldUp = float3(0,1,0);
public:
	static std::vector<Object*> allObjects;
	virtual Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const {
		return Rayhit();
	}
	Transform transform;
	Material material;
	Object() {
		allObjects.push_back(this);
		myIndex = allObjects.size() - 1;
	}
	virtual ~Object() {
		allObjects.erase(allObjects.begin() + myIndex);
	}
};
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
		Material& outColor,
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

		outColor = material;
		return true;
	}
public:
	Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const override {

		Material color;
		float3 normal, hitPoint;
		float distance;
		Rayhit hitResults;

		if (line_sphere_intersection(rayOrigin, rayDir, transform.position, radius, color, normal, hitPoint, distance)) {
			hitResults.material = color;
			hitResults.normal = normal;
			hitResults.point = hitPoint;
			hitResults.distance = distance;
			hitResults.valid = true;
		}
		return hitResults;
	}
};

float3 GetRandomNormalOrientedHemisphere(const float3& normal) {
	float3 sr;
	do {
		sr.x = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr.y = ((float)rand() / RAND_MAX - 0.5f) * 2;
		sr.z = ((float)rand() / RAND_MAX - 0.5f) * 2;
	} while (sr.Magnitude() > 1);
	if (float3::Dot(sr, normal) < 0) {
		sr = sr * -1;
	}
	return sr.Normalized();
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
Rayhit GetClosestObject(const float3& rayOrigin, const float3& rayDirection) {
	int objectCount = Object::allObjects.size();
	Rayhit returnResults;
	float shortestDistance = std::numeric_limits<float>().infinity();
	for (size_t i = 0; i < objectCount; i++)
	{
		Rayhit hitResults = Object::allObjects[i]->Raytrace(rayOrigin, rayDirection);
		if (hitResults.valid) {

			if (hitResults.distance < shortestDistance) {
				shortestDistance = hitResults.distance;
				returnResults = hitResults;
			}
		}
	}
	return returnResults;
}
Color RaytraceScene(const float3& rayOrigin, const float3& rayDirection) {
	auto hit = GetClosestObject(rayOrigin, rayDirection);
	if (!hit.valid) {
		return GetEnvironmentColor(rayDirection);
	}

	Color specularColor;
	Color reflectionColor;

	Color incomingLight = hit.material.EmissiveColor;
	Color hitColor = hit.material.BaseColor;
	
	for (size_t i = 0; i < MAXBOUNCES; i++)
	{
		auto sray = GetRandomNormalOrientedHemisphere(hit.normal);
		hit = GetClosestObject(hit.point + hit.normal * .00001f, sray);
		if (!hit.valid) {
			incomingLight += GetEnvironmentColor(sray) * hitColor;
			break;
		}
		incomingLight += hit.material.EmissiveColor * hitColor;
		float lightStrength = float3::Dot(hit.normal, sray, true);
		hitColor *= hit.material.BaseColor * lightStrength;
	}
	Color diffusedColor = incomingLight;

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

	return diffusedColor;
}

#pragma region ThreadUniqueRandFunctions
std::uniform_int_distribution<int>* dice_distribution = new std::uniform_int_distribution<int>(0, INT_MAX);
bool suspendAllThreads = false;
auto renderWorkers = new std::thread * [THREADS];
auto waitWorkerFlags = new bool[THREADS];
auto randFunctions = new std::_Binder<std::remove_cv<std::_Unforced>::type, std::uniform_int_distribution<int>&, std::mt19937&>*[THREADS];
#pragma endregion

void renderArea(
	unsigned int index,
	unsigned int minX, unsigned int maxX, 
	unsigned int minY, unsigned int maxY, 
	const Transform* camera) {

	int steps = 1 / (float)SCREEN_SCALE;

	//min inclusive
	//max exclusive
	while (!suspendAllThreads) {
		if (waitWorkerFlags[index] == false) {
			Color color;
			for (size_t i = minX; i < maxX; i += steps)
			{
				for (size_t j = minY; j < maxY; j += steps)
				{
					float3 rayDirection = GetRayDirection(*camera, i, j);
					color = RaytraceScene((*camera).position, rayDirection);

					for (size_t i1 = 0; i1 < steps; i1++)
					{
						for (size_t j1 = 0; j1 < steps; j1++)
						{
							SetScreenPixel(i + i1, j + j1, color);
						}
					}
				}
			}
			waitWorkerFlags[index] = true;
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
	SDL_Window* window = SDL_CreateWindow("Raytracer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_VULKAN);
	renderTarget = SDL_GetWindowSurface(window);
	SDL_FillRect(renderTarget, 0, Color(0, 0, 0, 0));
	auto cursor = SDL_GetDefaultCursor();
	SDL_SetCursor(cursor);
	long long frametime = 0;
	float fps = 0;
	SDL_SetRelativeMouseMode(SDL_TRUE);
#pragma endregion

#pragma region SceneSetup
	Transform camera;
	camera.position = float3(0, 0, 0);
	camera.RotateAboutAxis(0, WORDLRIGHT);
	Object* d;
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
			), Color(0, 0, 0), 0.5, 0);
		}
	}
	d = new Sphere(1, float3(0, 0, 5));
	d->material = Material(Color(1, 1, 1), Color(0, 0, 0), 0.5, 0);
	//d = new Sphere(2, float3(4, 4, 8));
	//d->material = Material(Color(1, 1, 1), Color(50, 50, 50), 0.5, 0);
	d = new Sphere(1000, float3(0, -1001.2f, 5));
	d->material = Material(Color(1, 1, 1), Color(0, 0, 0), 0.5, 0);

#pragma endregion

#pragma region ThreadSetup
	int div = round(SCREEN_WIDTH / (THREADS));
	for (size_t i = 0; i < THREADS; i++)
	{
		waitWorkerFlags[i] = false;
		auto mt1 = std::mt19937();
		mt1.seed(i);
		auto d = std::bind(*dice_distribution, mt1);

		randFunctions[i] = new std::_Binder<std::remove_cv<std::_Unforced>::type, std::uniform_int_distribution<int>&, std::mt19937&>(d);
		int initialX = div * i;
		renderWorkers[i] = new std::thread(renderArea, i, initialX, initialX + div, 0, SCREEN_HEIGHT, &camera);
	}
#pragma endregion

	bool pause = false;
	bool doSetFrame = false;
	float delta = 0;
	int mouseX = 0;
	int mouseY = 0;
	float mouseSpeed = .05;
	float moveSpeed = 1;

	std::chrono::steady_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	while (true) {
		SDL_PumpEvents();
		const Uint8* state = SDL_GetKeyboardState(NULL);
		if (state[SDL_SCANCODE_ESCAPE]) {
			break;
		}

		if (pause) {
			continue;
		}

		SDL_GetRelativeMouseState(&mouseX, &mouseY);
		if (mouseX != 0 || mouseY != 0) {
			doSetFrame = true;
		}
		camera.RotateAboutAxis(mouseX * mouseSpeed * delta, WORLDUP);
		camera.RotateAboutAxis(mouseY * mouseSpeed * delta, camera.right);

		bool falseCheck = false;
		for (size_t i = 0; i < THREADS; i++)
		{
			if (waitWorkerFlags[i] == false) {
				falseCheck = true;
				break;
			}
		}

		if (falseCheck) {
			continue;
		}

		float speed = moveSpeed * delta;
		float3 previousPos = camera.position;
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
		if (camera.position != previousPos) {
			doSetFrame = true;
		}


		SDL_UpdateWindowSurface(window);
		auto t2 = std::chrono::high_resolution_clock::now();
		frametime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		delta = frametime / 1000.0f;
		fps = 1000.0f / frametime;
		cout << "frametime: " << frametime << "\n";
		cout << "delta: " << delta << "\n";
		cout << "fps: " << fps << "\n";
		cout << "ACCUMULATIONFRAMES: " << ACCUMULATIONFRAMES << "\n";
		t1 = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < THREADS; i++)
		{
			waitWorkerFlags[i] = false;
		}
		if (doSetFrame) {
			setFrame = 0;
			frameCount = 0;
			doSetFrame = false;
		}
		else {
			frameCount++;
			setFrame = 1;
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
	delete[] waitWorkerFlags;
	delete[] randFunctions;

	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}