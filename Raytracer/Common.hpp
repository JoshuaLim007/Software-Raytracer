#ifndef COMMON_H
#define COMMON_H

#define WORDLRIGHT float3(1,0,0)
#define WORLDUP float3(0,1,0)
#define WORLDFORWARD float3(0,0,1)
#define max(a,b) ((((a) > (b)) ? (a) : (b)))
#define min(a,b) ((((a) < (b)) ? (a) : (b)))

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
	float3 operator-() const {
		return float3(x * -1, y * -1, z * -1);
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
	Color& operator *= (const Color& o) {
		*this = *this * o;
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
	Color& operator /= (const Color& o) {
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
	float SpecularAmount;
	Color BaseColor;
	Color EmissiveColor;
	Color SpecularColor = Color(1,1,1);
	Material(Color BaseColor, Color Emission, float smoothness, float metalness) {
		this->BaseColor = BaseColor;
		EmissiveColor = Emission;
		Smoothness = smoothness;
		SpecularAmount = metalness;
	}
	Material(Color BaseColor, Color SpecularColor, Color Emission, float smoothness, float metalness) {
		this->BaseColor = BaseColor;
		this->SpecularColor = SpecularColor;
		EmissiveColor = Emission;
		Smoothness = smoothness;
		SpecularAmount = metalness;
	}
	Material() {
		Smoothness = 0.5f;
		SpecularAmount = 0.0f;
		BaseColor = Color(1, 1, 1);
		EmissiveColor = Color(0, 0, 0);
	}
};
struct Rayhit {
	bool valid = false;
	float3 normal;
	float3 point;
	float distance;
};

#endif