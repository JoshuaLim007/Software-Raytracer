#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <vector>
#include "Common.hpp"
#include "json.hpp"
#include "imgui.h"
#include <math.h>

using json = nlohmann::json;
class Serializable abstract{
public:
	virtual json ToJSON() = 0;
};
class SerializeGUI abstract {
public:
	virtual void OnGUI() = 0;
};
class Object : public Serializable, public SerializeGUI {
public:
	virtual Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const {
		return Rayhit();
	}
	Transform transform;
	Material material;
	std::string name;
	json ToJSON() override {
		json data;
		data["Name"] = name;
		auto position = transform.position;
		data["Position"] = json::array({ position.x, position.y, position.z });
		data["Material"]["Smoothness"] = material.Smoothness;
		data["Material"]["Metalness"] = material.SpecularAmount;
		auto color = material.BaseColor;
		data["Material"]["Color"] = json::array({ color.r, color.g, color.b });
		color = material.EmissiveColor;
		data["Material"]["Emissive"] = json::array({ color.r, color.g, color.b });
		color = material.SpecularColor;
		data["Material"]["SpecularColor"] = json::array({ color.r, color.g, color.b });
		data["Material"]["SpecularAmount"] = material.SpecularAmount;
		data["Renderer"] = { { "Type", "None" } };
		return data;
	}
	void OnGUI() override {
		char nameBuffer[256];
		strcpy_s(nameBuffer, name.c_str());
		ImGui::InputText("Name", nameBuffer, 256);
		name = std::string(nameBuffer);
		float imGuiFloat3[3];
		imGuiFloat3[0] = transform.position.x;
		imGuiFloat3[1] = transform.position.y;
		imGuiFloat3[2] = transform.position.z;
		ImGui::DragFloat3("Position", imGuiFloat3, 0.1f);
		transform.position.x = imGuiFloat3[0];
		transform.position.y = imGuiFloat3[1];
		transform.position.z = imGuiFloat3[2];

		if (ImGui::CollapsingHeader("Base Color")) {
			imGuiFloat3[0] = material.BaseColor.r;
			imGuiFloat3[1] = material.BaseColor.g;
			imGuiFloat3[2] = material.BaseColor.b;
			ImGui::ColorPicker3("Color", imGuiFloat3);
			material.BaseColor.r = imGuiFloat3[0];
			material.BaseColor.g = imGuiFloat3[1];
			material.BaseColor.b = imGuiFloat3[2];
		}

		imGuiFloat3[0] = material.EmissiveColor.r;
		imGuiFloat3[1] = material.EmissiveColor.g;
		imGuiFloat3[2] = material.EmissiveColor.b;
		ImGui::InputFloat3("Emissive Color", imGuiFloat3);
		material.EmissiveColor.r = imGuiFloat3[0];
		material.EmissiveColor.g = imGuiFloat3[1];
		material.EmissiveColor.b = imGuiFloat3[2];
		ImGui::SliderFloat("Smoothness", &material.Smoothness, 0, 1);
		ImGui::SliderFloat("Specular Amount", &material.SpecularAmount, 0, 1);
		ImGui::NewLine();
	}
};

struct RayHitObject {
	Rayhit rayHit;
	Object* objectReference;
};

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
	json ToJSON() override{
		json data = Object::ToJSON();
		data["Renderer"] = { { "Type", "Sphere" }, {"Radius", radius}};
		return data;
	}
	void OnGUI() override {
		Object::OnGUI();
		ImGui::InputFloat("Sphere Radius", &radius);
		ImGui::NewLine();
	}
	Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const override {

		Material color;
		float3 normal, hitPoint;
		float distance;
		Rayhit hitResults;

		if (line_sphere_intersection(rayOrigin, rayDir, transform.position, radius, color, normal, hitPoint, distance)) {
			hitResults.normal = normal;
			hitResults.point = hitPoint;
			hitResults.distance = distance;
			hitResults.valid = true;
		}
		return hitResults;
	}
};

class Box : public Object {
private:
	// Box: https://www.shadertoy.com/view/ld23DV
	float iBox(float3 ro, float3 rd, float3 distBound, float3& normal,
		float3 boxSize) const {
		float3 m = sign(rd) / float3::Max(abs(rd), 1e-8);
		float3 n = m * ro;
		float3 k = abs(m) * boxSize;

		float3 t1 = -n - k;
		float3 t2 = -n + k;
		float tN = max(max(t1.x, t1.y), t1.z);
		float tF = min(min(t2.x, t2.y), t2.z);

		if (tN > tF || tF <= 0.) {
			return std::numeric_limits<float>().max();
		}
		else {
			if (tN >= distBound.x && tN <= distBound.y) {
				normal = -sign(rd) * step(t1.yzx(), t1) * step(t1.zxy(), t1);
				return tN;
			}
			else if (tF >= distBound.x && tF <= distBound.y) {
				normal = -sign(rd) * step(t1.yzx(), t1) * step(t1.zxy(), t1);
				return tF;
			}
			else {
				return std::numeric_limits<float>().max();
			}
		}
	}

public:
	float3 size;
	Box(float3 size) {
		this->size = size;
	}
	void OnGUI() override {
		Object::OnGUI();
		float temp[3];
		temp[0] = size.x;
		temp[1] = size.y;
		temp[2] = size.z;
		ImGui::DragFloat3("Cube Size", temp, 0.1f);
		size.x = temp[0];
		size.y = temp[1];
		size.z = temp[2];
	}
	json ToJSON() override {
		json data = Object::ToJSON();
		data["Renderer"] = { {"Type", "Cube"}, {"Size", json::array({size.x,size.y,size.z})}};
		return data;
	}

	Rayhit Raytrace(const float3& rayOrigin, const float3& rayDir) const override {
		float3 normal;
		float dist = iBox(rayOrigin - transform.position, rayDir, float3(0.01,10000,0), normal, size);
		Rayhit ret;
		ret.normal = normal;
		ret.point = rayOrigin + rayDir * dist;
		ret.distance = dist;
		ret.valid = dist == std::numeric_limits<float>().max() ? false : true;
		return ret;
	}
};

#endif