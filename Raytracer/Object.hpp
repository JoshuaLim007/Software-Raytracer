#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <vector>
#include "Common.hpp"
#include "json.hpp"
#include "imgui.h"


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
		ImGui::InputFloat3("Position", imGuiFloat3);
		transform.position.x = imGuiFloat3[0];
		transform.position.y = imGuiFloat3[1];
		transform.position.z = imGuiFloat3[2];
		imGuiFloat3[0] = material.BaseColor.r;
		imGuiFloat3[1] = material.BaseColor.g;
		imGuiFloat3[2] = material.BaseColor.b;
		ImGui::ColorPicker3("Color", imGuiFloat3);
		material.BaseColor.r = imGuiFloat3[0];
		material.BaseColor.g = imGuiFloat3[1];
		material.BaseColor.b = imGuiFloat3[2];
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

#endif