#ifndef PBR_H
#define PBR_H

#include <math.h>
#include "Common.hpp"
#include <limits>
#include <Windows.h>

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

#endif