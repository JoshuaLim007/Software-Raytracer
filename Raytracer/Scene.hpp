#ifndef SCENE_HPP
#define SCENE_HPP

#include <vector>
#include "Common.hpp"
#include "json.hpp"
#include <fstream>
#include <exception>

using json = nlohmann::json;

class Scene sealed {
private:
	std::string fileName;
	std::vector<Object*> sceneObjects;
public:
	std::string sceneName;
	Scene(std::string file) {
		fileName = file;
	}
	std::vector<Object*> GetObjects() {
		return sceneObjects;
	}
	std::string GetFilePath() {
		return fileName;
	}
	void Load() {
		sceneObjects.clear();
		std::ifstream f(fileName);
		if (!f.good()) {
			return;
		}
		try {
			json data = json::parse(f);
			sceneName = data["SceneName"];
			json sceneObject = data["SceneObjects"];
			
			for (auto i = sceneObject.begin(); i < sceneObject.end(); i++) {
				json value = i.value();
				json position = value["Position"];

				Object* obj;
				if (value["Renderer"]["Type"] == "Sphere") {
					obj = new Sphere(value["Renderer"]["Radius"], float3(position[0], position[1], position[2]));
				}
				else {
					obj = new Object();
					obj->transform.position = float3(position[0], position[1], position[2]);
				}

				if (value.contains("Material")) {
					json material = value["Material"];
					obj->material.Smoothness = material.contains("Smoothness") ? (float)value["Material"]["Smoothness"] : 0.5f;
					obj->material.SpecularAmount = material.contains("SpecularAmount") ? (float)value["Material"]["SpecularAmount"] : 0.1f;
					json color = material.contains("SpecularColor") ? value["Material"]["SpecularColor"] : json::array({1,1,1});
					obj->material.SpecularColor = Color(color[0], color[1], color[2]);
					color = material.contains("Color") ? value["Material"]["Color"] : json::array({ 1,1,1 });
					obj->material.BaseColor = Color(color[0], color[1], color[2]);
					color = material.contains("Emissive") ? value["Material"]["Emissive"] : json::array({ 0,0,0 });
					obj->material.EmissiveColor = Color(color[0], color[1], color[2]);
				}

				obj->name = value["Name"];
				sceneObjects.push_back(obj);
			}
		}
		catch (std::exception e) {
			std::cout << e.what()  << std::endl;
		}
		f.close();
		return;
	}
	void Unload() {
		for (size_t i = 0; i < sceneObjects.size(); i++)
		{
			delete sceneObjects[i];
		}
		sceneObjects.clear();
	}
	void Save() {
		json data;
		data["SceneName"] = sceneName;
		data["SceneObjects"] = json::array();
		for (size_t i = 0; i < sceneObjects.size(); i++)
		{
			data["SceneObjects"][i] = sceneObjects[i]->ToJSON();
		}
		std::string res = data.dump(4);
		std::ofstream f(fileName);
		f << res;
		f.close();
	}
	void SaveAs(std::string fileName) {
		this->fileName = fileName;
		Save();
	}
	void AddObject(Object* object) {
		sceneObjects.push_back(object);
	}
	void RemoveObject(Object* object) {
		for (auto i = sceneObjects.begin(); i < sceneObjects.end(); i++)
		{
			if ((*i) == object) {
				sceneObjects.erase(i);
			}
		}
	}
	~Scene() {
		Unload();
	}
};

#endif