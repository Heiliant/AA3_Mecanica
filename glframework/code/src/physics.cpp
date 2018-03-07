#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm/glm.hpp>

#define NPARTICLES 252
const float quadScale = 0.5f;
using v3 = glm::vec3;
#define MASS 1.F//F
struct particle {
	v3 P;
	v3 Po;

};


float arrayParticles[NPARTICLES * 3];
particle arrayStructParticles[NPARTICLES];

namespace ClothMesh {
	extern void setupClothMesh();
	extern void cleanupClothMesh();
	extern void updateClothMesh(float* array_data);
	extern void drawClothMesh();
}

bool show_test_window = false;
void GUI() {
	bool show = true;
	ImGui::Begin("Physics Parameters", &show, 0);

	// Do your GUI code here....
	{	
		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);//FrameRate
		
	}
	// .........................
	
	ImGui::End();

	// Example code -- ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
	if(show_test_window) {
		ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
		ImGui::ShowTestWindow(&show_test_window);
	}
}

void PhysicsInit() {
	int localZ = -5;
	for (int i = 0; i < NPARTICLES; ++i) {
		for (int j = 0; j < 3; ++j) {
			switch (j) {
			case 0:	arrayParticles[i * 3 + j] = -3.f+(i%14)*quadScale;
				arrayStructParticles[i].P.x = -3.f + (i % 14)*quadScale;
				if (i % 14 == 0)
					++localZ;
				break;
			case 1:arrayParticles[i * 3 + j] = 9;
				arrayStructParticles[i].P.y = 9;
				break;
			case 2:arrayParticles[i * 3 + j] = -2.f+(localZ)*quadScale;
				arrayStructParticles[i].P.z = -2.f + (localZ)*quadScale;;
				break;
			}
		}
		arrayStructParticles[i].Po = arrayStructParticles[i].P;
	}
	ClothMesh::updateClothMesh(arrayParticles);
}

void PhysicsUpdate(float dt) {

	for (int i = 1; i < NPARTICLES; ++i) {
		v3 aux = arrayStructParticles[i].P;
		if (i!=13){
		for (int j = 0; j < 3; ++j) {
			v3 sumFuerzas = { 0, MASS*(-9.81), 0 };
			switch (j) {
			case 0:
				arrayStructParticles[i].P.x += (arrayStructParticles[i].P.x - arrayStructParticles[i].Po.x) + (sumFuerzas.x / MASS)*glm::pow(dt, 2);
				arrayParticles[i * 3 + j] = arrayStructParticles[i].P.x;
				break;
			case 1:
				arrayStructParticles[i].P.y += (arrayStructParticles[i].P.y - arrayStructParticles[i].Po.y) + (sumFuerzas.y / MASS)*glm::pow(dt, 2);
				arrayParticles[i * 3 + j] = arrayStructParticles[i].P.y;
				break;
			case 2:
				arrayStructParticles[i].P.z += (arrayStructParticles[i].P.z - arrayStructParticles[i].Po.z) + (sumFuerzas.z / MASS)*glm::pow(dt, 2);
				arrayParticles[i * 3 + j] = arrayStructParticles[i].P.z;
				break;
			}
		}
		arrayStructParticles[i].Po = aux;
	}
	}
	ClothMesh::updateClothMesh(arrayParticles);
}

void PhysicsCleanup() {
	// Do your cleanup code here...
	// ............................
}