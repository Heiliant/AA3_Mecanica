#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm/glm.hpp>
#include <iostream>

#define NPARTICLES 252
const float quadScale = 0.5f;
using v3 = glm::vec3;
float ke = 1000;
float kd = 17;
float OD = 0.35f;
float gravity = -9.85f;
float elasticityCof = 0.5;
float frictionCof = 0;
#define MASS 1.f//F
struct particle {
	v3 P;
	v3 Po;
	v3 V;
	v3 Vo;
	v3 F;
};

//PLANOS
const glm::vec3 planes[6] = {
	{ 0, 1, 0 },
	{ 0, -1, 0 },
	{ 1, 0, 0 },
	{ -1, 0, 0 },
	{ 0, 0, 1 },
	{ 0, 0, -1 }
};

const glm::vec3 planePoint[6]{
	{ 0, 0, 0 },
	{ 0, 10, 0 },
	{ -5, 5, 0 },
	{ 5, 5, 0 },
	{ 0, 5, -5 },
	{ 0, 5, 5 }
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

v3 spring(particle P1, particle P2, float originalD=OD) {
	v3 Vector12 = P1.P - P2.P;
	return -(ke*(glm::length(Vector12) - originalD) + kd*glm::dot((P1.V - P2.V), glm::normalize(Vector12)))*glm::normalize(Vector12);
}

//COLISIONES PLANO
float planeD(glm::vec3 normal, glm::vec3 puntoP) {
	return -glm::dot(normal, puntoP);
}

bool hasCollided(particle particula, glm::vec3 normal, float d) {
	return ((glm::dot(normal, particula.Po) + d)*(glm::dot(normal, particula.P) + d) <= 0);
}

void rebote(particle &particula, glm::vec3 normal, glm::vec3 planeSpot) {
	float d = planeD(normal, planeSpot);

	particula.P = particula.P - (1+elasticityCof) * (glm::dot(normal, particula.P) + d)*normal;

	particula.V = particula.V - (1+elasticityCof) * glm::dot(normal, particula.V)*normal;

	glm::vec3 velocidadNormal = (normal*particula.Vo)*normal;
	glm::vec3 velocidadTangencial = particula.Vo - velocidadNormal;

	particula.V -= frictionCof*velocidadTangencial;
}


void FuncionUpdate(float dt) {
	for (int i = 1; i < NPARTICLES; ++i) {
		v3 aux = arrayStructParticles[i].P;
		if (i != 13) {
			for (int j = 0; j < 3; ++j) {
				arrayStructParticles[i].F += v3{ 0, MASS*(gravity), 0 };
				switch (j) {
				case 0:
					arrayStructParticles[i].P.x += (arrayStructParticles[i].P.x - arrayStructParticles[i].Po.x) + (arrayStructParticles[i].F.x / MASS)*glm::pow(dt, 2);
					arrayParticles[i * 3 + j] = arrayStructParticles[i].P.x;
					arrayStructParticles[i].V.x = (arrayStructParticles[i].P.x - arrayStructParticles[i].Po.x) / dt;
					break;
				case 1:
					arrayStructParticles[i].P.y += (arrayStructParticles[i].P.y - arrayStructParticles[i].Po.y) + (arrayStructParticles[i].F.y / MASS)*glm::pow(dt, 2);
					arrayParticles[i * 3 + j] = arrayStructParticles[i].P.y;
					arrayStructParticles[i].V.y = (arrayStructParticles[i].P.y - arrayStructParticles[i].Po.y) / dt;
					break;
				case 2:
					arrayStructParticles[i].P.z += (arrayStructParticles[i].P.z - arrayStructParticles[i].Po.z) + (arrayStructParticles[i].F.z / MASS)*glm::pow(dt, 2);
					arrayParticles[i * 3 + j] = arrayStructParticles[i].P.z;
					arrayStructParticles[i].V.z = (arrayStructParticles[i].P.z - arrayStructParticles[i].Po.z) / dt;
					break;
				}
			}
			arrayStructParticles[i].Po = aux;
		}
		//Colisiones
		int aux2 = 0;
		for (glm::vec3 a : planes) {
			if (hasCollided(arrayStructParticles[i], a, planeD(a, planePoint[aux2]))) {
				rebote(arrayStructParticles[i], a, planePoint[aux2]);
			}
			++aux2;
		}
	}

}

void Structural() {
	////////////////////////////////////////////////////STRUCTURAL////////////////////////////////////////////////////
	for (int i = 0; i < NPARTICLES; ++i) {
		if (i < 14) {
			if (i == 0) {
				arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i + 1]) + spring(arrayStructParticles[i], arrayStructParticles[i + 14]);
			}
			else if (i == 13) {
				arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i - 1]) + spring(arrayStructParticles[i], arrayStructParticles[i + 14]);
			}
			else {
				arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i + 14]) + spring(arrayStructParticles[i], arrayStructParticles[i + 1])
					+ spring(arrayStructParticles[i], arrayStructParticles[i - 1]);
			}
		}
		else if (i >= (NPARTICLES - 14)) {
			if (i == NPARTICLES - 14) {
				arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i + 1]) + spring(arrayStructParticles[i], arrayStructParticles[i - 14]);
			}
			else if (i == NPARTICLES - 1) {
				arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i - 1]) + spring(arrayStructParticles[i], arrayStructParticles[i - 14]);
			}
			else {
				arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i - 14]) + spring(arrayStructParticles[i], arrayStructParticles[i + 1])
					+ spring(arrayStructParticles[i], arrayStructParticles[i - 1]);
			}
		}
		else if (i % 14 == 0) {
			arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i - 14]) + spring(arrayStructParticles[i], arrayStructParticles[i + 1])
				+ spring(arrayStructParticles[i], arrayStructParticles[i + 14]);
		}
		else if (i % 14 == 13) {
			arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i - 14]) + spring(arrayStructParticles[i], arrayStructParticles[i - 1])
				+ spring(arrayStructParticles[i], arrayStructParticles[i + 14]);
		}
		else {
			arrayStructParticles[i].F = spring(arrayStructParticles[i], arrayStructParticles[i - 14]) + spring(arrayStructParticles[i], arrayStructParticles[i + 1])
				+ spring(arrayStructParticles[i], arrayStructParticles[i + 14]) + spring(arrayStructParticles[i], arrayStructParticles[i - 1]);
		}
	}
}

void Shear() {
	////////////////////////////////////////////////////SHEAR////////////////////////////////////////////////////
	for (int i = 0; i < NPARTICLES; ++i) {
		if (i < 14) {
			if (i == 0) {
				arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
			}
			else if (i == 13) {
				arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
			}
			else {
				arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
					+ spring(arrayStructParticles[i], arrayStructParticles[i + 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
			}
		}
		else if (i >= (NPARTICLES - 14)) {
			if (i == NPARTICLES - 14) {
				arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i - 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
			}
			else if (i == NPARTICLES - 1) {
				arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i - 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
			}
			else {
				arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i - 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
					+ spring(arrayStructParticles[i], arrayStructParticles[i - 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
			}
		}
		else if (i % 14 == 0) {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i - 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
				+ spring(arrayStructParticles[i], arrayStructParticles[i + 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
		}
		else if (i % 14 == 13) {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
				+ spring(arrayStructParticles[i], arrayStructParticles[i - 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
		}
		else {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
				+ spring(arrayStructParticles[i], arrayStructParticles[i - 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
				+ spring(arrayStructParticles[i], arrayStructParticles[i - 13], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)))
				+ spring(arrayStructParticles[i], arrayStructParticles[i + 15], glm::sqrt(glm::pow(OD, 2) + glm::pow(OD, 2)));
		}
	}
}

void Bending() {
	////////////////////////////////////////////////////BENDING////////////////////////////////////////////////////
	for (int i = 0; i < NPARTICLES; ++i) {
		if (i % 14 == 0 || i % 14 == 1) {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 2], 2 * OD);
		}
		else if (i % 14 == 13 || i % 14 == 12) {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i - 2], 2 * OD);
		}
		else {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 2], 2 * OD) + spring(arrayStructParticles[i], arrayStructParticles[i - 2], 2 * OD);
		}
	}

	for (int i = 0; i < NPARTICLES; ++i) {
		if (static_cast<int>(i / 14) == 0 || static_cast<int>(i / 14) == 1) {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 28], 2 * OD);
		}
		else if (static_cast<int>(i / 14) == 17 || static_cast<int>(i / 14) == 16) {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i - 28], 2 * OD);
		}
		else {
			arrayStructParticles[i].F += spring(arrayStructParticles[i], arrayStructParticles[i + 28], 2 * OD) + spring(arrayStructParticles[i], arrayStructParticles[i - 28], 2 * OD);
		}
		if (i % 14 == 13) {
			std::cout << std::endl;
		}
	}
}


void PhysicsInit() {
	int localZ = -5;
	for (int i = 0; i < NPARTICLES; ++i) {
		for (int j = 0; j < 3; ++j) {
			switch (j) {
			case 0:	arrayParticles[i * 3 + j] = -OD +(i%14)*OD;
				arrayStructParticles[i].P.x = -OD + (i % 14)*OD;
				if (i % 14 == 0)
					++localZ;
				break;
			case 1:arrayParticles[i * 3 + j] = 5;
				arrayStructParticles[i].P.y = 5;
				break;
			case 2:arrayParticles[i * 3 + j] = -OD +(localZ)*OD;
				arrayStructParticles[i].P.z = -OD + (localZ)*OD;;
				break;
			}
		}
		arrayStructParticles[i].Po = arrayStructParticles[i].P;
		arrayStructParticles[i].F = { 0, 0, 0 };
	}
	ClothMesh::updateClothMesh(arrayParticles);
}

void PhysicsUpdate(float dt) {
	for (int i = 0; i < 7; ++i) {
		FuncionUpdate(dt/7);
		Structural();
		Shear();
		Bending();

	}

	ClothMesh::updateClothMesh(arrayParticles);
}

void PhysicsCleanup() {
	// Do your cleanup code here...
	// ............................
}