#include <imgui\imgui.h>
#include <imgui\imgui_impl_sdl_gl3.h>
#include <glm/glm.hpp>
#include <iostream>

#define NPARTICLES 252
const float quadScale = 0.5f;
using v3 = glm::vec3;
float ke = 1000;
float kd = 17;
float OD = 0.4f;
float gravity = -9.85f;

namespace GUIvars {
	bool show_test_window = false;
	bool PlaySimulation = true;
	bool show = false;
	float ResetTime = 5;
	float GX = 0.0f;
	float GY = -9.81f;
	float GZ = 0.0f;
	float Gravity[3] = { GX, GY, GZ };
	float stretch_ke = 0;
	float stretch_kd = 0;
	float K_stretch[2] = { stretch_ke, stretch_kd };
	float shear_ke = 0;
	float shear_kd = 0;
	float K_shear[2] = { shear_ke, shear_kd };
	float bend_ke = 0;
	float bend_kd = 0;
	float K_bend[2] = { bend_ke, bend_kd };
	float ParticleLinl = 0;
	bool useCollisions = true;
	bool useSphereCollider = true;
	float eCo = 0.5;
	float fCo = 0.5;


	float SphX = 1.0f;
	float SphY = 1.0f;
	float SphZ = 1.0f;
	float SphPosition[3] = { SphX, SphY, SphZ };
	float SphRadius = 1.f;
	bool useGravity = true;
	float AX = 0.0f;
	float AY = -9.81f;
	float AZ = 0.0f;
	static float Gacceleration[3] = { AX, AY, AZ };
	bool colisionEsfera = true;

	//float SpherePos[3]{ 0, 2, 0 };
	//float SphereRad = 2.5f;
	extern bool renderSphere;
	extern bool renderCapsule;
}

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

namespace Sphere {
	extern void setupSphere(glm::vec3 pos = { GUIvars::SphPosition[0], GUIvars::SphPosition[1], GUIvars::SphPosition[2] }, float radius = GUIvars::SphRadius);
	extern void cleanupSphere();
	extern void updateSphere(glm::vec3 pos = glm::vec3(&GUIvars::SphPosition[0], &GUIvars::SphPosition[1], &GUIvars::SphPosition[2]), float radius = GUIvars::SphRadius);
	extern void drawSphere();
}

bool show_test_window = false;
void GUI() {

	ImGui::Begin("Physics Parameters", &GUIvars::show, 0);

	{
		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);//FrameRate
		ImGui::Checkbox("Play Simulation", &GUIvars::PlaySimulation);//CHECKBOX
		if (ImGui::Button("Reset Simulation")) {
		}
		ImGui::InputFloat("Reset time", &GUIvars::ResetTime);
		ImGui::InputFloat3("Gravity acceleration", GUIvars::Gravity);

		if (ImGui::TreeNode("Spring Parameters")) {     //FORCES
			ImGui::InputFloat2("K_stretch", GUIvars::K_stretch);
			ImGui::InputFloat2("K_shear", GUIvars::K_shear);
			ImGui::InputFloat2("K_bend", GUIvars::K_bend);
			ImGui::InputFloat("Particle linl", &GUIvars::ParticleLinl);
			ImGui::TreePop();
		}

		if (ImGui::TreeNode("Collisions")) {     //FORCES
			ImGui::Checkbox("Use Collisions", &GUIvars::useCollisions);
			ImGui::Checkbox("Use Sphere Collisions", &GUIvars::useSphereCollider);
			ImGui::TreePop();
			ImGui::InputFloat("Elasticity Coef", &GUIvars::eCo);
			ImGui::InputFloat("Friction Coef", &GUIvars::fCo);
		}

		//	if (ImGui::TreeNode("Emitter")) {
		//		ImGui::Text("Emitter Rate");
		//		ImGui::InputInt("Emitter Rate", &GUIvars::emitRate);
		//		//if (PlaySimulation) {		comprobacio de si cambia, spoiler: si que ho fa
		//		//}
		//		ImGui::InputInt("Particle Life", &GUIvars::particleLife);
		//		//if (PlaySimulation) {   comprobacio de si cambia, spoiler: si que ho fa
		//		//}

		//		ImGui::RadioButton("Fountain", &GUIvars::modo, 0); ImGui::SameLine(); ImGui::RadioButton("Cascade", &GUIvars::modo, 1); //FOUNTAIN OR CASCADE
		//		if (GUIvars::modo == 0) {
		//			ImGui::Text("Fountain Pos"); //BOX (VEC3) 
		//			ImGui::InputFloat3("Fountain Pos", GUIvars::Fposicioinicial);  //Tots els parametres de InputFloat3 o Input Int t els canvia desde la UI
		//			if (ImGui::InputFloat3("Fountain Dir", GUIvars::Fvelocitatinicial))
		//				ImGui::InputFloat("Fountain Angle", &GUIvars::FountainAngle);
		//		}
		//		if (GUIvars::modo == 1) {
		//			ImGui::Text("Cascade Pos"); //BOX (VEC3) 
		//			ImGui::InputFloat3("Cascade PosA", GUIvars::CposicioinicialA);
		//			ImGui::InputFloat3("Cascade PosB", GUIvars::CposicioinicialB);
		//			ImGui::InputFloat3("Cascade Velocity", GUIvars::Cvelocitatinicial);
		//		}
		//		ImGui::InputFloat("Fountain Angle", &GUIvars::FountainAngle);
		//		ImGui::TreePop();
		//	}


		//	if (ImGui::TreeNode("Elasticity & Friction")) {			//ELASTICITY &FRICTION
		//		ImGui::InputFloat("Elasticity Coefficient", &GUIvars::elasticityCof);
		//		ImGui::InputFloat("Friction Coefficient", &GUIvars::frictionCof);
		//		ImGui::TreePop();
		//	}

		//	if (ImGui::TreeNode("Colliders")) {						//COLLIDERS
		//		ImGui::Checkbox("Use Sphere Collider", &GUIvars::renderSphere);

		//		if (ImGui::InputFloat3("Sphere Position", GUIvars::SpherePos) || ImGui::InputFloat("Sphere Radius", &GUIvars::SphereRad)) {
		//			Sphere::setupSphere({ GUIvars::SpherePos[0], GUIvars::SpherePos[1], GUIvars::SpherePos[2] }, GUIvars::SphereRad);
		//		}
		//		ImGui::Checkbox("Use Capsule Collider", &GUIvars::renderCapsule);
		//		if (ImGui::InputFloat3("Capsule Position A", GUIvars::CapposicioinicialA) || ImGui::InputFloat3("Capsule Position B", GUIvars::CapposicioinicialB) || ImGui::InputFloat("Capsule Radius", &GUIvars::CapRadius)) {
		//			Capsule::setupCapsule({ GUIvars::CapposicioinicialA[0], GUIvars::CapposicioinicialA[1], GUIvars::CapposicioinicialA[2] },
		//			{ GUIvars::CapposicioinicialB[0], GUIvars::CapposicioinicialB[1] ,GUIvars::CapposicioinicialB[2] }, GUIvars::CapRadius);
		//		}

		//		ImGui::TreePop();
		//	}

		//	if (ImGui::TreeNode("Forces")) {     //FORCES
		//		ImGui::Checkbox("Use Gravity", &GUIvars::useGravity);
		//		ImGui::InputFloat3("Gravity Accel", GUIvars::Gacceleration);
		//		ImGui::TreePop();
		//	}
		//}
		//// .........................

		ImGui::End();
	}
}

v3 spring(particle P1, particle P2, float originalD = OD) {
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
	normal = glm::normalize(normal);
	float d = planeD(normal, planeSpot);

	particula.Po = particula.Po - (1 + GUIvars::eCo) * (glm::dot(normal, particula.Po) + d)*normal;

	particula.P = particula.P - (1 + GUIvars::eCo) * (glm::dot(normal, particula.P) + d)*normal;

	//particula.V = particula.V - (1+ GUIvars::elasticityCof) * glm::dot(normal, particula.V)*normal;

	//glm::vec3 velocidadNormal = (normal*particula.Vo)*normal;
	//glm::vec3 velocidadTangencial = particula.Vo - velocidadNormal;

	//particula.V -= GUIvars::frictionCof*velocidadTangencial;
}


//COLISIONES ESFERA
glm::vec3 colisionSpot(particle particula, glm::vec3 SpherePosition, float Sphradius) {
	float c =
		glm::pow(SpherePosition.x, 2) - 2 * particula.Po.x*SpherePosition.x + glm::pow(particula.Po.x, 2) +
		glm::pow(SpherePosition.y, 2) - 2 * particula.Po.y*SpherePosition.y + glm::pow(particula.Po.y, 2) +
		glm::pow(SpherePosition.z, 2) - 2 * particula.Po.z*SpherePosition.z + glm::pow(particula.Po.z, 2)
		- glm::pow(Sphradius, 2);

	c = glm::pow(glm::pow(particula.P.x - particula.Po.x, 2) + glm::pow(particula.P.y - particula.Po.y, 2) + glm::pow(particula.P.z - particula.Po.z, 2), 2);

	float b =
		2 * (SpherePosition.x*particula.P.x - SpherePosition.x*particula.Po.x - particula.Po.x*particula.P.x + glm::pow(particula.Po.x, 2)) +
		2 * (SpherePosition.y*particula.P.y - SpherePosition.y*particula.Po.y - particula.Po.y*particula.P.y + glm::pow(particula.Po.y, 2)) +
		2 * (SpherePosition.z*particula.P.z - SpherePosition.z*particula.Po.z - particula.Po.z*particula.P.z + glm::pow(particula.Po.z, 2));
	b = 2 * ((particula.P.x - particula.Po.x) + (particula.P.y - particula.Po.y) + (particula.P.z - particula.Po.z)*-1);


	float a =
		glm::pow(particula.P.x, 2) - 2 * particula.P.x*particula.Po.x + glm::pow(particula.Po.x, 2) +
		glm::pow(particula.P.y, 2) - 2 * particula.P.y*particula.Po.y + glm::pow(particula.Po.y, 2) +
		glm::pow(particula.P.z, 2) - 2 * particula.P.z*particula.Po.z + glm::pow(particula.Po.z, 2);

	a = glm::pow((particula.Po.x - SpherePosition.x), 2) + pow((particula.Po.y - SpherePosition.y), 2) + pow((particula.Po.z - SpherePosition.z), 2) -
		glm::pow(Sphradius, 2) + 2 * ((particula.Po.x - SpherePosition.x) + (particula.Po.y - SpherePosition.y) + (particula.Po.z - SpherePosition.z));

	double cuerpo = glm::pow(b, 2) - 4 * a*c;

	if (cuerpo<0) {
		cuerpo *= -1;
	}
	double alphaOne = (-b + glm::sqrt(cuerpo)) / (2 * a);
	double alphaTwo = (-b - glm::sqrt(cuerpo)) / (2 * a);

	if (alphaOne<alphaTwo) {
		return
		{ particula.Po.x + (-particula.P.x + particula.Po.x)*alphaOne,
			particula.Po.y + (-particula.P.y + particula.Po.y)*alphaOne,
			particula.Po.z + (-particula.P.z + particula.Po.z)*alphaOne
		};

	}
	else {
		return
		{ particula.Po.x + (-particula.P.x + particula.Po.x)*alphaTwo,
			particula.Po.y + (-particula.P.y + particula.Po.y)*alphaTwo,
			particula.Po.z + (-particula.P.z + particula.Po.z)*alphaTwo
		};
	}

}

glm::vec3 RSphere(particle particula, float alpha) {
	//Calculo la R que no es R de radio, creo que es R de recta, pero no estoy muy seguro, yo solo sigo lo que hice sobre papel.
	glm::vec3 R;
	R.x = (particula.Po.x + (particula.P.x - particula.Po.x)*alpha);
	R.y = (particula.Po.y + (particula.P.y - particula.Po.y)*alpha);
	R.z = (particula.Po.z + (particula.P.z - particula.Po.z)*alpha);
	return R;
}

glm::vec3 NormalSphere(glm::vec3 R, glm::vec3 SpherePosition) {
	//Calculo la normal del plano del rebote de la esfera.
	glm::vec3 N;
	float modulo;
	N.x = R.x - SpherePosition.x;
	N.y = R.y - SpherePosition.y;
	N.z = R.z - SpherePosition.z;
	modulo = sqrt(pow(N.x, 2) + pow(N.y, 2) + pow(N.z, 2));
	N.x /= modulo;
	N.y /= modulo;
	N.z /= modulo;
	return N;
}

float comprobanteEsfera(glm::vec3 R, glm::vec3 esfera, float alpha) {
	return sqrt(pow((R.x - esfera.x), 2) + pow((R.y - esfera.y), 2) + pow((R.z - esfera.z), 2));
}

bool hasCollidedSphere(particle particula, glm::vec3 SpherePosition, float Sphradius) {
	//Con esta función compruebo si la partícula ha colisionado con la esfera.
	return (glm::distance(particula.P, SpherePosition) - Sphradius) <= 0;
}

void ColisionesCubo(int i) {
	int aux2 = 0;
	for (glm::vec3 a : planes) {
		if (hasCollided(arrayStructParticles[i], a, planeD(a, planePoint[aux2]))) {
			rebote(arrayStructParticles[i], a, planePoint[aux2]);
		}
		++aux2;
	}
}

void ColisionesEsfera(int i) {
	if (GUIvars::renderSphere) {
		if (hasCollidedSphere(arrayStructParticles[i], { GUIvars::SphPosition[0], GUIvars::SphPosition[1], GUIvars::SphPosition[2] }, GUIvars::SphRadius)) {
			glm::vec3 collisionSpot = colisionSpot(arrayStructParticles[i], { GUIvars::SphPosition[0], GUIvars::SphPosition[1], GUIvars::SphPosition[2] }, GUIvars::SphRadius);

			glm::vec3 planeNormal = (collisionSpot - glm::vec3{ GUIvars::SphPosition[0], GUIvars::SphPosition[1], GUIvars::SphPosition[2] });
			rebote(arrayStructParticles[i], planeNormal, collisionSpot);
		}
	}
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
		ColisionesCubo(i);
		ColisionesEsfera(i);
		
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
			case 0:	arrayParticles[i * 3 + j] = -OD + (i % 14)*OD;
				arrayStructParticles[i].P.x = -OD + (i % 14)*OD;
				if (i % 14 == 0)
					++localZ;
				break;
			case 1:arrayParticles[i * 3 + j] = 5;
				arrayStructParticles[i].P.y = 5;
				break;
			case 2:arrayParticles[i * 3 + j] = -OD + (localZ)*OD;
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
		FuncionUpdate(dt / 7);
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