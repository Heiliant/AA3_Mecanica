#include <GL\glew.h>
#include <glm\gtc\type_ptr.hpp>
#include <glm\gtc\matrix_transform.hpp>


namespace GUIvars {
	extern bool PlaySimulation;
	extern int emitRate;

	//extern float SpherePos[3];
	//extern float SphereRad;
	extern float SphPosition[3];
	extern float SphRadius;

	bool renderSphere = true;
	bool renderCapsule = false;
	bool renderParticles = false;
	bool renderCloth = true;
	bool renderCube = false;
}

// Boolean variables allow to show/hide the primitives
//bool renderSphere = true;
//bool renderCapsule = false;
//bool renderParticles = false;
//bool renderCloth = true;
//bool renderCube = false;


namespace Sphere {
	extern void setupSphere(glm::vec3 pos = { GUIvars::SphPosition[0],  GUIvars::SphPosition[1], GUIvars::SphPosition[2] }, float radius = GUIvars::SphRadius);
	extern void cleanupSphere();
	extern void updateSphere(glm::vec3 pos = glm::vec3(&GUIvars::SphPosition[0], &GUIvars::SphPosition[1], &GUIvars::SphPosition[2]), float radius = GUIvars::SphRadius);
	extern void drawSphere();
}
namespace Capsule {
	extern void setupCapsule(glm::vec3 posA = glm::vec3(-3.f, 2.f, -2.f), glm::vec3 posB = glm::vec3(-4.f, 2.f, 2.f), float radius = 1.f);
	extern void cleanupCapsule();
	extern void updateCapsule(glm::vec3 posA, glm::vec3 posB, float radius = 1.f);
	extern void drawCapsule();
}
namespace LilSpheres {
	extern const int maxParticles;
	extern void setupParticles(int numTotalParticles, float radius = 0.05f);
	extern void cleanupParticles();
	extern void updateParticles(int startIdx, int count, float* array_data);
	extern void drawParticles(int startIdx, int count);
}
namespace ClothMesh {
	extern void setupClothMesh();
	extern void cleanupClothMesh();
	extern void updateClothMesh(float* array_data);
	extern void drawClothMesh();
}

namespace Cube {
	extern void setupCube();
	extern void cleanupCube();
	extern void updateCube(const glm::mat4& transform);
	extern void drawCube();
}

void setupPrims() {
	Sphere::setupSphere();
	Capsule::setupCapsule();
	LilSpheres::setupParticles(LilSpheres::maxParticles);
	ClothMesh::setupClothMesh();
	Cube::setupCube();
}
void cleanupPrims() {
	Sphere::cleanupSphere();
	Capsule::cleanupCapsule();
	LilSpheres::cleanupParticles();
	ClothMesh::cleanupClothMesh();
	Cube::cleanupCube();
}

void renderPrims() {
	if (GUIvars::renderSphere)
		Sphere::drawSphere();
	if (GUIvars::renderCapsule)
		Capsule::drawCapsule();

	if (GUIvars::renderParticles) {
		// You may need to rethink this piece of code...
		int startDrawingFromParticle = 0;
		int numParticlesToDraw = LilSpheres::maxParticles;
		LilSpheres::drawParticles(startDrawingFromParticle, numParticlesToDraw);
		// .............................................
	}
	
	if (GUIvars::renderCloth)
		ClothMesh::drawClothMesh();

	if (GUIvars::renderCube)
		Cube::drawCube();
}
