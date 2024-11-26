#include "CanvasTriangle.h"
#include <Utils.h>
#include <DrawingWindow.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <LoadObject.h>
#include <LoadMaterial.h>
#include <RayTriangleIntersection.h>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#define WIDTH 640
#define HEIGHT 480
// #define WIDTH 600
// #define HEIGHT 600
// #define WIDTH 1000
// #define HEIGHT 1000

std::vector<ModelTriangle> objectModelTriangles;
TextureMap texture;	
//task 4.5
glm::vec3 cameraPosition (0.0,0.0,4.0);
glm::vec3 centre(0.0,0.0,0.0);
float focalLength = 400;
float imageScalingPlane=1;
//task 5.3 camera orientation
//1st vec = right, 2nd vec =  up, 3rd vec = forward
glm::mat3 cameraOrientation (glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,1.0,0.0),glm::vec3(0.0,0.0,1.0));
//task 5.4 orbit toggle
bool orbitToggleLeft = false;
bool orbitToggleRight = false;
//task 6.5 light source position - little above the camera
glm::vec3 lightSource;
float lightRadius; //how many
//task 6.7 render mode
enum renderMode {wireframe, rasterised, rayTraced};
renderMode mode;
int renderModeInt;
//Global var
int max = 1000000;
float specularPower; //specular intensity
float specularShininess; //specular shininess
float sourceIntensity;
float ambient;
int shading; //0 flat 1 gouraud 2 phong
int shadow; //1 hard shadow 2 soft shadow



// Other functions  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//task 3.3 get points for triangle
CanvasTriangle randomStrokedTriangle(){
	
	CanvasPoint v0 (rand()%WIDTH,rand()%HEIGHT);
	CanvasPoint v1 (rand()%WIDTH,rand()%HEIGHT);
	CanvasPoint v2 (rand()%WIDTH,rand()%HEIGHT);
	CanvasTriangle randomStrokedTriangle(v0,v1,v2);
	return randomStrokedTriangle;
}


// Camera functions  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//task 5.2 camera translation matrix
//up down left right translation
void leftTranslation(glm::vec3 &cameraPosition){
	glm::vec3 translation(0.1,0.0,0.0);
	cameraPosition += translation;
}
void rightTranslation(glm::vec3 &cameraPosition){
	glm::vec3 translation(-0.1,0.0,0.0);
	cameraPosition += translation;
}
void upTranslation(glm::vec3 &cameraPosition){
	glm::vec3 translation(0.0,-0.1,0.0);
	cameraPosition += translation;
}
void downTranslation(glm::vec3 &cameraPosition){
	glm::vec3 translation(0.0,0.1,0.0);
	cameraPosition += translation;
}
//rotate camera on axis
void xRotateUp(glm::vec3 &cameraPosition,float theta){
	glm::mat3 rotationMatrix (glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,cos(theta/180),sin(theta/180)),glm::vec3(0.0,-sin(theta/180),cos(theta/180)));
	cameraPosition = rotationMatrix*cameraPosition;	
}
void xRotateDown(glm::vec3 &cameraPosition,float theta){
	glm::mat3 rotationMatrix (glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,cos(-theta/180),sin(-theta/180)),glm::vec3(0.0,-sin(-theta/180),cos(-theta/180)));
	cameraPosition = rotationMatrix*cameraPosition;
}
void yRotateLeft(glm::vec3 &cameraPosition,float theta){
	glm::mat3 rotationMatrix (glm::vec3(cos(theta/180),0.0,-sin(theta/180)),glm::vec3(0.0,1.0,0.0),glm::vec3(sin(theta/180),0.0,cos(theta/180)));
	cameraPosition = rotationMatrix*cameraPosition;
}
void yRotateRight(glm::vec3 &cameraPosition,float theta){
	glm::mat3 rotationMatrix (glm::vec3(cos(-theta/180),0.0,-sin(-theta/180)),glm::vec3(0.0,1.0,0.0),glm::vec3(sin(-theta/180),0.0,cos(-theta/180)));
	cameraPosition = rotationMatrix*cameraPosition;
}
void zRotateLeft(glm::vec3 &cameraPosition,float theta){
	glm::mat3 rotationMatrix (glm::vec3(cos(theta/180),sin(theta/180),0.0),glm::vec3(-sin(theta/180),cos(theta/180),0.0),glm::vec3(0.0,0.0,1.0));
	cameraPosition = rotationMatrix*cameraPosition;
}
void zRotateRight(glm::vec3 &cameraPosition,float theta){
	glm::mat3 rotationMatrix (glm::vec3(cos(-theta/180),sin(-theta/180),0.0),glm::vec3(-sin(-theta/180),cos(-theta/180),0.0),glm::vec3(0.0,0.0,1.0));
	cameraPosition = rotationMatrix*cameraPosition;
}
//task 5.3 camera orientation
void tiltup(glm::mat3 &cameraOrientation,float theta){
	glm::mat3 rotationMatrix (glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,cos(theta/180),sin(theta/180)),glm::vec3(0.0,-sin(theta/180),cos(theta/180)));
	cameraOrientation = rotationMatrix*cameraOrientation;
}
void tiltdown(glm::mat3 &cameraOrientation,float theta){
	glm::mat3 rotationMatrix (glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,cos(-theta/180),sin(-theta/180)),glm::vec3(0.0,-sin(-theta/180),cos(-theta/180)));
	cameraOrientation = rotationMatrix*cameraOrientation;
}
void panleft(glm::mat3 &cameraOrientation,float theta){
	glm::mat3 rotationMatrix (glm::vec3(cos(theta/180),0.0,-sin(theta/180)),glm::vec3(0.0,1.0,0.0),glm::vec3(sin(theta/180),0.0,cos(theta/180)));
	cameraOrientation =rotationMatrix*cameraOrientation;
}
void panright(glm::mat3 &cameraOrientation,float theta){
	glm::mat3 rotationMatrix (glm::vec3(cos(-theta/180),0.0,-sin(-theta/180)),glm::vec3(0.0,1.0,0.0),glm::vec3(sin(-theta/180),0.0,cos(-theta/180)));
	cameraOrientation = rotationMatrix*cameraOrientation;
}
//task 5.5 look at 
//alters 3 vectors, forward, up and right
//forward is zoomed in on the centre, up is aligned with the world up, right is perpendicular to forward and up
void lookAt(glm::vec3 &cameraPosition, glm::vec3 centre, glm::mat3 &cameraOrientation){
	//normalise vectors so they are directional vectors
	//find forward
	glm::vec3 forward = glm::normalize(cameraPosition-centre);
	//find right cross product vertical and forward
	glm::vec3 vertical(0.0,1.0,0.0);
	glm::vec3 right = glm::normalize(glm::cross(vertical,forward));
	//find up cross product forward and right normalised
	glm::vec3 up = glm::cross(forward,right);

	//change camera orientation
	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
	
}
//task 5.4 orbits
//orbits the position of the camera about the y axis, around centre	
void orbit(glm::vec3 &cameraPosition,float theta){

	//add toggle
	if (orbitToggleLeft==true){
		glm::mat3 rotationMatrix (glm::vec3(cos(theta/180),0.0,-sin(theta/180)),glm::vec3(0.0,1.0,0.0),glm::vec3(sin(theta/180),0.0,cos(theta/180)));
		//rotate from left to right
		cameraPosition = rotationMatrix*cameraPosition;

		//rotatefrom right to left
		//cameraPosition = cameraPosition*rotationMatrix;

		//look at to alter orientation with camera position
		glm::vec3 centre(0.0,0.0,0.0);
		lookAt(cameraPosition,centre,cameraOrientation);
	}
}
void orbitOther(glm::vec3 &cameraPosition,float theta){

	//add toggle
	if (orbitToggleRight==true){

		glm::mat3 rotationMatrix (glm::vec3(cos(theta/180),0.0,-sin(theta/180)),glm::vec3(0.0,1.0,0.0),glm::vec3(sin(theta/180),0.0,cos(theta/180)));
		//rotate from left to right
		//cameraPosition = rotationMatrix*cameraPosition;

		//rotatefrom right to left
		cameraPosition = cameraPosition*rotationMatrix;

		//look at to alter orientation with camera position
		glm::vec3 centre(0.0,0.0,0.0);
		lookAt(cameraPosition,centre,cameraOrientation);
	}
}



// Lighting + colour functions  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

//multipoint light source
//turn lihgt position into sphere
//find points on surface of sphere
std::vector <glm::vec3> lightSourceSphere(glm::vec3 lightPosition,float radius){

	//create sphere of light source
	int section = 10;	

	std::vector<std::vector<float>>  lightSourceSphere;
	for (float i = 0; i < section; i++)
	{
		float theta = i*2*M_PI/section;
		for (float j = 0; j < section; j++)
		{
			// code for float phi x y z were guided by https://www.geeksforgeeks.org/equation-of-a-sphere/
			float phi = j*M_PI/section;
			float x = lightPosition.x + radius*sin(phi)*cos(theta);
			float y = lightPosition.y + radius*sin(phi)*sin(theta);
			float z = lightPosition.z + radius*cos(phi);

			//ensure no duplicates
			std::vector<float> lightPoint = {x,y,z};
			//find if duplicate
			bool duplicate = false;
			for (std::vector<float> point : lightSourceSphere)
			{
				if (point==lightPoint){
					duplicate = true;
					break;
				}
			}

			if (duplicate==false){
				lightSourceSphere.push_back(lightPoint);
			}
			
		}
	}

	//print size of vecotr
	//std::cout << "Size of light source sphere: " << lightSourceSphere.size() << std::endl;
	std::vector<glm::vec3> lightSourceSphereGlm;
	for (std::vector<float> point : lightSourceSphere)
	{
		lightSourceSphereGlm.push_back(glm::vec3(point[0],point[1],point[2]));
	}

	return lightSourceSphereGlm;
}

//task get normal of triangle -UNUSED
//paramater vertices of the triangle
//return normal of the triangle as a vector
glm::vec3 getNormal(const glm::vec3 &v0,const  glm::vec3 &v1,const  glm::vec3 &v2){
	//find the direction from light to surface
	glm::vec3 e0 = v1 - v0;
	glm::vec3 e1 = v2 - v0;
	//cross product
	glm::vec3 normal = glm::cross(e0,e1);
	//normalise
	normal = glm::normalize(normal);
	return normal;
}

//task 6.6 calculate barycentric coordinates returns u v w for the triangle intersection
//consider making this more versitile
//paramater intersection point
//return barycentric coordinates vector
glm::vec3 barycentricCoordinates(const glm::vec3 &v0,const  glm::vec3 &v1,const  glm::vec3 &v2,const  glm::vec3 &p){
	//find the direction from light to surface
	glm::vec3 e0 = v1 - v0;
	glm::vec3 e1 = v2 - v0;
	glm::vec3 e2 = p - v0;

	//dot product
	float dot00 = glm::dot(e0,e0);
	float dot01 = glm::dot(e0,e1);
	float dot02 = glm::dot(e0,e2);
	float dot11 = glm::dot(e1,e1);
	float dot12 = glm::dot(e1,e2);

	//denominator
	float denominator = (dot00*dot11 - dot01*dot01);

	//barycentric coordinates
	float v = (dot11*dot02 - dot01*dot12)/denominator;
	float w = (dot00*dot12 - dot01*dot02)/denominator;
	float u = 1-w-v;

	return glm::vec3(u,v,w);
}

//bary cords for canvaspoint
glm::vec3 barycentricCoordinatesCanvasPoint(const CanvasPoint &v0, const CanvasPoint &v1, const CanvasPoint &v2 ,const  CanvasPoint &p){

	//find the direction from light to surface
	glm::vec3 e0 = glm::vec3(v1.x,v1.y,v1.depth) - glm::vec3(v0.x,v0.y,v0.depth);
	glm::vec3 e1 = glm::vec3(v2.x,v2.y,v2.depth) - glm::vec3(v0.x,v0.y,v0.depth);
	glm::vec3 e2 = glm::vec3(p.x,p.y,p.depth) - glm::vec3(v0.x,v0.y,v0.depth);

	//dot product
	float dot00 = glm::dot(e0,e0);
	float dot01 = glm::dot(e0,e1);
	float dot02 = glm::dot(e0,e2);
	float dot11 = glm::dot(e1,e1);
	float dot12 = glm::dot(e1,e2);

	//denominator
	float denominator = (dot00*dot11 - dot01*dot01);

	//barycentric coordinates
	float v = (dot11*dot02 - dot01*dot12)/denominator;
	float w = (dot00*dot12 - dot01*dot02)/denominator;
	float u = 1-w-v;

	return glm::vec3(u,v,w);
}
	
//get vertex normals of the intersected triangle
//paramater objectModelTriangles, intersection point
//return vertex normal as a vector
glm::vec3 getVertexNormal(const std::vector<ModelTriangle> &objectModelTriangles,const glm::vec3 &intersectionPointVertex, const RayTriangleIntersection &intersection){

	//find face normal of the intersected triangle
	//find faces that share the vertex
	std::vector<glm::vec3> faceNormals;
	float count = 0;
	for (ModelTriangle triangle : objectModelTriangles)
	{

		//if same triangele
		if (triangle.vertices==intersection.intersectedTriangle.vertices){
			count++;
			faceNormals.push_back(triangle.normal);
		}
		else{
					
			//same vex
			if (triangle.vertices[0]==intersectionPointVertex ){
				count++;
				faceNormals.push_back(triangle.normal);
			}
			else if (triangle.vertices[1]==intersectionPointVertex ){
				count++;
				faceNormals.push_back(triangle.normal);
			}
			else if (triangle.vertices[2]==intersectionPointVertex ){
				count++;
				faceNormals.push_back(triangle.normal);
			}
		}
	}
	//find the average of the face normals
	glm::vec3 vertexNormal = glm::vec3(0.0,0.0,0.0);
	for (glm::vec3 faceNormal : faceNormals)
	{
		vertexNormal+=faceNormal;
	}
	vertexNormal = vertexNormal/count;
	//normalise
	vertexNormal = glm::normalize(vertexNormal);
	return vertexNormal;
	
}

//Proximity: the closer a point on a surface is to a light, the brighter it will appear.
//brightness values should be within the range 0.0 (total darkness) to 1.0 (fully illuminated)
//no normalisation as need distance to light
//paramater light position, canvas position
//return brightness multiplier as a float
float proximityLighting (const glm::vec3 &lightPosition, const glm::vec3 &canvasPosition){

	//find distance from light to surface vec
	float lightToSurface = glm::length(lightPosition-canvasPosition);

	float brightness = 1/(4*M_PI*lightToSurface*lightToSurface);
	//brightness = std::max(0.0f,std::min(brightness,1.0f));	
	return brightness;

}

//Angle of incidence: the angle between the surface normal and the direction from the surface to the light source
//if the angle is 1, the light is parrallel to the surface and the surface will be fully illuminated. 
//If the angle is 0, the light is perpendicular to the surface and the surface will be in total darkness.
//paramater light position, canvas position, normal, source intensity
//return angle of incidence multiplier as a float
float angleOfIncidenceLighting (const glm::vec3 &lightPosition, const glm::vec3 &canvasPosition, const glm::vec3 &normal, float sourceIntensity){
	
	//find the direction from canvas to light
	glm::vec3 surfaceToLightVector = lightPosition - canvasPosition;

	//normalise since normal is a unit vector
	surfaceToLightVector = glm::normalize(surfaceToLightVector);
	//normal = glm::normalize(normal);

	//find the angle between the surface normal and surface to the light source
	float angleOfIncidence = glm::dot(normal,surfaceToLightVector);

	//make sure the angle not greater than 90 degrees
	angleOfIncidence = std::max(0.0f,angleOfIncidence);

	return angleOfIncidence;

}

//specular lighting: the angle between the direction from the surface to the camera and the direction from the surface to the light source
//if r reflection vector is the same as the direction from the surface to the camera, 
//dot product if dot product is 1 then the reflection vector is the same as the direction from the surface to the camera
//higher specular power = the shinier the surface but the smaller the area of the highlight
//paramater light position, canvas position, normal, camera position, specular power, source intensity
//return specular multiplier as a float
float specularLighting (const glm::vec3 &lightPosition, const glm::vec3 &canvasPosition, const glm::vec3 &normal, const glm::vec3 &cameraPosition, 
						float specularPower, float specularshininess){

	//vector from surface to camera
	glm::vec3 surfaceToCamera = cameraPosition - canvasPosition;
	//normalise
	surfaceToCamera = glm::normalize(surfaceToCamera);

	//rI vector of angle of incidence
	//vector of light to surface
	glm::vec3 rayIncidence = canvasPosition - lightPosition;
	//normalise
	rayIncidence = glm::normalize(rayIncidence);
	//find reflection vector
	//r reflection vector = r incident - 2(n dot r incident)
	glm::vec3 rayReflection = rayIncidence - 2.0f * normal * glm::dot(normal,rayIncidence);
	//normalise
	rayReflection = glm::normalize(rayReflection);

	//find the angle between the reflection vector and the surface to camera vector
	float specular = glm::dot(surfaceToCamera,rayReflection);
	//make sure the angle not greater than 90 degrees
	specular = std::max(0.0f,specular);

	//to the power of specular power
	specular = pow(specular,specularshininess);
	specular = specular*specularPower;


	//specular = std::max(0.0f,std::min(specular,1.0f));
	//full brightness =1 , no brightness = 0
	return specular;


}

//light is the light that is already present in the scene proximity + angle of incidence  specular = ambient light
//paramater proximity, angle of incidence, specular, source intensity, ambient factor
//return light multiplier as a float
float calculateLighting(float proximity, float angleOfIncidence, float specular, float sourceIntensity, float ambientfacotr){

	float light = sourceIntensity*(proximity*(angleOfIncidence+specular))+ambientfacotr;
	light = std::max(0.0f,std::min(light,1.0f));
	return light;
}

//task 6.7 Gouraud shading
//paramater light position, intersection, source intensity, specular power, object model triangles, ambient factor
//return light multiplier as a float
float gouraudShading(const glm::vec3 &lightPosition, const RayTriangleIntersection &intersection, const std::vector<ModelTriangle> &objectModelTriangles,float ambientfactor){
	
	//find barycentric coordinates
	glm::vec3 barycentric = barycentricCoordinates(intersection.intersectedTriangle.vertices[0],intersection.intersectedTriangle.vertices[1],intersection.intersectedTriangle.vertices[2],intersection.intersectionPoint);

	//find vertex normals of the intersected triangle
	glm::vec3 vertexNormal0 = getVertexNormal(objectModelTriangles,intersection.intersectedTriangle.vertices[0],intersection);
	glm::vec3 vertexNormal1 = getVertexNormal(objectModelTriangles,intersection.intersectedTriangle.vertices[1],intersection);
	glm::vec3 vertexNormal2 = getVertexNormal(objectModelTriangles,intersection.intersectedTriangle.vertices[2],intersection);

	//calculate lighting for each vertex
	float proximity0 = proximityLighting(lightPosition,intersection.intersectedTriangle.vertices[0]);
	float proximity1 = proximityLighting(lightPosition,intersection.intersectedTriangle.vertices[1]);
	float proximity2 = proximityLighting(lightPosition,intersection.intersectedTriangle.vertices[2]);

	float angleOfIncidence0 = angleOfIncidenceLighting(lightPosition,intersection.intersectedTriangle.vertices[0],vertexNormal0,sourceIntensity);
	float angleOfIncidence1 = angleOfIncidenceLighting(lightPosition,intersection.intersectedTriangle.vertices[1],vertexNormal1,sourceIntensity);
	float angleOfIncidence2 = angleOfIncidenceLighting(lightPosition,intersection.intersectedTriangle.vertices[2],vertexNormal2,sourceIntensity);

	float specular0 = specularLighting(lightPosition,intersection.intersectedTriangle.vertices[0],vertexNormal0,cameraPosition,specularPower,specularShininess);
	float specular1 = specularLighting(lightPosition,intersection.intersectedTriangle.vertices[1],vertexNormal1,cameraPosition,specularPower,specularShininess);
	float specular2 = specularLighting(lightPosition,intersection.intersectedTriangle.vertices[2],vertexNormal2,cameraPosition,specularPower,specularShininess);

	//calculate lighting for each vertex
	float light0 = calculateLighting(proximity0,angleOfIncidence0,specular0,sourceIntensity,ambientfactor);
	float light1 = calculateLighting(proximity1,angleOfIncidence1,specular1,sourceIntensity,ambientfactor);
	float light2 = calculateLighting(proximity2,angleOfIncidence2,specular2,sourceIntensity,ambientfactor);

	//find the lighting for the intersection point
	float calculateLighting = barycentric.x*light0 + barycentric.y*light1 + barycentric.z*light2;

	return calculateLighting;
}

//task 6.8 Phong shading maybe bump map
//paramater light position, intersection, source intensity, specular power, object model triangles, ambient factor
//return light multiplier as a float
float phongShading(const glm::vec3 &lightPosition, const RayTriangleIntersection &intersection, const std::vector<ModelTriangle> &objectModelTriangles,float ambientfactor){
	
	//find barycentric coordinates
	glm::vec3 barycentric = barycentricCoordinates(intersection.intersectedTriangle.vertices[0],intersection.intersectedTriangle.vertices[1],intersection.intersectedTriangle.vertices[2],intersection.intersectionPoint);
	//find vertex normals of the intersected triangle
	glm::vec3 vertexNormal0 = getVertexNormal(objectModelTriangles,intersection.intersectedTriangle.vertices[0],intersection);
	glm::vec3 vertexNormal1 = getVertexNormal(objectModelTriangles,intersection.intersectedTriangle.vertices[1],intersection);
	glm::vec3 vertexNormal2 = getVertexNormal(objectModelTriangles,intersection.intersectedTriangle.vertices[2],intersection);

	//interpolate vertex normal for the intersection point 
	glm::vec3 interpolatedNormal = barycentric.x*vertexNormal0 + barycentric.y*vertexNormal1 + barycentric.z*vertexNormal2;

	//normalise
	interpolatedNormal = glm::normalize(interpolatedNormal);

	//calculate lighting for the intersection point
	float proximity = proximityLighting(lightPosition,intersection.intersectionPoint);
	float angleOfIncidence = angleOfIncidenceLighting(lightPosition,intersection.intersectionPoint,interpolatedNormal,sourceIntensity);
	float specular = specularLighting(lightPosition,intersection.intersectionPoint,interpolatedNormal,cameraPosition,specularPower,specularShininess);

	float light = calculateLighting(proximity,angleOfIncidence,specular,sourceIntensity,ambientfactor);

	return light;


}

//overall lighting function
//paramater light position, camera position, closest valid intersection, object model triangles, ambient factor, source intensity, specular power, shading, is sphere
//return light multiplier as a vector
std::vector <float> overallLighting(const glm::vec3 &lightPosition, const glm::vec3 &cameraPosition, const RayTriangleIntersection &closestValidIntersection, const std::vector<ModelTriangle> &objectModelTriangles,
						float ambient){

	//lighting
	float calculatedLight;
	float brightness = 0.0;
	float brightnessShadows = 0.0;
	std::vector<float> light;

	if (closestValidIntersection.intersectedTriangle.object=="sphere"){
		specularShininess = 256;
	}
	else{
		specularShininess = 1;
	}	

	//calculate flat shading
	if (shading == 0){
		//normal
		float proximity = proximityLighting(lightPosition,closestValidIntersection.intersectionPoint);
		float angle = angleOfIncidenceLighting(lightPosition,closestValidIntersection.intersectionPoint,closestValidIntersection.intersectedTriangle.normal,sourceIntensity);
		float specular = specularLighting(lightPosition, closestValidIntersection.intersectionPoint,closestValidIntersection.intersectedTriangle.normal,cameraPosition,specularPower,specularShininess);
		calculatedLight = calculateLighting(proximity,angle,specular,sourceIntensity,ambient);
		brightness = calculatedLight;

		if (closestValidIntersection.intersectedTriangle.object=="sphere"){
			brightnessShadows= calculatedLight;
		}
		else{
			brightnessShadows = ambient;
		}


	}
	//gourand shading
	else if (shading == 1){
		calculatedLight = gouraudShading(lightPosition, closestValidIntersection, objectModelTriangles,ambient);
		brightness = calculatedLight;
		if (closestValidIntersection.intersectedTriangle.object=="sphere"){
			brightnessShadows= calculatedLight;
		}
		else{
			brightnessShadows = ambient;
		}


	}
	//phong shading
	else if (shading == 2){
		calculatedLight = phongShading(lightPosition, closestValidIntersection, objectModelTriangles,ambient);
		brightness = calculatedLight;
		if (closestValidIntersection.intersectedTriangle.object=="sphere"){
			brightnessShadows= calculatedLight;
		}
		else{
			brightnessShadows = ambient;
		}

	}

	light.push_back(brightness);
	light.push_back(brightnessShadows);
	
	return light;
}



// Interpolation ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//task 2.1 single interpolated return vecotr
std::vector<float> interpolateSingleFloats(float from,float to,int numberOfValues){

	//make vector
	std::vector<float> value;


	//if number of values is 1
	if (numberOfValues == 1){
		value.push_back(from);
		return value;
	}else{
		//find step
		float step = (to-from)/(numberOfValues-1);

		//interpolate function
		for (int i = 0; i < numberOfValues-1; i++)
		{
			value.push_back(from + (i*step));
		}
		return value;

	}


	
}

//task 2.2 interpolate 3 
std::vector<glm::vec3> interpolateThreeElementValues(const glm::vec3 &from, const glm::vec3 &to, float numberOfValues){
	//make vector
	std::vector<glm::vec3> value3;

	//if number of values is 1
	if (numberOfValues == 1){
		value3.push_back(from);
		return value3;
	}
	else{

		//vec3 for number of values for interpolation
		glm::vec3 numOfVals (numberOfValues-1);

		//step for each vecotr
		glm::vec3 step = (to-from)/(numOfVals);

		//interpolate and fill
		for (int i = 0; i < numberOfValues-1; i++)
		{

			//round 
			value3.push_back(glm::vec3(from+(float(i)*step)));
			
		}
		return value3;
	}

}

//task 3.5 interpolate with texture
std::vector<CanvasPoint> interpolatePointAndTexture (const CanvasPoint &from, const CanvasPoint &to, int numberOfValues){

	//make vector
	std::vector<CanvasPoint> value;

	if (numberOfValues == 1){
		CanvasPoint point = CanvasPoint(from.x,from.y);
		point.texturePoint=TexturePoint(from.texturePoint.x,from.texturePoint.y);
		value.push_back(point);
		return value;
	}
	else{


		//find step
		float xStep = (to.x-from.x)/(numberOfValues-1);
		float yStep = (to.y-from.y)/(numberOfValues-1);
		float xTextureStep = (to.texturePoint.x-from.texturePoint.x)/(numberOfValues-1);
		float yTextureStep = (to.texturePoint.y-from.texturePoint.y)/(numberOfValues-1);

		//interpolate function
		for (int i = 0; i < numberOfValues-1; i++)
		{
			CanvasPoint point = CanvasPoint(from.x+(i*xStep),from.y+(i*yStep));
			point.texturePoint=TexturePoint(from.texturePoint.x+(i*xTextureStep),from.texturePoint.y+(i*yTextureStep));
			value.push_back(point);
		}
		return value;
	}
	

}

//interpolate with depth
std::vector<CanvasPoint> interpolatePointAndDepth(const CanvasPoint &from, const CanvasPoint &to, int numberOfValues){

	//make vector
	std::vector<CanvasPoint> value;

	//print marker for debug
	// std::cout << "interpolating" << std::endl;
	// //print from and to
	// std::cout << "from: " << from.x << " " << from.y << " " << from.depth << std::endl;
	// std::cout << "to: " << to.x << " " << to.y << " " << to.depth << std::endl;
	
	if (numberOfValues == 1){
		CanvasPoint point = CanvasPoint(from.x,from.y,from.depth);
		value.push_back(point);
		return value;
	}else{
		//find step
		float xStep = (to.x-from.x)/(numberOfValues-1);
		float yStep = (to.y-from.y)/(numberOfValues-1);
		float depthStep = (to.depth-from.depth)/(numberOfValues-1);
		//interpolate function
		for (int i = 0; i < numberOfValues-1; i++)
		{
			CanvasPoint point = CanvasPoint(from.x+(i*xStep),from.y+(i*yStep),from.depth+(i*depthStep));
			//print point for debug
			//std::cout << "point: " << point.x << " " << point.y << " " << point.depth << std::endl;
			value.push_back(point);
		}
		return value;
	}


	

}



// Wireframe & Rasterisation  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//task 3.4 sort verticies v0=top v2=bottom v1=middle
//return sorted triangle
CanvasTriangle sortVertices(CanvasTriangle triangle){

	//v0=top v1=middle v2=botto

	if (triangle.v0().y>triangle.v1().y){

		std::swap(triangle.v0(),triangle.v1());
	}

	if (triangle.v0().y>triangle.v2().y){
		std::swap(triangle.v2(),triangle.v0());
	}

	if (triangle.v1().y>triangle.v2().y){
		std::swap(triangle.v1(),triangle.v2());
	}

	//if y are eaqual
	if (triangle.v0().y==triangle.v1().y){
		if (triangle.v0().x<triangle.v1().x){
			std::swap(triangle.v0(),triangle.v1());
		}
	}
	if (triangle.v0().y==triangle.v2().y){
		if (triangle.v0().x>triangle.v2().x){
			std::swap(triangle.v0(),triangle.v2());
		}
	}
	if (triangle.v1().y==triangle.v2().y){
		if (triangle.v1().x>triangle.v2().x){
			std::swap(triangle.v1(),triangle.v2());
		}
	}

	return triangle;
}

//task 3.4 find extra point
CanvasPoint fourthPoint(CanvasTriangle &triangle){


	//using ratios y2/y1==x2/x1 ... x2=x1*y2/y1
	float x1 = (triangle.v0().x-triangle.v2().x);
	float y1 = (triangle.v0().y-triangle.v2().y);
	float y2 = (triangle.v1().y-triangle.v0().y);
	float x2 = triangle.v0().x+ x1*y2/y1;

	CanvasPoint middle (x2,triangle.v1().y);

	return middle;
	
}

//task 3.5 fourth point for texture
//can use barycentric coordinates instead
CanvasPoint fourthPointTexture(CanvasTriangle &triangle){

	//using ratios y2/y1==x2/x1 ... x2=x1*y2/y1
	float x1 = (triangle.v2().x-triangle.v0().x);
	float y1 = (triangle.v2().y-triangle.v0().y);
	float y2 = (triangle.v1().y-triangle.v0().y);
	float x = triangle.v0().x+ x1*y2/y1;

	// std::cout << "ratio:  " << y2/y1 << std::endl;

	//use ratios for texture
	float fourthTextureX = triangle.v0().texturePoint.x + y2/y1 * (triangle.v2().texturePoint.x-triangle.v0().texturePoint.x);
	float fourthTextureY = triangle.v0().texturePoint.y + y2/y1 * (triangle.v2().texturePoint.y-triangle.v0().texturePoint.y);
	
	CanvasPoint middle (x,triangle.v1().y);
	middle.texturePoint=TexturePoint(fourthTextureX,fourthTextureY);

	// std::cout << "fourth:  " << middle << std::endl;
	// std::cout << "fourthtp:  " << middle.texturePoint << std::endl;



	return middle;
	
}

//task 4.9 find fourth point and depth
CanvasPoint fourthPointDepth(CanvasTriangle triangle){

	//using ratios y2/y1==x2/x1 ... x2=x1*y2/y1
	float x1 = (triangle.v2().x-triangle.v0().x);
	float y1 = (triangle.v2().y-triangle.v0().y);
	float y2 = (triangle.v1().y-triangle.v0().y);
	float x2 = triangle.v0().x+ x1*y2/y1;

	//find 1/depth
	triangle.v0().depth = 1/triangle.v0().depth;
	triangle.v1().depth = 1/triangle.v1().depth;
	triangle.v2().depth = 1/triangle.v2().depth;

	//find depth
	float z1 = triangle.v2().depth-triangle.v0().depth;
	float z = triangle.v0().depth+ z1*y2/y1;
	CanvasPoint middle (x2,triangle.v1().y,1/z);

	return middle;
	
}

//task 4.5 projection in 3D to 2D
CanvasPoint projectVertexOnToCanvasPoint(const glm::vec3 &cameraPosition, float focalLength,const glm::vec3 &vertexPosition){
	//using similar triangles
	//hv/dv=hi/di
	//find u v from x y z
	//f is focal length, w is width, h is height
	//u=fx/z+W/2
	//v=fy/z+H/2
	//model vertices are in camera coordinates, camera is at 0,0,0, x y is plane, z distance from camera
	//+z towads me, -z away from me
	//start with camera coord 0,0,4.0
	//focal length 2.0
	//std::cout << "camera: " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
	//std::cout<< "vertex: " << vertexPosition.x << " " << vertexPosition.y << " " << vertexPosition.z << std::endl;
	//transpose model vertices to camera coordinates
	//individual vertex
	//flip
	//print marker
	//std::cout << "vertex position:  " << vertexPosition.x << " " << vertexPosition.y << " " << vertexPosition.z << std::endl;
	
	glm::vec3 vertexPosition1;
	vertexPosition1 = vertexPosition - cameraPosition;

	glm::vec3 vert = vertexPosition1*cameraOrientation;

	float u = (-((focalLength*vert.x/vert.z)*imageScalingPlane)+WIDTH/2);
	float v = ((focalLength*vert.y/vert.z)*imageScalingPlane)+HEIGHT/2;

	//print u v
	//std::cout << "u: " << u << " v: " << v << std::endl;
	//add deoth buffer to z
	//interpolate with 1/z
	float depth = -1/vert.z;

	//round to int
	u = round(u);
	v = round(v);

	//std::cout << "u: " << u << " v: " << v << std::endl;
	//print x y z
	// std::cout << "x: " << vertexPosition.x << " y: " << vertexPosition.y << " z: " << vertexPosition.z << std::endl;
	//std::cout << "x: " << u << " y: " << v << " z: " << depth << std::endl;
	return CanvasPoint(u,v,depth);
}

//task 4.9 depth buffer 2D array
std::vector<std::vector<float>> depthBufferArray(DrawingWindow &window){
	//create 2d depth buffer
	std::vector<std::vector<float>> depthBuffer;
	for (size_t i = 0; i < window.width; i++)
	{
		std::vector<float> row;
		for (size_t j = 0; j < window.height; j++)
		{
			row.push_back(0);
		}
		depthBuffer.push_back(row);
	}
	return depthBuffer;
}



// RayTracing functions  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//task 6.2 ray tracing
//find intersection of ray with triangle
//ray start point s
//driection vector d
//positon on ray = startpoint + scalar * direction or r=s+td
//points on triangle p0 p1 p2 edge vectors e0 e1 e2
//r = p0 + u*e0 + v*e1 = s + td
//solve for t u v
// [t u v] = [d e0 e1]^-1 (s-p0)
//find ray direction from camera position to point on canvas
//returns closest intersection point
RayTriangleIntersection getClosestValidIntersection(const glm::vec3 &cameraPosition, const glm::vec3 &rayDirection, const std::vector<ModelTriangle> &objectModelTriangles, const TextureMap &texture){
	//for each triangle find if intersection with ray
	
	std::vector<RayTriangleIntersection> triangleIntersections;
	for (ModelTriangle triangle :objectModelTriangles)
	{
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		//possible solution 
		//t=absolute distance from camera to interesection
		//u= proportional distance along the triangle's first edge that the intersection
		//v= proportional distance along the triangle's second edge that the intersection
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		//find r intersection point on triangle relative to p0
		//approch 1 r = s + t*d
		//glm::vec3 r2 = cameraPosition + rayDirection * t;
		//approach 2 r = p0 + u*e0 + v*e1
		glm::vec3 r = triangle.vertices[0] + possibleSolution.y * e0 + possibleSolution.z * e1;
		// //print r1 and r2
		// std::cout << "r1:  " << r1.x << " " << r1.y << " " << r1.z << std::endl;
		// std::cout << "r2:  " << r2.x << " " << r2.y << " " << r2.z << std::endl;

		// //task 6.3 validate intersection
		//check if intersection is on triangle
		//make sure mirror isnt seeing self
		if (u>=0.0 && v>=0.0 && u+v<=1.0 && t>=0.0){

				//store triangle intersections
				RayTriangleIntersection intersection;
				intersection.intersectionPoint = r;
				intersection.distanceFromCamera = t;
				intersection.intersectedTriangle = triangle;
				//counter on how many intersections at point
				intersection.triangleIndex = triangleIntersections.size();

				//if colour 
				if (triangle.isTextured){

					//find textue point
					float xTexture =(1-u-v)*triangle.texturePoints[0].x + u*triangle.texturePoints[1].x + v*triangle.texturePoints[2].x;
					float yTexture =(1-u-v)*triangle.texturePoints[0].y + u*triangle.texturePoints[1].y + v*triangle.texturePoints[2].y;
					//round first
					xTexture = xTexture*texture.width;
					yTexture = yTexture*texture.height;
					//cover to colour
					uint32_t convCol = texture.pixels[round(xTexture)+(round(yTexture)*texture.width)];
					intersection.intersectedTriangle.colour = Colour((convCol & 0x00ff0000) >> 16, (convCol & 0x0000ff00) >> 8, (convCol & 0x000000ff));

				}
				triangleIntersections.push_back(intersection);


		}
			
	}
	RayTriangleIntersection closestIntersection;
	//may change to a large number
	closestIntersection.distanceFromCamera = max;

	for (RayTriangleIntersection intersection : triangleIntersections)
	{

		//find closest intersection that is not behind camera

		if (intersection.distanceFromCamera<closestIntersection.distanceFromCamera){
			closestIntersection = intersection;
		}

		//print marker
		//std::cout << "closest intersection:  " << closestIntersection.intersectionPoint.x << " " << closestIntersection.intersectionPoint.y << " " << closestIntersection.intersectionPoint.z << std::endl;
	}
	
	//return closest intersection
	//print marker
	//std::cout << "closest intersection:  " << closestIntersection.intersectionPoint.x << " " << closestIntersection.intersectionPoint.y << " " << closestIntersection.intersectionPoint.z << std::endl;


	return closestIntersection;
}

//task 6.4 convert 2d pixel to 3d ray direction
glm::vec3 rayDirectionVector(const glm::vec3 &cameraPosition, const glm::vec2 &pixel){

	float u = pixel.x;
	float v = pixel.y;
	float x = (-u+WIDTH/2)/imageScalingPlane;
	float y = (v-HEIGHT/2)/imageScalingPlane;
	float z = focalLength;

	glm::vec3 ray = cameraPosition- glm::vec3(x,y,z);
	//normalise
	glm::vec3 rayDirection = ray.x*cameraOrientation[0] + ray.y*cameraOrientation[1] + ray.z*cameraOrientation[2];
	rayDirection = glm::normalize(rayDirection);
	return rayDirection;

}

//task 6.5 draw ray tracing in the scene with shadows
glm::vec3 rayDirectionVectorShadows(const glm::vec3 &surfacePosition, const glm::vec3 &lightPosition){

	glm::vec3 rayDirection = surfacePosition - lightPosition;
	//normalise
	rayDirection = glm::normalize(rayDirection);

	return rayDirection;
}

//task 6.5 draw ray tracing in the scene with shadows
//can point see light??
//send ray from surface to light
bool shadowRay (RayTriangleIntersection intersection, glm::vec3 lightPosition, const std::vector<ModelTriangle> &objectModelTriangles){

	//find ray direction
	glm::vec3 rayDirection = lightPosition - intersection.intersectionPoint;
	//normalise
	rayDirection = glm::normalize(rayDirection);

	//distance from intersection to light
	float distance = glm::length(lightPosition - intersection.intersectionPoint);

	//print distance
	//std::cout << "distance:  " << distance << std::endl;

	bool isShadow=false;
	
	//doe the ray intersect with any other triangles
	std::vector<RayTriangleIntersection> triangleIntersections;
	for (ModelTriangle triangle :objectModelTriangles)
	{
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = intersection.intersectionPoint - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		//possible solution 
		//t=absolute distance from camera to interesection
		//u= proportional distance along the triangle's first edge that the intersection
		//v= proportional distance along the triangle's second edge that the intersection
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;

		//print t u v
		//std::cout << "t:  " << t << " u: " << u << " v: " << v << std::endl;

		//find r intersection point on triangle relative to p0
		//approch 1 r = s + t*d
		//glm::vec3 r2 = cameraPosition + rayDirection * t;
		//approach 2 r = p0 + u*e0 + v*e1
		glm::vec3 r = triangle.vertices[0] + possibleSolution.y * e0 + possibleSolution.z * e1;
		// //print r1 and r2
		// std::cout << "r1:  " << r1.x << " " << r1.y << " " << r1.z << std::endl;
		// std::cout << "r2:  " << r2.x << " " << r2.y << " " << r2.z << std::endl;

		// //task 6.3 validate intersection
		//check if intersection is on triangle
		if (u>=0.0 && v>=0.0 && u+v<=1.0 && t>0.01 && t<distance){

			//store triangle intersections
			RayTriangleIntersection intersection;
			intersection.intersectionPoint = r;
			intersection.distanceFromCamera = t;
			intersection.intersectedTriangle = triangle;
			//counter on how many intersections at point
			intersection.triangleIndex = triangleIntersections.size();
			//print 
			//std::cout << "shadow:  " << isShadow << std::endl;
			//triangleIntersections.push_back(intersection);
			isShadow = true;
			break;

		}
	}

	return isShadow;

}

//mirror
RayTriangleIntersection mirrorRay (RayTriangleIntersection intersection, glm::vec3 rayDirection, const std::vector<ModelTriangle> &objectModelTriangles){


	//print distance
	//std::cout << "distance:  " << distance << std::endl;

	std::vector<RayTriangleIntersection> triangleIntersections;
	for (ModelTriangle triangle :objectModelTriangles)
	{
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = intersection.intersectionPoint - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		//possible solution 
		//t=absolute distance from camera to interesection
		//u= proportional distance along the triangle's first edge that the intersection
		//v= proportional distance along the triangle's second edge that the intersection
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;

		//print t u v
		//std::cout << "t:  " << t << " u: " << u << " v: " << v << std::endl;

		//find r intersection point on triangle relative to p0
		//approch 1 r = s + t*d
		//glm::vec3 r2 = cameraPosition + rayDirection * t;
		//approach 2 r = p0 + u*e0 + v*e1
		glm::vec3 r = triangle.vertices[0] + possibleSolution.y * e0 + possibleSolution.z * e1;
		// //print r1 and r2
		// std::cout << "r1:  " << r1.x << " " << r1.y << " " << r1.z << std::endl;
		// std::cout << "r2:  " << r2.x << " " << r2.y << " " << r2.z << std::endl;

		// //task 6.3 validate intersection
		//check if intersection is on triangle
		if (u>=0.0 && v>=0.0 && u+v<=1.0 && t>0.0001 && t<100){

			//store triangle intersections
			RayTriangleIntersection intersection;
			intersection.intersectionPoint = r;
			intersection.distanceFromCamera = t;
			intersection.intersectedTriangle = triangle;
			//counter on how many intersections at point
			intersection.triangleIndex = triangleIntersections.size();
			//print 
			//std::cout << "shadow:  " << isShadow << std::endl;


			//if texture
			if (triangle.isTextured){

				//find textue point
				float xTexture =(1-u-v)*triangle.texturePoints[0].x + u*triangle.texturePoints[1].x + v*triangle.texturePoints[2].x;
				float yTexture =(1-u-v)*triangle.texturePoints[0].y + u*triangle.texturePoints[1].y + v*triangle.texturePoints[2].y;
				//round first
				xTexture = xTexture*texture.width;
				yTexture = yTexture*texture.height;
				//cover to colour
				uint32_t convCol = texture.pixels[round(xTexture)+(round(yTexture)*texture.width)];
				intersection.intersectedTriangle.colour = Colour((convCol & 0x00ff0000) >> 16, (convCol & 0x0000ff00) >> 8, (convCol & 0x000000ff));

			}
			
			triangleIntersections.push_back(intersection);

		}

	}	
	RayTriangleIntersection closestIntersection;
	closestIntersection.distanceFromCamera = max;
	
	for (RayTriangleIntersection intersection : triangleIntersections)
	{

		//find closest intersection that is not behind camera

		if (intersection.distanceFromCamera<closestIntersection.distanceFromCamera){
			closestIntersection = intersection;
		}

		//print marker
		//std::cout << "closest intersection:  " << closestIntersection.intersectionPoint.x << " " << closestIntersection.intersectionPoint.y << " " << closestIntersection.intersectionPoint.z << std::endl;
	}

	return closestIntersection;
}



// Drawing functions  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//task 3.2 draw line susing canvas points interpolate
void drawLine(DrawingWindow &window, CanvasPoint &from, CanvasPoint &to, Colour &colour){

	//get colour
	uint32_t convCol = (255<<24) + (colour.red<<16) + (colour.green<<8) +colour.blue;

	//find number of values
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	//find largest number of values
	int numberOfSteps = std::max(abs(xDiff),abs(yDiff));
	//step size
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;


	for (int i = 0; i < numberOfSteps; i++)
	{
		float x = from.x+(i*xStepSize);
		float y = from.y+(i*yStepSize);
		//round x and y
		x = round(x);
		y = round(y);

		if (x>=0 && x<window.width && y>=0 && y<window.height){
			window.setPixelColour(x,y,convCol);
		}

	}

}

//task 3.5 draw textured line
void drawTextureLine(DrawingWindow &window, CanvasPoint &from, CanvasPoint &to, TextureMap &texture){

	//similar to draw line
	//find number of values
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;

	//same for texture
	float xDiffTexture = to.texturePoint.x-from.texturePoint.x;
	float yDiffTexture = to.texturePoint.y-from.texturePoint.y;

	//find largest number of values
	int numberOfSteps = std::max(abs(xDiff),abs(yDiff));

	//step size
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	//step for texture
	float xStepSizeTexture = xDiffTexture/numberOfSteps;
	float yStepSizeTexture = yDiffTexture/numberOfSteps;

	//step texture to line
	for (int i = 0; i < numberOfSteps; i++)
	{
		float x = from.x+(i*xStepSize);
		float y = from.y+(i*yStepSize);

		if (x>=0 && x<window.width && y>=0 && y<window.height){

			float xTexture = from.texturePoint.x+(i*xStepSizeTexture);
			float yTexture = from.texturePoint.y+(i*yStepSizeTexture);


			//convert texture to colour
			uint32_t convCol = texture.pixels[round(xTexture)+(round(yTexture)*texture.width)];



			window.setPixelColour(x,y,convCol);
		}
	}

}

//task 4.9 draw line with depth buffer
void drawDepthLine(DrawingWindow &window, CanvasPoint &from, CanvasPoint &to, Colour &colour, std::vector<std::vector<float>> &depthBuffer){

	//get colour
	uint32_t convCol = (255<<24) + (colour.red<<16) + (colour.green<<8) +colour.blue;

	//print x y z
	// std::cout << "from:  " << from.x << " " << from.y << " " << from.depth << std::endl;
	// std::cout << "to:  " << to.x << " " << to.y << " " << to.depth << std::endl;
	//inverse depth for interpolation
	//from.depth = from.depth;
	//to.depth = to.depth;


	//find number of values
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
	//find largest number of values
	int numberOfSteps = std::max(abs(xDiff),abs(yDiff));
	//step size
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	float zStepSize = zDiff/numberOfSteps;


	for ( int i = 0; i < numberOfSteps; i++)
	{
		float x = from.x+(i*xStepSize);
		float y = from.y+(i*yStepSize);

		//round x and y
		x = round(x);
		y = round(y);

		if (x>=0 && x<window.width && y>=0 && y<window.height){

			float z = from.depth+(i*zStepSize);

			//check depth buffer
			if (z>depthBuffer[x][y]){
				//update depth buffer
				depthBuffer[x][y] = z;
;


				window.setPixelColour(x,y,convCol);
			}
		}

	}
}

//extra
void drawDepthTextureLine(DrawingWindow &window, CanvasPoint &from, CanvasPoint &to, TextureMap &texture, std::vector<std::vector<float>> &depthBuffer){

	//similar to draw line
	//find number of values
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;

	//same for texture
	float xDiffTexture = to.texturePoint.x-from.texturePoint.x;
	float yDiffTexture = to.texturePoint.y-from.texturePoint.y;

	//find largest number of values
	int numberOfSteps = std::max(abs(xDiff),abs(yDiff));

	//step size
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	//step for texture
	float xStepSizeTexture = xDiffTexture/numberOfSteps;
	float yStepSizeTexture = yDiffTexture/numberOfSteps;

	//step texture to line
	for (int i = 0; i < numberOfSteps; i++)
	{
		float x = from.x+(i*xStepSize);
		float y = from.y+(i*yStepSize);

		if (x>=0 && x<window.width && y>=0 && y<window.height){

			float xTexture = from.texturePoint.x+(i*xStepSizeTexture);
			float yTexture = from.texturePoint.y+(i*yStepSizeTexture);

			//convert texture to colour
			uint32_t convCol = texture.pixels[round(xTexture)+(round(yTexture)*texture.width)];
			//print marker
			//std::cout << "x: " << x << " y: " << y << " z: " << from.depth << std::endl;
			//std::cout << "xTexture: " << xTexture << " yTexture: " << yTexture << std::endl;

			//check depth buffer
			if (from.depth>depthBuffer[x][y]){
				//update depth buffer
				depthBuffer[x][y] = from.depth;

				window.setPixelColour(x,y,convCol);
			}
		}
	}

}

//task 3.3 draw stroked triangle
void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour){
	
	//draw triangle
	drawLine(window,triangle.v0(),triangle.v1(),colour);
	drawLine(window,triangle.v0(),triangle.v2(),colour);
	drawLine(window,triangle.v1(),triangle.v2(),colour);
}

//task 3.4 draw filled triangle
void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, Colour outlineColour){

	//sort vertices
	CanvasTriangle sortedTriangle = sortVertices(triangle);
	//find forth point
	CanvasPoint fourth = fourthPoint(sortedTriangle);
	//rename all points
	CanvasPoint top = sortedTriangle.v0();
	CanvasPoint left = fourth;
	CanvasPoint right = sortedTriangle.v1();
	CanvasPoint bottom = sortedTriangle.v2();
	if (left.x>right.x){
		std::swap(left,right);
	}

	//round to canvas point
	top = CanvasPoint(round(top.x),round(top.y));
	left = CanvasPoint(round(left.x),round(left.y));
	right = CanvasPoint(round(right.x),round(right.y));
	bottom = CanvasPoint(round(bottom.x),round(bottom.y));


	// std::cout << "top:  " << top << std::endl;
	// std::cout << "left:  " << left << std::endl;
	// std::cout << "right:  " << right << std::endl;
	// std::cout << "bottom:  " << bottom << std::endl; 

	//rows
	float topTriangleRows = std::abs(top.y -left.y)+1;
	float bottomTriangleRows = std::abs(left.y - bottom.y)+1;

	// std::cout << "top rows:  " << topTriangleRows << std::endl;
	// std::cout << "bottom rows:  " << bottomTriangleRows << std::endl;

	//interpolate top
	glm::vec3 fromT (top.x,top.x,top.y);
	//std::cout << "from:  " << from.x << " " << from.y << " " << from.z << std::endl;
	glm::vec3 toT (left.x,right.x,left.y);
	//std::cout << "to:  " << to.x << " " << to.y << " " << to.z << std::endl;
	std::vector<glm::vec3> topTriangle = interpolateThreeElementValues(fromT,toT,topTriangleRows);
	for(size_t i=0; i<topTriangle.size(); i++) 
	{
		glm::vec3& vec=topTriangle[i];
		//std::cout << "(" <<vec.x << ","<<vec.y<<","<<vec.z << " )";
		//std::cout << std::endl;

		//round to canvas point
		CanvasPoint from = CanvasPoint(round(vec.x),round(vec.z));
		CanvasPoint to = CanvasPoint(round(vec.y),round(vec.z));

		drawLine(window,from,to,colour);
	}

	//interpolate bottom
	//glm vector from
	glm::vec3 fromB (left.x,right.x,left.y);
	//std::cout << "from:  " << fromB.x << " " << fromB.y << " " << fromB.z << std::endl;
	glm::vec3 toB (bottom.x,bottom.x,bottom.y);
	//std::cout << "to:  " << toB.x << " " << toB.y << " " << toB.z << std::endl;
	std::vector<glm::vec3> bottomTriangle = interpolateThreeElementValues(fromB,toB,bottomTriangleRows);
	for(size_t i=0; i<bottomTriangle.size(); i++) 
	{
		glm::vec3& vec=bottomTriangle[i];
		//std::cout << "(" <<vec.x << ","<<vec.y<<","<<vec.z << " )";
		//std::cout << std::endl;

		//round to canvas point
		CanvasPoint from = CanvasPoint(round(vec.x),round(vec.z));
		CanvasPoint to = CanvasPoint(round(vec.y),round(vec.z));


		drawLine(window,from,to,colour);
	}


	
	//draw edge
	drawLine(window,top,left,outlineColour);
	drawLine(window,top,right,outlineColour);
	drawLine(window,left,bottom,outlineColour);
	drawLine(window,right,bottom,outlineColour);
	
	//print stoed triangle

}

//task 3.5 draw textured triangle from ppm
void drawFilledTextureTriangle (DrawingWindow &window, CanvasTriangle &triangle, TextureMap &texture){
	
	//std::cout << "h" << texture.height << std::endl;
	//std::cout << "w" << texture.width << std::endl;	

	//sort vertices
	CanvasTriangle sortedTriangle = sortVertices(triangle);
	//find forth point
	CanvasPoint fourth = fourthPointTexture(sortedTriangle);

	//check tp swapped
	//std::cout << "sorted:  " << sortedTriangle << std::endl;
	//std::cout << "fourth:  " << fourth << std::endl;
	//std::cout << "ratio:  " << ratio << std::endl;
	//std::cout << "fourthtp:  " << fourth.texturePoint << std::endl;

	//sort vertices
	CanvasPoint top = sortedTriangle.v0();
	CanvasPoint left = fourth;
	CanvasPoint right = sortedTriangle.v1();
	CanvasPoint bottom = sortedTriangle.v2();
	if (left.x>right.x){
		std::swap(left,right);
	}

	// std::cout << "top:  " << top << std::endl;
	// std::cout << "TOP TEXTURE:  " << top.texturePoint << std::endl;
	// std::cout << "left:  " << left << std::endl;
	// std::cout << "LEFT TEXTURE:  " << left.texturePoint << std::endl;
	// std::cout << "right:  " << right << std::endl;
	// std::cout << "RIGHT TEXTURE:  " << right.texturePoint << std::endl;
	// std::cout << "bottom:  " << bottom << std::endl;
	// std::cout << "BOTTOM TEXTURE:  " << bottom.texturePoint << std::endl;


	//find number of rows
	float topTriangleRows = abs(top.y-left.y)+1;
	float bottomTriangleRows = abs(left.y-bottom.y)+1;

	//interpolate top with texture
	std::vector<CanvasPoint> topToLeftTriangle = interpolatePointAndTexture(top,left,topTriangleRows);
	std::vector<CanvasPoint> topToRightTriangle = interpolatePointAndTexture(top,right,topTriangleRows);
	for (size_t i = 0; i < topTriangleRows ; i++)
	{
		//draw line
		CanvasPoint from = CanvasPoint(round(topToLeftTriangle[i].x),round(topToLeftTriangle[i].y));
		from.texturePoint=TexturePoint(topToLeftTriangle[i].texturePoint.x,topToLeftTriangle[i].texturePoint.y);
		CanvasPoint to = CanvasPoint(round(topToRightTriangle[i].x),round(topToRightTriangle[i].y));
		to.texturePoint=TexturePoint(topToRightTriangle[i].texturePoint.x,topToRightTriangle[i].texturePoint.y);

		//print texture point
		// std::cout << "from:  " << from.texturePoint.x << " " << from.texturePoint.y << std::endl;
		// std::cout << "to:  " << to.texturePoint.x << " " << to.texturePoint.y << std::endl;


		drawTextureLine(window,from,to,texture);

	}

	//interpolate bottom with texture
	std::vector<CanvasPoint> leftToBottomTriangle = interpolatePointAndTexture(left,bottom,bottomTriangleRows);
	std::vector<CanvasPoint> rightToBottomTriangle = interpolatePointAndTexture(right,bottom,bottomTriangleRows);
	for (size_t i = 0; i < bottomTriangleRows ; i++)
	{
		//draw line
		CanvasPoint from = CanvasPoint(round(leftToBottomTriangle[i].x),round(leftToBottomTriangle[i].y));
		from.texturePoint=TexturePoint(leftToBottomTriangle[i].texturePoint.x,leftToBottomTriangle[i].texturePoint.y);
		CanvasPoint to = CanvasPoint(round(rightToBottomTriangle[i].x),round(rightToBottomTriangle[i].y));
		to.texturePoint=TexturePoint(rightToBottomTriangle[i].texturePoint.x,rightToBottomTriangle[i].texturePoint.y);
		drawTextureLine(window,from,to,texture);

	}

	//draw edge
	// Colour white = Colour(255,255,255);
	// drawStrokedTriangle(window, sortedTriangle,white);

	//draw texture line edge
	drawTextureLine(window,top,left,texture);
	drawTextureLine(window,top,right,texture);
	drawTextureLine(window,left,bottom,texture);
	drawTextureLine(window,right,bottom,texture);


}

//task 4.9 draw filled triangle with depth buffer
void drawFilledTriangleDepthBuffer(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, Colour outlineColour){

	//sort vertices
	CanvasTriangle sortedTriangle = sortVertices(triangle);

	//find forth point
	CanvasPoint fourth = fourthPointDepth(sortedTriangle);
	//rename all points
	CanvasPoint top = sortedTriangle.v0();
	CanvasPoint left = fourth;
	CanvasPoint right = sortedTriangle.v1();
	CanvasPoint bottom = sortedTriangle.v2();	

	if (left.x>right.x){
		std::swap(left,right);
	}
	//round to canvas point
	top = CanvasPoint(round(top.x),round(top.y),top.depth);
	left = CanvasPoint(round(left.x),round(left.y), left.depth);
	right = CanvasPoint(round(right.x),round(right.y), right.depth);
	bottom = CanvasPoint(round(bottom.x),round(bottom.y),bottom.depth);
	//rows
	float topTriangleRows = std::abs(top.y -left.y)+1;
	float bottomTriangleRows = std::abs(left.y - bottom.y)+1;

	//interpolate top
	std::vector<CanvasPoint> topToLeft=interpolatePointAndDepth(top,left,topTriangleRows);
	std::vector<CanvasPoint> topToRight=interpolatePointAndDepth(top,right,topTriangleRows);


	for(size_t i=0; i<topToLeft.size(); i++)
	{
		//round x and y for canvas point
		CanvasPoint from = CanvasPoint(round(topToLeft[i].x),round(topToLeft[i].y),topToLeft[i].depth);
		CanvasPoint to = CanvasPoint(round(topToRight[i].x),round(topToRight[i].y),topToRight[i].depth);

		drawDepthLine(window,from,to,colour,depthBuffer);
	
	}

	

	//interpolate bottom
	std::vector<CanvasPoint> leftToBottom=interpolatePointAndDepth(left,bottom,bottomTriangleRows);
	std::vector<CanvasPoint> rightToBottom=interpolatePointAndDepth(right,bottom,bottomTriangleRows);


	for(size_t i=0; i<leftToBottom.size(); i++)
	{
		//round to canvas point
		CanvasPoint from = CanvasPoint(round(leftToBottom[i].x),round(leftToBottom[i].y),leftToBottom[i].depth);
		CanvasPoint to = CanvasPoint(round(rightToBottom[i].x),round(rightToBottom[i].y),rightToBottom[i].depth);

		drawDepthLine(window,from,to,colour,depthBuffer);

	}


	//draw edge
	drawDepthLine(window,top,left,outlineColour,depthBuffer);
	drawDepthLine(window,top,right,outlineColour,depthBuffer);
	drawDepthLine(window,left,bottom,outlineColour,depthBuffer);
	drawDepthLine(window,right,bottom,outlineColour,depthBuffer);




}

//task 4.6 pointcloud render
void drawPointCloud(DrawingWindow &window, const glm::vec3 &cameraPosition, float focalLength,const std::vector<ModelTriangle> &objectModelTriangles){
	//for each triangle
	for (ModelTriangle triangle : objectModelTriangles)
	{
		//transform to canvas points
		CanvasPoint v0 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[0]);
		CanvasPoint v1 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[1]);
		CanvasPoint v2 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[2]);

		//draw on canvas
		uint32_t colour = (255<<24) + (triangle.colour.red<<16) + (triangle.colour.green<<8) +triangle.colour.blue;
		window.setPixelColour(v0.x,v0.y,colour);
		window.setPixelColour(v1.x,v1.y,colour);
		window.setPixelColour(v2.x,v2.y,colour);
	}
		
}

//task 4.7 wireframe render
void drawWireframe(DrawingWindow &window, const glm::vec3 &cameraPosition, float focalLength,const std::vector<ModelTriangle> &objectModelTriangles){

	//for each triangle
	for (ModelTriangle triangle : objectModelTriangles)
	{
		//transform to canvas points
		CanvasPoint v0 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[0]);
		CanvasPoint v1 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[1]);
		CanvasPoint v2 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[2]);

		// if (triangle.isTextured){
		// 	//if there are texture points on the triangle add them
		// 	v0.texturePoint = triangle.texturePoints[0];
		// 	v1.texturePoint = triangle.texturePoints[1];
		// 	v2.texturePoint = triangle.texturePoints[2];
		// 	CanvasTriangle sortedTriangle = sortVertices(CanvasTriangle(v0,v1,v2));
		// 	drawTextureLine(window,sortedTriangle.v0(),sortedTriangle.v1(),texture);
		// } else{
		// 	//draw on canvas
		// 	drawStrokedTriangle(window, CanvasTriangle(v0,v1,v2),triangle.colour);
		// }
		//white outline
		Colour white = Colour(255,255,255);
		drawStrokedTriangle(window, CanvasTriangle(v0,v1,v2),white);
	}
}

//use barycentric coordinates to find if point is inside triangle
void drawTextureMapping(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2 , TextureMap &texture, std::vector<std::vector<float>> &depthBuffer){


	// //sort canvas points
	CanvasTriangle sortedTriangle = sortVertices(CanvasTriangle(v0,v1,v2));

	//find forth point
	CanvasPoint fourth = fourthPoint(sortedTriangle);
	//find barrycentric coordinates
	glm::vec3 forthbarycentric = barycentricCoordinatesCanvasPoint(sortedTriangle.v0(),sortedTriangle.v1(),sortedTriangle.v2(),fourth);
	//find depth
	fourth.depth= forthbarycentric.x*sortedTriangle.v0().depth + forthbarycentric.y*sortedTriangle.v1().depth + forthbarycentric.z*sortedTriangle.v2().depth;
	//find texture point
	fourth.texturePoint.x = forthbarycentric.x*sortedTriangle.v0().texturePoint.x + forthbarycentric.y*sortedTriangle.v1().texturePoint.x + forthbarycentric.z*sortedTriangle.v2().texturePoint.x;
	fourth.texturePoint.y = forthbarycentric.x*sortedTriangle.v0().texturePoint.y + forthbarycentric.y*sortedTriangle.v1().texturePoint.y + forthbarycentric.z*sortedTriangle.v2().texturePoint.y;

	//find top lef t right bottom
	CanvasPoint top = sortedTriangle.v0();
	CanvasPoint left = fourth;
	CanvasPoint right = sortedTriangle.v1();
	CanvasPoint bottom = sortedTriangle.v2();
	if (left.x>right.x){
		std::swap(left,right);
	}
	top.x = round(top.x);
	top.y = round(top.y);
	left.x = round(left.x);
	left.y = round(left.y);
	right.x = round(right.x);
	right.y = round(right.y);
	bottom.x = round(bottom.x);
	bottom.y = round(bottom.y);


	//interpolate top
	float topTriangleRows = std::abs(top.y -left.y)+1;
	std::vector<CanvasPoint> topToLeft=interpolatePointAndDepth(top,left,topTriangleRows);
	std::vector<CanvasPoint> topToRight=interpolatePointAndDepth(top,right,topTriangleRows);
	for(int i = 0; i<topToLeft.size(); i++)
	{
		//find length for each rounded i row
		topToLeft[i].x = round(topToLeft[i].x);
		topToLeft[i].y = round(topToLeft[i].y);
		int length = std::abs(std::round(topToLeft[i].x - topToRight[i].x));
	

		for (int j=0; j< length; j++){

			CanvasPoint point = CanvasPoint(round(topToLeft[i].x)+j,round(topToLeft[i].y),topToLeft[i].depth);
			//check point on window
			if (point.x>=0 && point.x<window.width && point.y>=0 && point.y<window.height){
				//use barycentric coordinates
				glm::vec3 barycentric = barycentricCoordinatesCanvasPoint(top,left,right,point);
				//use for finding texture point
				point.texturePoint.x = barycentric.x*top.texturePoint.x + barycentric.y*left.texturePoint.x + barycentric.z*right.texturePoint.x;
				point.texturePoint.y = barycentric.x*top.texturePoint.y + barycentric.y*left.texturePoint.y + barycentric.z*right.texturePoint.y;	

				//check depth buffer
				if (point.depth>=depthBuffer[point.x][point.y]){
					//update depth buffer
					depthBuffer[point.x][point.y] = point.depth;
					//convert texture to colour
					uint32_t convCol = texture.pixels[round(point.texturePoint.x)+(round(point.texturePoint.y)*texture.width)];
					//set pixel colour
					window.setPixelColour(point.x,point.y,convCol);
				}
			}
		}
	}

	// //interpolate bottom
	float bottomTriangleRows = std::abs(left.y - bottom.y)+1;
	std::vector<CanvasPoint> leftToBottom=interpolatePointAndDepth(left,bottom,bottomTriangleRows);
	std::vector<CanvasPoint> rightToBottom=interpolatePointAndDepth(right,bottom,bottomTriangleRows);

	for (int i=0 ; i<leftToBottom.size(); i++)
	{
		//find length for each rounded i row
		leftToBottom[i].x = round(leftToBottom[i].x);
		leftToBottom[i].y = round(leftToBottom[i].y);
		int length = std::abs(std::round(leftToBottom[i].x - rightToBottom[i].x));
		for (int j=0; j< length; j++){

			CanvasPoint point = CanvasPoint(round(leftToBottom[i].x)+j,round(leftToBottom[i].y),leftToBottom[i].depth);
			//check point on window
			if (point.x>=0 && point.x<window.width && point.y>=0 && point.y<window.height){
				//use barycentric coordinates
				glm::vec3 barycentric = barycentricCoordinatesCanvasPoint(left,right,bottom,point);
				//use for finding texture point
				point.texturePoint.x = barycentric.x*left.texturePoint.x + barycentric.y*right.texturePoint.x + barycentric.z*bottom.texturePoint.x;
				point.texturePoint.y = barycentric.x*left.texturePoint.y + barycentric.y*right.texturePoint.y + barycentric.z*bottom.texturePoint.y;	

				//check depth buffer
				if (point.depth>=depthBuffer[point.x][point.y]){
					//update depth buffer
					depthBuffer[point.x][point.y] = point.depth;
					//convert texture to colour
					uint32_t convCol = texture.pixels[round(point.texturePoint.x)+(round(point.texturePoint.y)*texture.width)];
					//set pixel colour
					window.setPixelColour(point.x,point.y,convCol);
				}
			}
		}
	}


}

//task 4.8 rasterize render
void drawRasterisedScene(DrawingWindow &window, const glm::vec3 &cameraPosition, float focalLength,const std::vector<ModelTriangle> &objectModelTriangles){

	//for each triangle
	for (ModelTriangle triangle : objectModelTriangles)
	{
		//transform to canvas points
		CanvasPoint v0 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[0]);
		CanvasPoint v1 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[1]);
		CanvasPoint v2 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[2]);

		if (triangle.isTextured){
			//if there are texture points on the triangle add them
			//nothing
		} else{

			//print marker
		//std::cout << "here" << std::endl;
		//draw on canvas
		drawFilledTriangle(window,CanvasTriangle(v0,v1,v2),triangle.colour,triangle.colour);
		}
	}
}

//task 4.9 wireframe render with depth buffer
void drawRasterisedSceneDepth(DrawingWindow &window, const glm::vec3 &cameraPosition, float focalLength, const std::vector<ModelTriangle> &objectModelTriangles){

	std::vector<std::vector<float>> depthBuffer = depthBufferArray(window);


	for (ModelTriangle triangle : objectModelTriangles)
	{	

		CanvasPoint v0 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[0]);
		CanvasPoint v1 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[1]);
		CanvasPoint v2 = projectVertexOnToCanvasPoint(cameraPosition,focalLength,triangle.vertices[2]);

		if (triangle.isTextured){
			
			v0.texturePoint.x = triangle.texturePoints[0].x*texture.width;
			v0.texturePoint.y = triangle.texturePoints[0].y*texture.height;
			v1.texturePoint.x = triangle.texturePoints[1].x*texture.width;
			v1.texturePoint.y = triangle.texturePoints[1].y*texture.height;
			v2.texturePoint.x = triangle.texturePoints[2].x*texture.width;
			v2.texturePoint.y = triangle.texturePoints[2].y*texture.height;


			drawTextureMapping(window,v0,v1,v2,texture,depthBuffer);
		
			
		} else{

			drawFilledTriangleDepthBuffer(window,CanvasTriangle(v0,v1,v2),triangle.colour,depthBuffer,triangle.colour);
		}
	}

}

//task 6.4 draw ray tracing in the scene
void drawRayTrace(DrawingWindow &window, const glm::vec3 &cameraPosition, const std::vector<ModelTriangle> &ModelTriangles,const glm::vec3 &lightPosition,
								float ambient, float sourceIntensity, float specularPower, int shading){

	//for each pixel on canvas
	for (size_t y = 0; y < window.height; y++)
	{
		for (size_t x = 0; x < window.width; x++)
		{
			//find ray direction convert 2d to 3d
			//2d pixel
			glm::vec2 pixel(x,y);	
			glm::vec3 rayDirection = rayDirectionVector(cameraPosition,pixel);
			
				//find closest intersection
			RayTriangleIntersection closestValidIntersection = getClosestValidIntersection(cameraPosition,rayDirection,ModelTriangles,texture);



			//draw on canvas
			if (closestValidIntersection.distanceFromCamera<max){

				Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red,closestValidIntersection.intersectedTriangle.colour.green,closestValidIntersection.intersectedTriangle.colour.blue);
			
				uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
				window.setPixelColour(x,y,colour);


			}
			else{
				//draw black
				uint32_t colour = (255<<24) + (0<<16) + (0<<8) +0;
				window.setPixelColour(x,y,colour);
			}
		}
	}

}

//task 6.5 draw ray tracing in the scene with shadows
//if light to surface ray intersects with same  triangle
void drawRayTraceWithShadow(DrawingWindow &window, const glm::vec3 &cameraPosition, const std::vector<ModelTriangle> &ModelTriangles, const glm::vec3 &lightPosition,
								float ambient, float sourceIntensity, float specularPower, int shading){

	//flat shading
	if(shadow==1){
		for (size_t y = 0; y < window.height; y++)
		{
			for (size_t x = 0; x < window.width; x++)
			{
				//find ray direction convert 2d to 3d
				//2d pixel
				glm::vec2 pixel(x,y);	
				glm::vec3 rayDirection = rayDirectionVector(cameraPosition,pixel);
				//find closest intersection
				RayTriangleIntersection closestValidIntersection = getClosestValidIntersection(cameraPosition,rayDirection,ModelTriangles,texture);
				//find if intersection has shadow
				bool isShadow = shadowRay(closestValidIntersection,lightPosition,ModelTriangles);
				std::vector<float> lightIntensity = overallLighting(lightPosition,cameraPosition,closestValidIntersection,ModelTriangles,ambient);
				float brightness = lightIntensity[0];
				float brightnessShadows = lightIntensity[1];

				//draw on canvas
				if (closestValidIntersection.distanceFromCamera<max && isShadow==false){
					Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red*brightness,closestValidIntersection.intersectedTriangle.colour.green*brightness,closestValidIntersection.intersectedTriangle.colour.blue*brightness);
					uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
					window.setPixelColour(x,y,colour);
				}
				else{
					//draw black
					Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red*brightnessShadows,closestValidIntersection.intersectedTriangle.colour.green*brightnessShadows,closestValidIntersection.intersectedTriangle.colour.blue*brightnessShadows);
					uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
					window.setPixelColour(x,y,colour);
				}
			}
		}
	}

	//for gourand and phong
	else{
	std::vector <glm::vec3> lightSphere = lightSourceSphere(lightPosition,lightRadius);
	float lightSSize = lightSphere.size();
	

		for (int y = 0; y < window.height; y++)
		{
			for(int x = 0; x < window.width; x++)
			{
				
				//find ray direction convert 2d to 3d
				//2d pixel
				glm::vec2 pixel(x,y);	
				glm::vec3 rayDirection = rayDirectionVector(cameraPosition,pixel);


				//find closest intersection
				RayTriangleIntersection closestValidIntersection = getClosestValidIntersection(cameraPosition,rayDirection,ModelTriangles,texture);

		
				RayTriangleIntersection reflectionIntersection;

				//if itersection is mirror, find reflection
				if (closestValidIntersection.intersectedTriangle.colour.name=="Mirror"){
					//find reflection
					glm::vec3 rayIncidence = normalize(rayDirection);
					//print out normals
					//std::cout << "normal:  " << closestValidIntersection.intersectedTriangle.normal.x << " " << closestValidIntersection.intersectedTriangle.normal.y << " " << closestValidIntersection.intersectedTriangle.normal.z << std::endl;
					//find reflection vector
					glm::vec3 normal = normalize(closestValidIntersection.intersectedTriangle.normal);
					//print out normal
					//std::cout << "normal:  " << normal.x << " " << normal.y << " " << normal.z << std::endl;
					glm::vec3 rayReflection = rayIncidence - (2.0f * normal * glm::dot(rayIncidence,normal));
					//normalise
					rayReflection = glm::normalize(rayReflection);
					//find intersection point
					reflectionIntersection = mirrorRay(closestValidIntersection,rayReflection,ModelTriangles);
					//store colour

					closestValidIntersection.intersectedTriangle.colour = reflectionIntersection.intersectedTriangle.colour;


					//if colour is black print mirro object name
					// if (closestValidIntersection.intersectedTriangle.colour.red==0 && closestValidIntersection.intersectedTriangle.colour.green==0 && closestValidIntersection.intersectedTriangle.colour.blue==0){
					// 	//print marker
					// 	std::cout << "mirror object:  " << reflectionIntersection.intersectedTriangle.object << std::endl;
					// 	//print distance
					// 	std::cout << "distance:  " << reflectionIntersection.distanceFromCamera << std::endl;
					// }

					//find if reflection intersects has shadow
					bool isShadow = shadowRay(reflectionIntersection,lightPosition,ModelTriangles);
					bool hardShadow = true;
					int countshadow = 0;
					if (isShadow){
						for (glm::vec3 light : lightSphere){
							//check if shadow
							if (shadowRay(reflectionIntersection,light,ModelTriangles)==false && light!=lightPosition){
								hardShadow = false;
								countshadow++;
							}
						}
					}

					//find light intensity
					std::vector<float> lightIntensity = overallLighting(lightPosition,cameraPosition,reflectionIntersection,ModelTriangles,ambient);
					float brightness = lightIntensity[0];
					float brightnessShadows = lightIntensity[1];

					if (isShadow==false){
						Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red*brightness,closestValidIntersection.intersectedTriangle.colour.green*brightness,closestValidIntersection.intersectedTriangle.colour.blue*brightness);
						uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
						window.setPixelColour(x,y,colour);
					}
					else{
						float lightSoft = (countshadow/lightSSize);

						if (lightSoft<=ambient){
							//brightnessShadows = ambient;
						}
						else if (lightSoft>=brightness){
							brightnessShadows = brightness;
						}
						else{
							brightnessShadows = lightSoft;
						}

						Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red*brightnessShadows,closestValidIntersection.intersectedTriangle.colour.green*brightnessShadows,closestValidIntersection.intersectedTriangle.colour.blue*brightnessShadows);
						uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
						window.setPixelColour(x,y,colour);
					}


				
				}

				else if (closestValidIntersection.distanceFromCamera<max){

					//print object colour naem
					//std::cout << "object:  " << closestValidIntersection.intersectedTriangle.colour.name << std::endl;
					// if (closestValidIntersection.intersectedTriangle.colour.name=="White"){
					// 	//print marker
					// 	std::cout << "white object" << std::endl;
					// 	//print vertices
					// 	std::cout << "v0:  " << closestValidIntersection.intersectedTriangle.vertices[0].x << " " << closestValidIntersection.intersectedTriangle.vertices[0].y << " " << closestValidIntersection.intersectedTriangle.vertices[0].z << std::endl;
					// 	//print pixel
					// 	std::cout << "pixel:  " << x << " " << y << std::endl;
					// }

					RayTriangleIntersection shadowIntersection;
					bool isShadow;
					bool hardShadow = true;
					int countshadow = 0;
					
					if (closestValidIntersection.intersectedTriangle.colour.name=="Mirror"){

						isShadow=false;

					}
					else{
						shadowIntersection = closestValidIntersection;
						isShadow =  shadowRay(shadowIntersection,lightPosition,ModelTriangles);

					}
					//change for soft shadows

					if (isShadow){

						for (glm::vec3 light : lightSphere){
							//check if shadow
							if (shadowRay(shadowIntersection,light,ModelTriangles)==false && light!=lightPosition){
								hardShadow = false;
								//print marker
								//std::cout << "soft shadow" << std::endl;
								//shadowMapSoft[closestValidIntersection.intersectedTriangle.object].push_back(closestValidIntersection.intersectionPoint);
								//break;
								countshadow++;
							}
						}

						// //print hardsadow bool
						// //std::cout << "hard shadow:  " << hardShadow << std::endl;

						// if (hardShadow==false){
						// 	shadowMapSoft[closestValidIntersection.intersectedTriangle.object].push_back(closestValidIntersection.intersectionPoint);
						// }
						// else{
						// 	shadowMapHard[closestValidIntersection.intersectedTriangle.object].push_back(closestValidIntersection.intersectionPoint);
						// }
					}


					//find light intensity
					std::vector<float> lightIntensity = overallLighting(lightPosition,cameraPosition,shadowIntersection,ModelTriangles,ambient);
					float brightness = lightIntensity[0];
					float brightnessShadows = lightIntensity[1];

					
					

					if (isShadow==false){


						Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red*brightness,closestValidIntersection.intersectedTriangle.colour.green*brightness,closestValidIntersection.intersectedTriangle.colour.blue*brightness);
						uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
						window.setPixelColour(x,y,colour);
					}
					else{
						//std::cout << "count:  " << countshadow << std::endl;
						float lightSoft = (countshadow/lightSSize);

						if (lightSoft<=ambient){
							//brightnessShadows = ambient;
						}
						else if (lightSoft>=brightness){
							brightnessShadows = brightness;
						}
						else{
							brightnessShadows = lightSoft;
						}


						//print brightness
						//std::cout << "brightness:  " << brightnessShadows << std::endl;
						Colour co = Colour(closestValidIntersection.intersectedTriangle.colour.red*brightnessShadows,closestValidIntersection.intersectedTriangle.colour.green*brightnessShadows,closestValidIntersection.intersectedTriangle.colour.blue*brightnessShadows);
						uint32_t colour = (255<<24) + (co.red<<16) + (co.green<<8) +co.blue;
						window.setPixelColour(x,y,colour);

					}
				}

				else{
					//draw black
					uint32_t colour = (255<<24) + (0<<16) + (0<<8) +0;
					window.setPixelColour(x,y,colour);

				}
			
				
			}

		}
	}

	// //pt light pixel on canvas
	// CanvasPoint lightPixel = projectVertexOnToCanvasPoint(cameraPosition,2,lightPosition);
	// //print pixel
	// //set pixel to bllue
	// uint32_t colour = (255<<24) + (255<<16) + (255<<8) +255;
	// //within 5 pixels
	// for (int i = -5; i < 5; i++)
	// {
	// 	for (int j = -5; j < 5; j++)
	// 	{
	// 		window.setPixelColour(10,50,colour);
	// 	}
	// }
	// std::cout << "light pixel:  " << lightPixel.x << " " << lightPixel.y << std::endl;
	// //draw light sour sphere on canvas
	// for (glm::vec3 light : lightSphere){
	// 	CanvasPoint lightPixel = projectVertexOnToCanvasPoint(cameraPosition,2,light);
	// 	//print pixel
	// 	//set pixel to bllue
	// 	uint32_t colour = (255<<24) + (0<<16) + (0<<8) +0;
	// 	//within 5 pixels
	// 	for (int i = -5; i < 5; i++)
	// 	{
	// 		for (int j = -5; j < 5; j++)
	// 		{
	// 			window.setPixelColour(lightPixel.x+i,lightPixel.y+j,colour);
	// 		}
	// 	}
	// // 	//std::cout << "light pixel:  " << lightPixel.x << " " << lightPixel.y << std::endl;
	// }
}

//draw wireframe, rasterised and raytraced scene
void draw(DrawingWindow &window) {
	window.clearPixels();

	/*
	// //task 2.3
	// //grey scale
	// std::vector<float> greyscale;
	// greyscale = interpolateSingleFloats(256,0,window.width);
	// for (size_t y = 0; y < window.height; y++) {
	// 	for (size_t x = 0; x < window.width; x++) {
	// 		//mx 256 is white min 0 is black
	// 		//float red = rand() % 256;
	// 		float red = greyscale[x];
	// 		float green = greyscale[x];
	// 		float blue = greyscale[x];
	// 		uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	// 		window.setPixelColour(x, y, colour);
	// 	}
	// }

	// //task 2.5
	// glm::vec3 topLeft(255, 0, 0);        // red 
	// glm::vec3 topRight(0, 0, 255);       // blue 
	// glm::vec3 bottomRight(0, 255, 0);    // green 
	// glm::vec3 bottomLeft(255, 255, 0);   // yellow
	// //left colum
	// std::vector <glm::vec3> leftCol;
	// leftCol = interpolateThreeElementValues(topLeft,bottomLeft,window.height);
	// //right colum
	// std::vector <glm::vec3> rightCol;
	// rightCol = interpolateThreeElementValues(topRight,bottomRight,window.height);
	// //window
	// //for each row 
	// for (float i = 0; i < window.height; i++)
	// {
	// 	//all coumns interpolate eachcrow
	// 	std::vector <glm::vec3> allCol;
	// 	allCol = interpolateThreeElementValues(leftCol[i],rightCol[i],window.width);
	// 	for (float j = 0; j < window.width; j++)
	// 	{
	// 		int32_t colour = (255 << 24) + (int(allCol[j].x) << 16) + (int(allCol[j].y) << 8) + int(allCol[j].z);
	// 		window.setPixelColour(j,i,colour);
	// 	}	
	// }
	*/
	orbit(cameraPosition,M_PI/30);

			
	//task 6.7 render function
	switch (mode)
	{
	case wireframe:
		
		drawWireframe(window,cameraPosition,focalLength,objectModelTriangles);
		break;

	case rasterised:

		drawRasterisedSceneDepth(window,cameraPosition,focalLength,objectModelTriangles);

		break;
		
	case rayTraced:

		drawRayTraceWithShadow(window,cameraPosition,objectModelTriangles,lightSource,ambient,sourceIntensity,specularPower,shading);
		break;
	
	default:
		break;
	}
		
	

}

//cameria animation
void cameraAnimation(DrawingWindow &window){

	//orbit(cameraPosition,M_PI);
	//cameraPosition = {1,1,4.0};
	specularPower = 0.5; //specular intensity
	specularShininess = 1; //specular shininess
	sourceIntensity = 10;
	ambient = 0.2;
	lightSource = {0,0.75,0.5};
	lightRadius = 0.1; //how many
	//cameraPosition= {1.41142, -0.088799, 3.01}	;
	for (int i = 0; i < 260; i++)
	{	
		//std::cout << "i:  " << i << std::endl;
		//save all frames of the loop
		std::stringstream ss;
		ss  <<"frame" << std::setw(5) << std::setfill('0') << i << ".ppm";
		window.savePPM(ss.str());
		draw(window);
		window.renderFrame();
	
		// if (i<45){
		// 	mode = wireframe;
		// }
		// else if (i>44 && i<150){
		// 	mode = rasterised;
		// }
		// else if (i>149 && i<170){
		// 	mode = rayTraced;
		// 	//hard shadow
		// 	if (i>149 && i<155){
		// 		shadow = 1;
		// 		shading = 0;
		// 	}
		// 	//soft shadow
		// 	else if (i>154 && i<160){
		// 		shadow = 2;
		// 		shading = 0;
		// 	}
		// 	//gouraud shading
		// 	else if (i>159 && i<165){
		// 		shading = 1;
		// 	}
		// 	//phong shading
		// 	else if (i>164 && i<170){
		// 		shading = 2;
		// 	}
		// }
		// else{
		// 	mode = wireframe;
		// }
		// print 
		// std::cout << "i:  " << i << std::endl;

		// if (i<100){
		// 	mode = wireframe;
		// }
		// else if(i>78 && i<82){
		// 	mode = rayTraced;
		// 	shadow = 1;
		// 	shading = 0;
		// }
		// else{
		// 	mode = rayTraced;
		// 	shadow = 2;
		// 	shading = 2;
		// }

		// start loop around the circle
		int circleSteps = 25;
		float circleRadius = std::sqrt(2.0);
		float angleStep = M_PI / 2.0 /circleSteps;
		float z = 4.0;
		// Loop to calculate positions
		if (i<100){
			//angle for the current step
			float angle = i * angleStep;
			float x = circleRadius * cos(angle);
			float y = circleRadius* sin(angle);
			z = z-0.01*i;
			cameraPosition = {x, y, z};

			//print camera position
			//std::cout << "camera position:  " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
		}

		//turn on orbit 
		if (i>99 && i<115){		

			//print camera position
			// std::cout << "camera position:  " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
			// //print camera orientation
			// for (int i = 0; i < 3; i++)
			// {
			// 	std::cout << "camera orientation:  " << cameraOrientation[i].x << " " << cameraOrientation[i].y << " " << cameraOrientation[i].z << std::endl;
			// }

			if (cameraPosition.x>0.05){

				cameraPosition.x -=0.1;
			}
			else{
				cameraPosition.x = 0.0;
			}
			if (cameraPosition.y>0.05){
				cameraPosition.y -=0.1;
			}
			else{
				cameraPosition.y = 0.0;
			}
			cameraPosition.z -=0.02;
		}
		if (i>114 && i<150){
			// std::cout << "camera position:  " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
			// //print camera orientation
			// for (int i = 0; i < 3; i++)
			// {
			// 	std::cout << "camera orientation:  " << cameraOrientation[i].x << " " << cameraOrientation[i].y << " " << cameraOrientation[i].z << std::endl;
			// }
			orbitToggleRight = true;
			orbitOther(cameraPosition,M_PI/2);
			cameraPosition.z -=0.01;
			cameraPosition.y +=0.01;
		}
		if (i>149 && i<180){
			orbitToggleRight = false;
			cameraPosition.z -=0.01;
			cameraPosition.x +=0.02;
			cameraPosition.y +=0.005;
			lookAt(cameraPosition,glm::vec3(0,0,0),cameraOrientation);

		}
		if (i>179 && i<210){

			orbitToggleLeft = true;
			orbit(cameraPosition,M_PI/2);
			cameraPosition.z +=0.01;
			cameraPosition.x -=0.01;
			cameraPosition.y -=0.01;
		}
		if (i>209 && i<260){
			lookAt(cameraPosition,glm::vec3(0,0,0),cameraOrientation);
			orbitToggleLeft = false;
			cameraPosition.z +=0.01;

			//print x and y
			//std::cout << "x:  " << cameraPosition.x << " y:  " << cameraPosition.y << std::endl;

			if (cameraPosition.x<-0.01){
				cameraPosition.x +=0.01;
			}
			else if (cameraPosition.x>0.01){
				cameraPosition.x -=0.01;
			}
			else{
				cameraPosition.x = 0.0;
			}
			if (cameraPosition.y<-0.01){
				cameraPosition.y +=0.01;
			}
			else if (cameraPosition.y>0.01){
				cameraPosition.y -=0.01;
			}
			else{
				cameraPosition.y = 0.0;
			}

			

		}


	}

	
}



// Handle user input  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		//task 5.2 move camera
		float theta = M_PI;
		if (event.key.keysym.sym == SDLK_LEFT) //move left
		{
			//translate camera
		leftTranslation(cameraPosition);
			std::cout << "LEFT" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) //move right
		{
			//translate camera
			rightTranslation(cameraPosition);
			std::cout << "RIGHT" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_UP) //move up
		{
			//translate camera
			upTranslation(cameraPosition);
			std::cout << "UP" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_DOWN) //move down
		{
			//translate camera
			downTranslation(cameraPosition);
			std::cout << "DOWN" << std::endl;
		}
		
		else if (event.key.keysym.sym == SDLK_l) //task 3.3 draw stroked triangle
		{
			//random colour
			Colour ranColour(rand()%256,rand()%256,rand()%256);
			drawStrokedTriangle(window,randomStrokedTriangle(),ranColour);
			std::cout << "u" << std::endl;
		} 
		else if (event.key.keysym.sym == SDLK_f)//task 3.4 draw filled triangle
		{
			Colour ranCol(rand()%256,rand()%256,rand()%256);
			drawFilledTriangle(window,randomStrokedTriangle(),ranCol,Colour(255,255,255));
			std::cout << "f" << std::endl;
		}	
		else if (event.key.keysym.sym==SDLK_t) //task 3.5 draw textured triangle
		{

			CanvasPoint v0 = CanvasPoint(160,10);
			v0.texturePoint= TexturePoint(195,5);
			CanvasPoint v1 = CanvasPoint(300,230);
			v1.texturePoint= TexturePoint(395,380);
			CanvasPoint v2 = CanvasPoint(10,150);
			v2.texturePoint=TexturePoint(65,330);
			CanvasTriangle tri=  CanvasTriangle(v0,v1,v2);
			TextureMap texture = TextureMap("texture.ppm");
			Colour col = Colour(255,255,255);
			drawFilledTextureTriangle(window,tri ,texture);
			drawStrokedTriangle(window,tri,col);
			std::cout << "t" << std::endl;
		}
		
		
		else if (event.key.keysym.sym == SDLK_w) //wireframe
		{
			mode = wireframe;
			std::cout << "wireframe" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_e) //rasterised
		{
			mode = rasterised;
			std::cout << "rasterised" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_r) //ray traced
		{
			mode = rayTraced;
			std::cout << "ray traced" << std::endl;
		}

		else if (event.key.keysym.sym == SDLK_1) //rotate camera up on x axis
		{
			//rotate camera
			xRotateUp(cameraPosition,theta);
			std::cout << "X rotation up" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_2) //rotate camera down on x axis
		{
			//rotate camera
			xRotateDown(cameraPosition,theta);
			std::cout << "X rotation down" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_3) //rotate camera left on y axis
		{ 
			//rotate camera
			yRotateLeft(cameraPosition,theta);
			std::cout << "Y rotation left" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_4) //rotate camera right on y axis
		{
			//rotate camera
			yRotateRight(cameraPosition,theta);
			std::cout << "Y rotation right" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_5) //rotate camera clockwiseon z axis
		{
			//rotate camera
			zRotateRight(cameraPosition,theta);
			std::cout << "Z rotation" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_6) //rotate camera anticlockwise on z axis
		{
			//rotate camera
			zRotateLeft(cameraPosition,theta);
			std::cout << "Z rotation" << std::endl;
		}
		
		else if (event.key.keysym.sym == SDLK_x) //tilt camera
		{
			//rotate camera
			tiltup(cameraOrientation,theta);
			std::cout << "tilting" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_c) //tilt camera
		{
			//rotate camera
			tiltdown(cameraOrientation,theta);
			std::cout << "tilting" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_y) //pan camera
		{
			//rotate camera
			panleft(cameraOrientation,theta);
			std::cout << "panning" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_u) //pan camera
		{
			//rotate camera
			panright(cameraOrientation,theta);
			std::cout << "panning" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_o) //orbit camera
		{
			//toggle orbit on off
			orbitToggleLeft = !orbitToggleLeft;
			if (orbitToggleLeft){
				std::cout << "orbiting" << std::endl;
			}
			else{
				std::cout << "not orbiting" << std::endl;
			}
		}
		//look at camera
		else if (event.key.keysym.sym == SDLK_a) //look at camera
		{
			//look at camera
			lookAt(cameraPosition,glm::vec3(0.0,0.0,0.0),cameraOrientation);
			std::cout << "looking at" << std::endl;
		}


		else if (event.key.keysym.sym == SDLK_n) //reset camera position
		{
			//reset camera position
			cameraPosition = glm::vec3(0.0,0.0,4.0);
			cameraOrientation = {glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,1.0,0.0),glm::vec3(0.0,0.0,1.0)};
			std::cout << "RESET" << std::endl;
		}


	} else if (event.type == SDL_MOUSEBUTTONDOWN) {

		//debugg 
		int mouseX, mouseY;
		SDL_GetMouseState(&mouseX, &mouseY);
		std::cout << "MOUSE CLICK (" << mouseX << ", " << mouseY << ")" << std::endl;
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");

		//print camera position
		std::cout << "camera position:  " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
		//print camera orientation
		std::cout << "camera orientation:  " << cameraOrientation[0].x << " " << cameraOrientation[0].y << " " << cameraOrientation[0].z << std::endl;
	}

}



int main(int argc, char *argv[]) {
	
	
	DrawingWindow window = DrawingWindow(WIDTH,HEIGHT,false);
	SDL_Event event;
	/*

	// //test task2.2
	// std::vector<float> result;
	// result = interpolateSingleFloats(2.2, 8.5, 7);
	// //size_t unsigned int represent vecto
	// for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	// //char output and flushed newline
	// std::cout << std::endl;

	// //test task2.4
	// //test three element interpolation
	// glm::vec3 from(1.0, 4.0, 9.2);
	// glm::vec3 to(4.0, 1.0, 9.8);
	// std::vector<glm::vec3> result1;
	// result1 = interpolateThreeElementValues(from,to,4);
	// for(size_t i=0; i<result1.size(); i++) 
	// {
	// 	glm::vec3& vec=result1[i];
	// 	std::cout << "(" <<vec.x << ","<<vec.y<<","<<vec.z << " )";
	// 	std::cout << std::endl;
	// }

	// //test task 3.2
	// Colour colour = Colour(255,255,255);
	// CanvasPoint from (0,0);
	// CanvasPoint to (round(window.width/2),round(window.height/2));
	// CanvasPoint from1 (window.width-1,0);
	// CanvasPoint to1 (round(window.width/2),round(window.height/2));
	// CanvasPoint from2 (round(window.width/2),0);
	// CanvasPoint to2 (round(window.width/2),window.height-1);
	// CanvasPoint from3 (round(window.width/3),round(window.height/2));
	// CanvasPoint to3 (round(window.width*2/3),round(window.height/2));
	
	// //test task 4.3
	// //load material file
	// std::string materialFile = "cornell-box.mtl";
	// std::vector<Materials> materials = LoadMaterial::loadMaterialFile(materialFile);
	// // for (Materials material : materials)
	// // {
	// // 	std::cout << material.name << " " << material.colour << std::endl;
	// // }
	// //test task 4.2
	// //load object file
	// std::string filename = "cornell-box.obj";
	// std::vector<ModelTriangle> ModelTriangles = LoadObject::loadObjectFile(filename,materials);
	// //print to cout
	// // for (ModelTriangle triangle : ModelTriangles)
	// // {
	// // 	//if white print
	// // 	if (triangle.colour.red==255 && triangle.colour.green==255 && triangle.colour.blue==255){
	// // 		std::cout << triangle << std::endl;
	// // 	}
	// // }

	// // //task 4.6
	// //drawPointCloud(window,cameraPosition,focalLength,ModelTriangles);
	// // //task 4.7
	// // drawWireframe(window,cameraPosition,focalLength,ModelTriangles);
	// // //task 4.8
	// //drawRasterisedScene(window,cameraPosition,focalLength,ModelTriangles);
	// // //task 4.9
	//drawRasterisedSceneDepth(window,cameraPosition,focalLength,ModelTriangles);

	//week 6
	//drawRayTrace(window);
	//drawRayTraceWithShadow(window);
	//drawRasterisedScene(window,cameraPosition,focalLength,ModelTriangles);
	// drawRasterisedSceneDepth(window,cameraPosition,focalLength,ModelTriangles);
	// drawRayTrace(window,cameraPosition,focalLength,ModelTriangles,lightSource);
	// drawRayTraceWithShadow(window,cameraPosition,focalLength,ModelTriangles,lightSource);
	*/

	texture = TextureMap("texture.ppm");
	std::string materialFile = "cornell-box.mtl";
	std::vector<Materials> materials = LoadMaterial::loadMaterialFile(materialFile,texture);

	std::string filename = "cornell-box.obj";
	std::vector<ModelTriangle> box= LoadObject::loadObjectFile(filename,materials,0.5);

	std::string filename1 = "sphere.obj";	
	std::vector<ModelTriangle> sphere= LoadObject::loadObjectFile(filename1,materials,0.3);	

	//combine two objects
	objectModelTriangles.insert(objectModelTriangles.end(),box.begin(),box.end());
	objectModelTriangles.insert(objectModelTriangles.end(),sphere.begin(),sphere.end());

	/*int max = 1000000;
	bool isSphere = false; //false for triangle true for sphere this is for shadows appearing close to the sphere
	float specularPower = 1; //256 for sphrere
	float sourceIntensity = 10;
	float ambient = 0.2;
	int shading = 0; //0 flat 1 gouraud 2 phong
	float imageScalingPlane;
	bool isTexture=false; //false for no texture true for texture
	*/
	//test all draw funnctions
	//drawPointCloud(window,cameraPosition,focalLength,objectModelTriangles);
	//drawWireframe(window,cameraPosition,focalLength,objectModelTriangles);
	//drawRasterisedScene(window,cameraPosition,focalLength,objectModelTriangles);
	//drawRasterisedSceneDepth(window,cameraPosition,focalLength,objectModelTriangles);
	//drawRayTrace(window,cameraPosition,objectModelTriangles,lightSource,ambient,sourceIntensity,specularPower,shading,isSphere);
	//drawRayTraceWithShadow(window,cameraPosition,objectModelTriangles,lightSource,ambient,sourceIntensity,specularPower,shading,isSphere);
	// int remderrrrrr = 5;
	// cameraAnimation(window);
	// if (remderrrrrr==0){
	// 	mode = wireframe;
	// 	renderModeInt = 0;
	// 	cameraAnimation(window);
	// 	std::string ffmpeg_command = "ffmpeg -framerate 15 -i 0frame%05d.ppm -c:v libx264 -pix_fmt yuv420p outputwire.mp4";
	// 	int ret_code = system(ffmpeg_command.c_str());
	// 	if (ret_code == 0) {
	// 		std::cout << "Video created successfully as output.mp4" << std::endl;
	// 	} else {
	// 		std::cerr << "FFmpeg command failed with code: " << ret_code << std::endl;
	// 	}
	// 	remderrrrrr++;
	// }
	// if (remderrrrrr==1){
	// 	mode = rasterised;
	// 	renderModeInt = 1;
	// 	cameraAnimation(window);
	// 	std::string ffmpeg_command = "ffmpeg -framerate 15 -i 1frame%05d.ppm -c:v libx264 -pix_fmt yuv420p outputraster.mp4";
	// 	int ret_code = system(ffmpeg_command.c_str());
	// 	if (ret_code == 0) {
	// 		std::cout << "Video created successfully as outputrasterised.mp4" << std::endl;
	// 	} else {
	// 		std::cerr << "FFmpeg command failed with code: " << ret_code << std::endl;
	// 	}
	// 	remderrrrrr++;
	// }
	// if (remderrrrrr==2){
	// 	mode = rayTraced;
	// 	shadow = 1;
	// 	shading = 0;
	// 	renderModeInt = 2;
	// 	cameraAnimation(window);
	// 	std::string ffmpeg_command = "ffmpeg -framerate 15 -i 2frame%05d.ppm -c:v libx264 -pix_fmt yuv420p outputraytrace1.mp4";
	// 	int ret_code = system(ffmpeg_command.c_str());
	// 	if (ret_code == 0) {
	// 		std::cout << "Video created successfully as outputrasterised.mp4" << std::endl;
	// 	} else {
	// 		std::cerr << "FFmpeg command failed with code: " << ret_code << std::endl;
	// 	}
	// 	remderrrrrr++;
	// }
	// if (remderrrrrr==3){
	// 	mode = rayTraced;
	// 	shadow = 2;
	// 	shading = 0;
	// 	renderModeInt = 3;
	// 	cameraAnimation(window);
	// 	std::string ffmpeg_command = "ffmpeg -framerate 15 -i 3frame%05d.ppm -c:v libx264 -pix_fmt yuv420p outputraytrace2.mp4";
	// 	int ret_code = system(ffmpeg_command.c_str());
	// 	if (ret_code == 0) {
	// 		std::cout << "Video created successfully as outputrasterised.mp4" << std::endl;
	// 	} else {
	// 		std::cerr << "FFmpeg command failed with code: " << ret_code << std::endl;
	// 	}
	// 	remderrrrrr++;
	// }
	// if (remderrrrrr==4){
	// 	mode = rayTraced;
	// 	shadow = 2;
	// 	shading = 1;
	// 	renderModeInt = 4;
	// 	cameraAnimation(window);
	// 	std::string ffmpeg_command = "ffmpeg -framerate 15 -i 4frame%05d.ppm -c:v libx264 -pix_fmt yuv420p outputraytrace3.mp4";
	// 	int ret_code = system(ffmpeg_command.c_str());
	// 	if (ret_code == 0) {
	// 		std::cout << "Video created successfully as outputrasterised.mp4" << std::endl;
	// 	} else {
	// 		std::cerr << "FFmpeg command failed with code: " << ret_code << std::endl;
	// 	}
	// 	remderrrrrr++;
	// }
	// if (remderrrrrr==5){
	// 	mode = rayTraced;
	// 	shadow = 2;
	// 	shading = 2;
	// 	renderModeInt = 5;
	// 	remderrrrrr++;
	// }


	//settings for raytrace
	// specularPower = 0.5; //specular intensity
	// specularShininess = 1; //specular shininess
	// sourceIntensity = 10;
	// ambient = 0.2;
	// shading = 2; //0 flat 1 gouraud 2 phong
	// shadow = 2; //1 hard shadow 2 soft shadow
	// lightSource = {0,0.75,0.5};
	// lightRadius = 0.1; //how many


	while (true) {
		// We MUST poll for events - otherwise the window will freeze	
		if (window.pollForInputEvents(event)){		
			handleEvent(event, window);
			//std::cout << "camera position:  " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
		}
		/* //test week2 
		// draw(window);
		// //test week3 
		// drawLine(window,from,to,colour);
		// drawLine(window,from1,to1,colour);
		// drawLine(window,from2,to2,colour);
		// drawLine(window,from3,to3,colour);		
		// test week 5
		//draw(window);		
		// test week 6
		//drawRayTrace(window);
		//drawRayTraceWithShadow(window);*/		
		mode = rasterised;
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}



}
