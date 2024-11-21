//load object file 

//task4.1 read in file cornell-box.obj
//read in vertices and faces
//use ifstream make calls to getline, 
//split function in Utils to tokenise line
//store vertices and faces in vectors
//store object name and colour
//scale vertices
//load material file
//store object colour
//vector of triangles

#include "LoadObject.h"
#include "Utils.h"
#include <fstream>
#include <vector> 
#include <string>
#include <glm/glm.hpp>


std::vector<ModelTriangle> LoadObject::loadObjectFile(std::string &filename,std::vector<Materials> materials, float scaleFactor) {
  
    //read from file
    std::ifstream inputStream(filename);
    std::string nextLine;


    //for storing vertices and faces
    std::vector<ModelTriangle> triangles;
    std::vector<glm::vec3> vertices;
    std::vector<TexturePoint> verticesTexture;
    std::vector<glm::vec3> faces;
    //file also has ojbect name and colour
    std::string objectName;
    Colour objectColours;
    bool isTex = false;

    //group obects into vector
    std::vector<std::vector<ModelTriangle>> objects;

    //read in file
    while (std::getline(inputStream, nextLine)) {


        //tokenise line
        std::vector<std::string> tokens = split(nextLine, ' ');

        //print to cout
        //std::cout << nextLine << std::endl;

        //check if line is mtllib cornell-box.mtl
        if (tokens[0] == "mtllib") 
        {
            //store file header
            std::string mtllib = tokens[1];
            //print to cout
            //std::cout << "mtllib: " << mtllib << std::endl;
        }
        
        //check for object name
        else if (tokens[0] == "o") 
        {
            //store object name
            objectName = tokens[1];
            //print to cout
            //std::cout << "o: " << objectName << std::endl;
            
           
        }     
        //check for object colour
        else if (tokens[0] == "usemtl") {

            //store object colour
            for (Materials material : materials)
            {
                if (tokens[1]==material.name)
                {
                    if (material.isColour)
                    {
                        isTex = false;
                        objectColours.name = material.name;
                        objectColours.red = material.colour.red;
                        objectColours.green = material.colour.green;
                        objectColours.blue = material.colour.blue;

                        
                    }
                    else if (material.isTexture)
                    {
                        //set white if no colour
                        isTex = true;
                        objectColours.name = material.name;
                        objectColours.red = 255;
                        objectColours.green = 255;
                        objectColours.blue = 255;
                    }

                    
                }
            }
            //print to cout
           // std::cout << "usemtl: " << tokens[1] << std::endl;

        }
           
        //find vertices 
        else if (tokens[0] == "v")
        {
            //store vertices as vec3
            glm::vec3 vertex = glm::vec3(std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3]));
            //scale vertices
            vertex *= scaleFactor;
            vertices.push_back(vertex);
            //print to cout
            //std::cout << "Vertices: " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
        }

        //find texture points
        else if (tokens[0] == "vt")
        {
            //store texture points as vec3
            TexturePoint texturePoint = TexturePoint(std::stof(tokens[1]), std::stof(tokens[2]));
            verticesTexture.push_back(texturePoint);
            //print to cout
            //std::cout << "Texture Points: " << texturePoint.x << " " << texturePoint.y << std::endl;
        }


        //find faces
        else if (tokens[0] == "f")
        {


            if (isTex){
                std::vector<std::string> face1 = split(tokens[1], '/');
                std::vector<std::string> face2 = split(tokens[2], '/');
                std::vector<std::string> face3 = split(tokens[3], '/');
                //store faces as vec3
                glm::vec3 vert = glm::vec3(std::stoi(face1[0])-1, std::stoi(face2[0])-1, std::stoi(face3[0])-1);
                glm::vec3 vertT = glm::vec3(std::stoi(face1[1])-1, std::stoi(face2[1])-1, std::stoi(face3[1])-1);

                //print to cout
                //
                // std::cout << "Faces: " << vert.x << " " << vert.y << " " << vert.z << std::endl;
                // std::cout << "Faces Texture: " << vertT.x << " " << vertT.y << " " << vertT.z << std::endl;
                ///TODO FINNISH FUNCTIOn

                //get normal of triangle
                glm::vec3 e0 = vertices[vert.y] - vertices[vert.x];
                glm::vec3 e1 = vertices[vert.z] - vertices[vert.x];
                glm::vec3 normal = glm::cross(e0,e1);
                normal = glm::normalize(normal); //NOTE may be unecessary

                //print normals
                //std::cout << "Normal: " << normal.x << " " << normal.y << " " << normal.z << std::endl;

                //add vertices in glm to triangle
                ModelTriangle triangle = ModelTriangle();
                triangle.vertices[0] = vertices[vert.x];
                triangle.vertices[1] = vertices[vert.y];
                triangle.vertices[2] = vertices[vert.z];
                triangle.normal = normal;
                triangle.object = objectName;
                triangle.texturePoints[0] = verticesTexture[vertT.x];
                triangle.texturePoints[1] = verticesTexture[vertT.y];
                triangle.texturePoints[2] = verticesTexture[vertT.z];
                triangle.isTextured = true;

                // //print verti
                // std::cout << "Vertices: " << triangle.vertices[0].x << " " << triangle.vertices[0].y << " " << triangle.vertices[0].z << std::endl;
                // //prtin texture points
                // std::cout << "Texture Points: " << triangle.texturePoints[0] << " " << triangle.texturePoints[1] << " " << triangle.texturePoints[2] << std::endl;



                //print triangle with object colour
                //std::cout << triangle << std::endl;
                //std::cout << "Colour: " << triangle.colour << std::endl;
                triangles.push_back(triangle);
            }
            else{
                //store faces as vec3
                glm::vec3 vert = glm::vec3(std::stoi(tokens[1])-1, std::stoi(tokens[2])-1, std::stoi(tokens[3])-1); 
                //find normal of triangle
                glm::vec3 e0 = vertices[vert.y] - vertices[vert.x];
                glm::vec3 e1 = vertices[vert.z] - vertices[vert.x];
                glm::vec3 normal = glm::cross(e0,e1);
                normal = glm::normalize(normal); //NOTE may be unecessary

                //add vertices in glm to triangle
                ModelTriangle triangle = ModelTriangle();
                triangle.vertices[0] = vertices[vert.x];
                triangle.vertices[1] = vertices[vert.y];
                triangle.vertices[2] = vertices[vert.z];
                triangle.normal = normal;
                triangle.object = objectName;
                triangle.colour = objectColours;
                triangle.isTextured = false;

                triangles.push_back(triangle);
            }




        }


    }
    //close
    inputStream.close();
    //print to cout
    //print triangles
    // for (ModelTriangle triangle : triangles)
    // {
    //     std::cout << triangle << std::endl;
    // }

    // //print number of triangles
    // std::cout << "Number of triangles: " << triangles.size() << std::endl;
    // //count number of texture points
    // std::cout << "Number of texture points: " << verticesTexture.size() << std::endl;

    return triangles;

}
