#include "LoadMaterial.h"
#include <fstream>
#include <vector>
#include <string>
#include "Utils.h"
#include <glm/glm.hpp>

std::vector<Materials> LoadMaterial::loadMaterialFile(std::string &filename, TextureMap &texture)

{
    //read from file
    std::ifstream inputStream(filename);
    std::string nextLine;

    //store material
    std::string colourName;
    glm::vec3 colour;
    
    //material vector
    std::vector<Materials> allMaterial;
    Materials material;


    while (std::getline(inputStream, nextLine))
    {

        //tokenise line
        std::vector<std::string> tokens = split(nextLine, ' ');

        //print to cout
        //std::cout << nextLine << std::endl;


        //check for new material
        if (tokens[0] == "newmtl")
        {
            //store material name
            colourName = tokens[1];

        }
        else if (tokens[0] == "Kd")
        {
            //store material colour
            //covert int to colour
            float r = std::stof(tokens[1])*255;
            float g = std::stof(tokens[2])*255;
            float b = std::stof(tokens[3])*255;

            colour = glm::vec3(r,g,b);
            //print to cout
            //std::cout << "Kd: " << colour.x << " " << colour.y << " " << colour.z << std::endl;

            material.colour = Colour(colour.x,colour.y,colour.z);
            material.name = colourName;
            material.isColour = true;
            material.isTexture = false;
            allMaterial.push_back(material);
            

        }
        else if (tokens[0] == "map_Kd")
        {
            //store file name in colours
            material.textureMap = texture;
            material.name = colourName;
            material.isColour = false;
            material.isTexture = true;
            allMaterial.push_back(material);
            
        }
        
    }
    //close
    inputStream.close();

    // //print vector
    // for (Materials material : material)
    // {
    //     std::cout << material.name << " " << material.colour << std::endl;
    // }
    
    return allMaterial;

}
