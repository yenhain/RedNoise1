//class for materials
#pragma once
#include "Colour.h"
#include "TextureMap.h"


struct Materials {

        std::string name{};
        Colour colour{};
        TextureMap textureMap{};

       // Materials(const std::string name, Colour colour, TextureMap textureMap);
       bool isColour = false;
       bool isTexture = false;
       //bool isMirror = false;


};