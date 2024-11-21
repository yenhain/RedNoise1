//header to load material
#pragma once

#include "Colour.h"
#include "TextureMap.h"
#include "Materials.h"


class LoadMaterial {
    public:
     

       //load material
       static std::vector<Materials> loadMaterialFile(std::string &filename, TextureMap &texture);


        
};
