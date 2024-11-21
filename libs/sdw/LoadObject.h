//header for load object

#include "ModelTriangle.h"
#include <string>
#include <vector>
#include "LoadMaterial.h"
#include "TextureMap.h"
#include "TexturePoint.h"



class LoadObject {
    public:
     
        static std::vector<ModelTriangle> loadObjectFile(std::string &filename, std::vector<Materials> materials,float scaleFactor);
        
};
    


