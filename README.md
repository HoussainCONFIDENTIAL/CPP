# MNT

C++ project - FISE(2) ROB 2026

### Author
Hussein Rammal <hussein.rammal@ensta-bretagne.org> 
### Description
Computer project aimed at generating a marine background image from geographical data. This is a link to a fully descriptive video for the project: <https://youtu.be/Fhz4EDRbnbo>
### Explanation video [uploaded on YouTube]().

## Table of Contents
1. [Instructions](#instructions)
2. [How it works](#how-it-works)
3. [Work done](#work-done)
4. [Additional information](#additional-information)

## Instructions
- Open the main folder and run the build.sh file (using ./build.sh) which generates the documentation and the executable file  create_raster
- Move Guerledan_Feb19_50cm_wgs84.txt to the  ./files/ folder
- Go to a terminal in the ./build/ folder (using cd ./build/), and type: ./create_raster ../files/Guerledan_Feb19_50cm_wgs84.txt <number of pixels> 
- The created image will be found in the ./images/ folder
  
Summary diagram

   graph TD;
    build.sh --> ./files;
    build.sh --> ./doc;
    build.sh --> ./build;
    ./files --> fichiers.txt;
    ./doc --> documentation;
    ./build --> create_raster;
    
## How it works
The main steps of the program are:
1. Parsing the text file to retrieve and convert the data
2. Using the [delaunator.h](https://github.com/delfrrr/delaunator-cpp) file to perform a Delaunay triangulation with the given points
3. Iterating through the triangles to store them in a 32x32 grid (indexed by two integers)
4. Creating a color scale and the binary .ppm file
5. Iterating through the pixels
   1. Narrowing down the region it belongs to using binary search
   2. Searching among the triangles in the region for the one it belongs to
   3. Calculating the corresponding depth
   4. Converting to color and writing to the file

    graph TB;
     file.txt --> class_Point --conversion with proj--> vecteur_v -->delaunator --> quadrillage;
     class_Point --> map_depth;
     pixels --binary search--> quadrillage --> triangle --> map_depth --> color -->image.ppm;
     triangle --> shading --> image.ppm;
   
## Work done
Summary of the work done compared to the project requirements:
- The generated image is in binary (.ppm P6) format and in color
- The compilation is done with cmake
- The coordinates in the files are converted to Lambert93 (see proj parameters)
- The code is documented and I hope it's understandable
- An .html documentation file is generated
"# CPP" 
