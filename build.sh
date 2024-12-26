mkdir -p doc  # To store the documentation  
# Creation of the documentation
cd ./src
doxygen config
cd ../doc
rm -r latex # We keep the HTML version
cd ..
# Creation of the executable
mkdir -p images # To store the images 
mkdir -p files  # To store the text files
mkdir -p build
cd ./build
cmake ..
make
 
# The build.sh creates a doc folder in the folder of the project, runs doxygen config found in the src folder, removes the latex directory within the doc folder, keeping only the HTML directory, creates the images, files and build directories, then runs the CMakeLists.txt. After that, it compiles the project to create the excutable.
