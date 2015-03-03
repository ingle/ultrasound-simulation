g++ common/phantom.cpp create/createphantom.cpp -I common -I create -o createPhantom
g++ common/phantom.cpp compress/*.cpp -I common -I compress -o compressPhantom
g++ common/phantom.cpp rfData/*.cpp -I common -I rfData -O3 -o rfDataProgram
