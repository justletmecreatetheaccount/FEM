cd genmesh/build
cmake ..
make
./genmesh
cd ../../data
python3 fixMesh.py
cd ..