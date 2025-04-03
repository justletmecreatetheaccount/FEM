echo -e "\033[34m[INFO] Navigating to the build folder... \033[0m\n"
cd build
echo -e "\033[34m[INFO] Building project... \033[0m\n"
cmake ..
echo -e "\033[34m[INFO] Compiling... \033[0m\n"
make
echo -e "\033[34m[INFO] Running the solver... \033[0m\n"
./solver
cd ..