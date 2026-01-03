mkdir -p build
cd build

# Configure CMake
cmake ..

# Build the project
cmake --build .

# Run from the root directory so paths work correctly
./timbre example.wav
