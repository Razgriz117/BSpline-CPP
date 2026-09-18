cmake -S . -B build && cmake --build build -j6
cd build && ./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]' ; cd -