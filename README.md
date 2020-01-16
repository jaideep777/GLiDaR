# GLiDaR

Read and Visualize LAS files

## Compiling

0) Install openGL libraries
1) Install libglm from apt - `sudo apt-get install libglm-dev`
2) Compile QuickGL as a shared library - Clone https://github.com/jaideep777/quickGL and type `make`
3) Intall liblas from apt - `sudo apt-get install liblas-dev liblas-c-dev`
4) Edit path to libquickgl.so in `glidar/Makefile`, and set `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/quickgl/lib`
5) compile glidar - `make`

