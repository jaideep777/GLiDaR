#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec4 aCol;
out vec4 exCol;

uniform mat4 transform;
uniform float pointSize=2;

void main(){
   gl_Position = transform*vec4(aPos, 1.0);
   gl_PointSize = pointSize;
   exCol = aCol;

}

