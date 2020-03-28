#version 330 core

in vec4 exCol;
out vec4 FragColor;

void main(){
   FragColor = exCol;
   // FragColor = vec4(1.0f, 0.5f, 0.5f, 1.0f);
}


