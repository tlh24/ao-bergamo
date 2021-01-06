#version 150
in vec2 position;

uniform mat4 mvp;
uniform vec4 myColor; 
out vec4 fragColor;

void main() {
	//convert from pixel coordinates. 
	vec4 pos;
	pos.x = position.x/1024.0 - 1.0; 
	pos.y = 1.0 - position.y/1024.0; //flip to correspond to GL coords
	pos.z = 0.0; 
	pos.w = 1.0; 
	gl_Position = mvp * pos;
	fragColor = myColor; 
}
