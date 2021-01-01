#version 150
uniform sampler2D myTexture;
in vec2 position;
in vec2 vertTexCoord;
out vec2 fragTexCoord;
uniform mat4 mvp;
uniform vec4 flip; 

void main() {
	gl_Position = mvp * vec4(position, 0, 1);
// 	fragTexCoord = vec2(vertTexCoord.x , vertTexCoord.y );
	fragTexCoord = vec2(vertTexCoord.x * flip.x, vertTexCoord.y * flip.y);
}
