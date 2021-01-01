#version 150
in vec4 vertex; // <vec2 pos, vec2 tex>
out vec2 TexCoords;

uniform mat4 mvp;

void main()
{
    gl_Position = mvp * vec4(vertex.xy, 0.0, 1.0);
    TexCoords = vertex.zw;
}  

