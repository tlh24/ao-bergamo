#version 150
in vec2 fragTexCoord;
uniform sampler2D myTexture;
void main() {
	gl_FragColor = texture2D(myTexture, fragTexCoord);
	gl_FragColor.a = 1.0;
	//input is 8-bit grayscale in the red channel only.
	gl_FragColor.g = gl_FragColor.r;
	gl_FragColor.b = gl_FragColor.r;
}
