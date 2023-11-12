#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture;

void main()
{
    vec4 col = texture(screenTexture, TexCoords);
    FragColor = vec4(col.rgba / col.a);
	//FragColor = vec4(col.rgb, 0.0);
} 