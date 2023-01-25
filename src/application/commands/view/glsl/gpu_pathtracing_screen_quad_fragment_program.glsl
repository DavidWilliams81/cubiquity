#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture;

void main()
{
	ivec2 texSize = textureSize(screenTexture, 0);
	vec2 offset = vec2(1.0 / texSize.x, 1.0 / texSize.y);
	float count = 1.0;
    vec4 col = texture(screenTexture, TexCoords);
	
	vec4 sample = texture(screenTexture, TexCoords + vec2(-offset.x, -offset.y));
	if(abs(col.a - sample.a) < 0.5) { col.rgb += sample.rgb; count += 1.0;}
	
	sample = texture(screenTexture, TexCoords + vec2(-offset.x, +offset.y));
	if(abs(col.a - sample.a) < 0.5) { col.rgb += sample.rgb; count += 1.0;}
	
	sample = texture(screenTexture, TexCoords + vec2(+offset.x, -offset.y));
	if(abs(col.a - sample.a) < 0.5) { col.rgb += sample.rgb; count += 1.0;}
	
	sample = texture(screenTexture, TexCoords + vec2(+offset.x, +offset.y));
	if(abs(col.a - sample.a) < 0.5) { col.rgb += sample.rgb; count += 1.0;}
	
	FragColor = vec4(col.rgb / count, col.a);
}
