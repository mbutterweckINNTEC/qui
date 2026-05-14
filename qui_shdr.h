#ifndef QUI_SHDR_H
#define QUI_SHDR_H

extern int qui_shdr_po, qui_shdr_M;
extern int qui_shdr_qr_po[3], qui_shdr_qr_M[3], qui_shdr_qr_o[3], qui_shdr_qr_r[3], qui_shdr_qr_R[3];

int qui_shdr_mk();
int qui_shdr_rm();

#ifdef QUI_IMPL

int qui_shdr_po, qui_shdr_M;
int qui_shdr_qr_po[3], qui_shdr_qr_M[3], qui_shdr_qr_o[3], qui_shdr_qr_r[3], qui_shdr_qr_R[3];

static int shdr_bld(int n, char  const *shsrc[], int shtyp[], char *tfvar) {
	char log[4096];
	int po = 0, so[n] = {}, stts = 0;

	if (!(po = glCreateProgram()))
		goto end;

	for (int i = 0; i < n; ++i) {
		if (!(so[i] = glCreateShader(shtyp[i])))
			goto end;

		glShaderSource(so[i], 1, &shsrc[i], NULL);
		glCompileShader(so[i]);

		glGetShaderiv(so[i], GL_COMPILE_STATUS, &stts);
		if (!stts) {
			glGetShaderInfoLog(so[i], 4096, NULL, log);
			fprintf(stderr, "vsh: %d\n", stts);
			puts(log);
			goto end;
		}

		glAttachShader(po, so[i]);
	}

	if (tfvar)
		glTransformFeedbackVaryings(po, 1, (void*)&tfvar, GL_INTERLEAVED_ATTRIBS);

	glLinkProgram(po);

	glGetProgramiv(po, GL_LINK_STATUS, &stts);
	if (!stts) {
		glGetProgramInfoLog(po, 4096, NULL, log);
		printf("po: %d\n", stts);
		puts(log);
		goto end;
	}
end:
	if (!stts || 0 == po) {
		glDeleteProgram(po);
		po = 0;
	}

	for (int i = 0; i < n; ++i) {
		if (so) {
			glDeleteShader(so[i]);
		}
	}

        return po;
}

char *qui_vsh =
	"#version 440"				"\n"
	"in      vec4 v;"			"\n"
	"in      vec4 c;"			"\n"
	"uniform mat4 M;"			"\n"
	"out     vec4 k;"			"\n"
	""					"\n"
	"void main() {"				"\n"
	"	gl_Position = M * v;"		"\n"
	"	k = c;"				"\n"
	"}"					"\n";

char *qui_fsh =
	"#version 440"				"\n"
	"in  vec4 k;"				"\n"
	"out vec4 K;"				"\n"
	""					"\n"
	"void main() {"				"\n"
	"	K = k;"				"\n"
	"}"					"\n";

char *qui_vsh_qr = 
	"#version 440"				"\n"
	"in      vec4 v;"			"\n"
	"uniform mat4 M;"			"\n"
	""					"\n"
	"void main() {"				"\n"
	"	gl_Position = M * v;"		"\n"
	"}"					"\n";

char *qui_gsh_qr[3] = {
	"#version 440 core"					"\n"
	"layout(points) in;"					"\n"
	"layout(points, max_vertices = 1) out;"			"\n"
	""							"\n"
	"uniform vec3 o, r;"					"\n"
	"uniform float R;"					"\n"
	""							"\n"
	"out vec3 xv;"						"\n"
	""							"\n"
	"uniform mat4 M;"					"\n"
	"void main() {"						"\n"
	"	vec3 p = gl_in[0].gl_Position.xyz /"		"\n"
	"		gl_in[0].gl_Position.w;"		"\n"
	"	vec3 s = p - o;"				"\n"
	"	float l = dot(s, r);"				"\n"
	"	if (0 < l) {"					"\n"
	"		float d = length(s - l * r);"		"\n"
	"		if (d < R) {"				"\n"
	"			gl_Position = vec4(p, 1.0);"	"\n"
	"			xv = p;"			"\n"
	"			EmitVertex();"			"\n"
	"			EndPrimitive();"		"\n"
	"		}"					"\n"
	"	}"						"\n"
	"}"							"\n",

	"#version 440 core"					"\n"
	"layout(lines) in;"					"\n"
	"layout(line_strip, max_vertices = 2) out;"		"\n"
	""							"\n"
	"uniform vec3 o, r;"					"\n"
	"uniform float R;"					"\n"
	""							"\n"
	"out vec3 xv;"						"\n"
	""							"\n"
	"uniform mat4 M;"					"\n"
	"void main() {"						"\n"
	"	vec3 pa = gl_in[0].gl_Position.xyz /"		"\n"
	"		gl_in[0].gl_Position.w;"		"\n"
	"	vec3 pb = gl_in[1].gl_Position.xyz /"		"\n"
	"		gl_in[1].gl_Position.w;"		"\n"
	"	vec3 b = pb - pa;"				"\n"
	"	vec3 c = pa - o;"				"\n"
	"	vec3 rxb = cross(r, b);"			"\n"
	"	float d = length(dot(c, rxb)) / length(rxb);"	"\n"
	"	if (d < R) {"					"\n"
	"		vec3 rxrxb = cross(r, rxb);"		"\n"
	"		float s = -dot(c, rxrxb) / dot(b, rxrxb);\n"
	"		if (0.0 <= s && s <= 1.0) {"		"\n"
	"			gl_Position = vec4(pa, 1.0);"	"\n"
	"			xv = pa;"			"\n"
	"			EmitVertex();"			"\n"
	"			gl_Position = vec4(pb, 1.0);"	"\n"
	"			xv = pb;"			"\n"
	"			EmitVertex();"			"\n"
	"			EndPrimitive();"		"\n"
	"		}"					"\n"
	"	}"						"\n"
	"}"							"\n",
};

int qui_shdr_mk() {
	qui_shdr_po = shdr_bld(2, (char const *[]){ qui_vsh, qui_fsh}, (int[]){ GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, NULL);

	if (!qui_shdr_po)
		return -1;

	qui_shdr_M = glGetUniformLocation(qui_shdr_po, "M");

	if (qui_shdr_M == -1)
		return -1;

	for (int i = 0; i < 3; ++i) {
		qui_shdr_qr_po[i] = shdr_bld(2, 
				(char const *[]){ qui_vsh_qr, qui_gsh_qr[i]}, 
				(int[]){ GL_VERTEX_SHADER, GL_GEOMETRY_SHADER}, "xv"
		);

		if (!qui_shdr_qr_po[i])
			return -1;

		if (-1 == (qui_shdr_qr_M[i] = glGetUniformLocation(qui_shdr_qr_po[i], "M")))
			return -1;

		if (-1 == (qui_shdr_qr_o[i] = glGetUniformLocation(qui_shdr_qr_po[i], "o")))
			return -1;

		if (-1 == (qui_shdr_qr_r[i] = glGetUniformLocation(qui_shdr_qr_po[i], "r")))
			return -1;

		if (-1 == (qui_shdr_qr_R[i] = glGetUniformLocation(qui_shdr_qr_po[i], "R")))
			return -1;
	}

	return 0;
}

int qui_shdr_rm() {
	glDeleteProgram(qui_shdr_po);

	qui_shdr_po = 0;
	qui_shdr_M = 0;

	return 0;
}

#endif /* QUI_IMPL */
#endif /* QUI_SHDR_H */
