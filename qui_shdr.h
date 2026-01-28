static int shdr_mk(char  const *vshsrc, char  const *fshsrc) {
        char log[4096];
        int po = 0, vso = 0, fso = 0, stts = 0;

        vso = glCreateShader(GL_VERTEX_SHADER);
        if (!vso)
                goto end;

        glShaderSource(vso, 1, &vshsrc, NULL);
        glCompileShader(vso);

        fso = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fso)
                goto end;

        glShaderSource(fso, 1, (char const * const*)&fshsrc, NULL);
        glCompileShader(fso);

        glGetShaderiv(vso, GL_COMPILE_STATUS, &stts);
        if (!stts) {
                glGetShaderInfoLog(vso, 4096, NULL, log);
                printf("vsh: %d\n", stts);
                puts(log);
                goto end;
        }

        glGetShaderiv(fso, GL_COMPILE_STATUS, &stts);
        if (!stts) {
                glGetShaderInfoLog(fso, 4096, NULL, log);
                printf("fsh: %d\n", stts);
                puts(log);
                goto end;
        }

        po = glCreateProgram();
        if (!po)
                goto end;

        glAttachShader(po, vso);
        glAttachShader(po, fso);
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

        if (vso)
                glDeleteShader(vso);

        if (fso)
                glDeleteShader(fso);

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


struct qui_shdr {
	int po;

	/* uniform locations */
	int M;
};

int qui_shdr_mk(struct qui_shdr *qs) {
	if (!qs)
		return -1;

	qs->po = shdr_mk(qui_vsh, qui_fsh);

	if (!qs->po)
		return -1;

	qs->M = glGetUniformLocation(qs->po, "M");

	if (qs->M == -1)
		return -1;

	return 0;
}
