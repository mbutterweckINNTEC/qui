extern int qui_shdr_po, qui_shdr_M;

int qui_shdr_mk();
int qui_shdr_rm();

#ifdef QUI_IMPL

int qui_shdr_po, qui_shdr_M;

static int shdr_bld(char  const *vshsrc, char  const *fshsrc) {
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


int qui_shdr_mk() {
	qui_shdr_po = shdr_bld(qui_vsh, qui_fsh);

	if (!qui_shdr_po)
		return -1;

	qui_shdr_M = glGetUniformLocation(qui_shdr_po, "M");

	if (qui_shdr_M == -1)
		return -1;

	return 0;
}

int qui_shdr_rm() {
	glDeleteProgram(qui_shdr_po);

	qui_shdr_po = 0;
	qui_shdr_M = 0;

	return 0;
}

#endif