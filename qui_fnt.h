#define QUI_FNT_MX 0x110000

struct qui_fnt {
	int vao;
	int vbo;
	int ebo;

	struct {
		int e0;
		int en;
		float lsb;
		float rsb;
		float xadv;
	} glph[QUI_FNT_MX];	/* michal@todo: we can reduce size by making it hashtable */
};

static int qui_fnt_ld(struct qui_fnt *fnt, char *pth /* .obj wavefront! */);
static int qui_fnt(struct qui_ctx *qs, struct qui_fnt *fnt);

static int qui_fnt_ld(struct qui_fnt *fnt, char *pth /* .obj wavefront! */) {
	char ln[4096];
	FILE *f;
	float *v = NULL;
	int *e = NULL;
	int vn = 0, en = 0;
	int vs = 0, es = 0;
	int g = 0, err = 0;
	char *a;

	if (!fnt)
		return -1;

	memset(fnt, 0, sizeof(struct qui_fnt));

	f = fopen(pth, "r");

	if (!f)
		return -1;

	while (fgets(ln, 4096, f)) {
		switch(ln[0]) {
		case 'o':
			fnt->glph[g].en = en - fnt->glph[g].e0;
			if (a = strstr(ln, "U+")) {
				g = strtol(a+2, NULL, 16);
				if (QUI_FNT_MX < g) {
					printf("err a\n");
					err = -1;
					break;
				}
			} else {
				err = -1;
				break;
			}

			fnt->glph[g].e0 = en;

			if (a = strstr(ln, "lsb")) {
				sscanf(a, "lsb %f", &fnt->glph[g].lsb);
			} else {
					printf("err b\n");
				err = -1;
				break;
			}

			if (a = strstr(ln, "rsb")) {
				sscanf(a, "rsb %f", &fnt->glph[g].rsb);
			} else {
					printf("err c\n");
				err = -1;
				break;
			}

			if (a = strstr(ln, "xadv")) {
				sscanf(a, "xadv %f", &fnt->glph[g].xadv);
			} else {
					printf("err d\n");
				err = -1;
				break;
			}

			break;
		case 'v':
			while (vs < vn + 2) {
				vs = vs ? vs * 2 : 4096;
				v = realloc(v, vs * 2 * sizeof(float));
			}
			if (2 == sscanf(ln, "v %f %f", v + vn, v + vn + 1)) {
				vn += 2;
			}
			break;
		case 'f':
			if (es < en + 3) {
				es = es ? es * 2 : 4096;
				e = realloc(e, es * 3 * sizeof(int));
			}
			if (3 == sscanf(ln, "f %d %d %d", e + en, e + en + 1, e + en + 2)) {
				--e[en++];
				--e[en++];
				--e[en++];
			}
			break;
		case '#':
		default:
			continue;
		};

		if (err)
			break;
	}

	fnt->glph[g].en = 3 * en - fnt->glph[g].e0;

	fclose(f);

	if (err) {
					printf("err e\n");
		free(v);
		free(e);

		return err;
	}

	glCreateBuffers(1, &fnt->vbo);
	glCreateBuffers(1, &fnt->ebo);

	glNamedBufferStorage(fnt->vbo, sizeof(float) * vn, v, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
	glNamedBufferStorage(fnt->ebo, sizeof(int) * en, e, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

	glGenVertexArrays(1, &fnt->vao);

	glBindVertexArray(fnt->vao);
	glBindBuffer(GL_ARRAY_BUFFER, fnt->vbo);
	glVertexAttribPointer(0, 2, GL_FLOAT, 0, 8, 0);
	glEnableVertexAttribArray(0); 
//	glVertexAttrib3f(1, 1.f, 0.5f, 0.f);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, fnt->ebo);
	glBindVertexArray(0);

	free(v);
	free(e);

	return 0;
}


static int qui_fnt(struct qui_ctx *qc, struct qui_fnt *fnt) {
	if (qc) {
		qc->fnt = fnt;
		return 0;
	}
	return -1;
}