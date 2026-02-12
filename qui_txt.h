#define QUI_TXT_MX 256

int qui_txt_ms = 0;

static int qui_txt(char *a, float44_t M, float4_t c) {
	int u, stts = 1;

	enum {
		STTS_STRT = 0x1
	};

	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);

	float44_t X = {
		1.f, 0.f, 0.f, 0.f,
		0.f, 1.f, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		0.f, 0.f, 0.f, 1.f
	};
	float44_t PVMX = identity_sc;

	if (!a || *a == '\0')
		return -1;

	if (!qui_fnt)
		return -1;

	if (qui_flgs & QUI_FLGS_AA) {
		glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
		glEnable(GL_BLEND);
		glEnable(GL_POLYGON_SMOOTH);
	}

	glUseProgram(qui_shdr_po);
	glBindVertexArray(qui_fnt->vao);

	glVertexAttrib4fv(1, (GLfloat*)&c);

	while (*a) {
		u = qui_utf8_pop(&a);

		if (QUI_FNT_MX <= u)
			continue;

		if (stts & STTS_STRT) {
			X.m[3][0] = qui_fnt->glph[u].lsb;
			stts ^= STTS_STRT;
		}

		PVMX = mul_float44(mul_float44(X, M), mul_float44(V, P));
		glUniformMatrix4fv(qui_shdr_M, 1, 0, &PVMX.m[0][0]);

		/* michal@todo: it may be good to change it to glMultiDrawElements and index uniform array through drawID */
		if (qui_fnt->glph[u].en && u != ' ')
			glDrawElements(GL_TRIANGLES, qui_fnt->glph[u].en, GL_UNSIGNED_INT, (void*)(qui_fnt->glph[u].e0 * sizeof(int)));

		X.m[3][0] += qui_fnt->glph[u].xadv;
	}

	if (qui_flgs & QUI_FLGS_AA) {
		glDisable(GL_POLYGON_SMOOTH);
		glDisable(GL_BLEND);
	}

	glBindVertexArray(0);

	return 0;
}
