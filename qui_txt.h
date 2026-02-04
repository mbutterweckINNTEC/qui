#define QUI_TXT_MX 256

int qui_txt_ms = 0;

static int qui_txt(struct qui_ctx *qc, char *a, float44_t M, float4_t c) {
	int u, stts = 1;

	enum {
		STTS_STRT = 0x1
	};

	float44_t X = {
		1.f, 0.f, 0.f, 0.f,
		0.f, 1.f, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		0.f, 0.f, 0.f, 1.f
	};
	float44_t MX;

	if (!qc || !a || *a == '\0')
		return -1;

	if (!qc->fnt)
		return -1;

	if (qc->flgs & QUI_FLGS_AA) {
		glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
		glEnable(GL_BLEND);
		glEnable(GL_POLYGON_SMOOTH);
	}

	glUseProgram(qc->po);
	glBindVertexArray(qc->fnt->vao);

	glVertexAttrib4fv(1, (GLfloat*)&c);

	while (*a) {
		u = qui_utf8_pop(&a);

		if (QUI_FNT_MX <= u)
			continue;

		if (stts & STTS_STRT) {
			X.m[3][0] = qc->fnt->glph[u].lsb;
			stts ^= STTS_STRT;
		}

		MX = mul_float44(X, M);
		glUniformMatrix4fv(qc->M, 1, 0, &MX.m[0][0]);

		/* michal@todo: it may be good to change it to glMultiDrawElements and index uniform array through drawID */
		glDrawElements(GL_TRIANGLES, qc->fnt->glph[u].en, GL_UNSIGNED_INT, (void*)(qc->fnt->glph[u].e0 * sizeof(int)));

		X.m[3][0] += qc->fnt->glph[u].xadv;
	}

	if (qc->flgs & QUI_FLGS_AA) {
		glDisable(GL_POLYGON_SMOOTH);
		glDisable(GL_BLEND);
	}

	glBindVertexArray(0);

	return 0;
}
