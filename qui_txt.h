#define QUI_TXT_MX 256

int qui_txt_utf8_pop(char **a) {
    int c = (unsigned char)(*a)[0];
    if (c < 0x80) {
		(*a)++;
		return c;
	}
    if ((c & 0xe0) == 0xc0) {
		c = ((c&31)<<6) | ((*a)[1]&0x3f);
		*a += 2;
		return c;
	}
    if ((c & 0xf0) == 0xe0) {
		c = ((c&0xf)<<0xC) | (((*a)[1]&0x3f)<<6) | ((*a)[2]&0x3f);
		*a+=3; return c;
	}
    c = ((c&7)<<0x12) | (((*a)[1]&0x3f)<<12) | (((*a)[2]&0x3f)<<6) | ((*a)[3]&0x3f);
	*a+=4;
    return c;
}

static int qui_txt(struct qui_ctx *qc, char *a, float44_t M) {
	int u;
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


	glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
	glEnable(GL_BLEND);
	glEnable(GL_POLYGON_SMOOTH);

	glUseProgram(qc->po);
	glBindVertexArray(qc->fnt->vao);

	X.m[3][0] += qc->fnt->glph[u].lsb;

	while (*a) {
		u = qui_txt_utf8_pop(&a);

		if (u < QUI_FNT_MX) {
			MX = mul_float44(X, M);
			glUniformMatrix4fv(qc->M, 1, 0, &MX.m[0][0]);

			/* todo:@michal: it may be good to change it to glMultiDrawElements and index uniform array through drawID */
			glDrawElements(GL_TRIANGLES, qc->fnt->glph[u].en, GL_UNSIGNED_INT, (void*)(qc->fnt->glph[u].e0 * sizeof(int)));

			X.m[3][0] += qc->fnt->glph[u].xadv;
		}
	}
	glDrawElements(GL_TRIANGLES, 106*3, GL_UNSIGNED_INT, NULL);
	glDisable(GL_POLYGON_SMOOTH);
	glDisable(GL_BLEND);
}
