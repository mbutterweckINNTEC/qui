#ifndef QUI_QR_H
#define QUI_QR_H

enum {
	QUI_QR_TYP_PNTS,
	QUI_QR_TYP_LNS,
	QUI_QR_TYP_TRIS,
};

int qui_qr_bgn(float3_t r, float3_t o, float R, int sz, int typ);
int qui_qr_psh_arr(float44_t M, int vbo, int vn, int vstrd, int voff);
int qui_qr_end();

extern float3_t *qui_qr[3];
extern int qui_qr_n[3];
static int const qui_qr_s[3] = { 1, 3, 4 };

enum {
	QUI_QR_FND_NRST_PNT = 0x1,
	QUI_QR_FND_NRST_LN = 0x2,
	QUI_QR_FND_NRST_TRI = 0x4,
	QUI_QR_FND_NRST = QUI_QR_FND_NRST_PNT | QUI_QR_FND_NRST_LN | QUI_QR_FND_NRST_TRI,
	QUI_QR_FND_BST_PNT = 0x8 | QUI_QR_FND_NRST_PNT,
	QUI_QR_FND_BST_LN = 0x8 | QUI_QR_FND_NRST_LN,
	QUI_QR_FND_BST_TRI = 0x8 | QUI_QR_FND_NRST_TRI,
	QUI_QR_FND_BST = 0x8 | QUI_QR_FND_NRST,

	/* internal */
	QUI_QR_FND_BST_ANY = 0x8,
};
int qui_qr_find(int fnd, float3_t o, float r, int *typ /*out*/);

#ifdef QUI_IMPL

int qui_qr_typ;
int qui_qr_xo, qui_qr_xbo[3], qui_qr_qo, qui_qr_vao, qui_qr_n[3], qui_qr_sz[3];
float3_t *qui_qr[3];
float3_t qui_qr_r, qui_qr_o;
float qui_qr_R;

int qui_qr_ini() {
	if (!qui_qr_qo)
		glGenQueries(1, &qui_qr_qo);

	if (!qui_qr_qo)
		goto ouch;

	if (!qui_qr_xo)
		glGenTransformFeedbacks(1, &qui_qr_xo);

	if (!qui_qr_xo)
		goto ouch;

	if (!qui_qr_xbo[0])
		glGenBuffers(3, qui_qr_xbo);

	if (!qui_qr_xbo[0])
		goto ouch;

	if (!qui_qr_vao)
		glGenVertexArrays(1, &qui_qr_vao);

	if (!qui_qr_vao)
		goto ouch;

	return 0;
ouch:
	if (qui_qr_xo)
		glDeleteQueries(1, &qui_qr_qo);

	if (qui_qr_xo)
		glDeleteTransformFeedbacks(1, &qui_qr_xo);

	if (qui_qr_xbo[0])
		glDeleteBuffers(3, qui_qr_xbo);

	if (qui_qr_vao)
		glDeleteVertexArrays(1, &qui_qr_vao);

	return -1;

}

int qui_qr_bgn(float3_t r, float3_t o, float R, int sz, int typ) {
	if (qui_qr_ini())
		return -1;

	if (qui_qr[typ]) {
		glUnmapNamedBuffer(qui_qr_xbo[typ]);
	}

	qui_qr_typ = typ;
	qui_qr_r = normal_float3(r);
	qui_qr_o = o;
	qui_qr_R = R;
	qui_qr_n[typ] = 0;

	if (qui_qr_sz[typ] != sz) {
		qui_qr_sz[typ] = sz;
		glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, qui_qr_xbo[typ]);
		glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, qui_qr_sz[typ] * sizeof(float) * 3, NULL, GL_DYNAMIC_READ);
	}

	glBeginQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, qui_qr_qo);

	glBindBufferRange(GL_TRANSFORM_FEEDBACK_BUFFER, 0, qui_qr_xbo[typ], 0, qui_qr_sz[typ] * sizeof(float) * 3);

	glUseProgram(qui_shdr_qr_po[typ]);
	glUniform3f(qui_shdr_qr_o[typ], qui_qr_o.x, qui_qr_o.y, qui_qr_o.z);
	glUniform3f(qui_shdr_qr_r[typ], qui_qr_r.x, qui_qr_r.y, qui_qr_r.z);
	glUniform1f(qui_shdr_qr_R[typ], qui_qr_R);

	glBeginTransformFeedback(GL_POINTS);

	glBindVertexArray(qui_qr_vao);

	return 0;
}

int qui_qr_psh_arr(float44_t M, int vbo, int vn, int vstrd, int voff) {
	glUniformMatrix4fv(qui_shdr_qr_M[qui_qr_typ], 1, 0, &M.m[0][0]);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, 0, vstrd, (char*)NULL + voff);
	glEnableVertexAttribArray(0);
	glEnable(GL_RASTERIZER_DISCARD);

	int gltyp;

	switch (qui_qr_typ) {
	case QUI_QR_TYP_PNTS: gltyp = GL_POINTS; break;
	case QUI_QR_TYP_LNS: gltyp = GL_LINES; break;
	case QUI_QR_TYP_TRIS: gltyp = GL_TRIANGLES; break;
	default:
		return -1;
	};

	glDrawArrays(gltyp, 0, vn);
	glDisable(GL_RASTERIZER_DISCARD);

	return 0;
}

int qui_qr_psh_elm(float44_t M, int vbo, int ebo, int en, int vstrd, int voff) {
	glUniformMatrix4fv(qui_shdr_qr_M[qui_qr_typ], 1, 0, &M.m[0][0]);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, 0, vstrd, (char*)NULL + voff);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

	int gltyp;

	switch (qui_qr_typ) {
	case QUI_QR_TYP_PNTS: gltyp = GL_POINTS; break;
	case QUI_QR_TYP_LNS: gltyp = GL_LINES; break;
	case QUI_QR_TYP_TRIS: gltyp = GL_TRIANGLES; break;
	default:
		return -1;
	};

	glEnable(GL_RASTERIZER_DISCARD);
	glDrawElements(gltyp, en, GL_UNSIGNED_INT, NULL);
	glDisable(GL_RASTERIZER_DISCARD);

	return 0;
}

int qui_qr_end() {
	glBindVertexArray(0);
	glEndTransformFeedback();
	glEndQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN);
	glGetQueryObjectiv(qui_qr_qo, GL_QUERY_RESULT, &qui_qr_n[qui_qr_typ]);
//	qui_qr_n[qui_qr_typ] *= qui_qr_typ + 1;
	qui_qr[qui_qr_typ] = glMapNamedBuffer(qui_qr_xbo[qui_qr_typ], GL_READ_ONLY);

	return 0;
}

int qui_qr_find(int fnd, float3_t o, float R, int *typ /*out*/) {
	float t = FLT_MAX;
	int k = -1;
	int q = 0;

	for (int f = QUI_QR_FND_NRST_PNT; q < 3; f <<= 1, ++q) {
		if (f & fnd) {
			for (int i = 0; i < qui_qr_n[q]; i += qui_qr_s[q]) {
				float l = length_float3(sub_float3(o, qui_qr[q][i + qui_qr_s[q] - 1]));
				if (l < t) {
					if (0 <= k && fnd & QUI_QR_FND_BST_ANY) {
						float r = length_float3(
							sub_float3(
								qui_qr[*typ][k + qui_qr_s[*typ] - 1],
								qui_qr[q][i + qui_qr_s[q] - 1]));
						if (r < R)
							continue;
					}
					t = l;
					k = i;
					*typ = q;
				}
			}
		}
	}

	return k;
}

#endif /* QUI_IMPL */
#endif /* QUI_QR_H */
