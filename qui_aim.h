#ifndef QUI_AIM_H
#define QUI_AIM_H

int qui_aim(float3_t *p, float3_t *n, float3_t *s, float3_t *t);

extern float qui_aim_R;
extern float qui_aim_lnwdth;

/* PRIV */

int qui_aim_mk();
int qui_aim_rm();

#ifdef QUI_IMPL

int qui_aim_vbo, qui_aim_vao;

float qui_aim_R = 0.03125f;
float qui_aim_lnwdth = 1;

#define QUI_AIM_N 32

#define QUI_AIM_O 0.03125f
#define QUI_AIM_L 0.0625f
#define QUI_AIM_R 0.0625f

#define QUI_AIM_X_CLR (float3_t) { 1.0f, 0.5f, 0.5f }
#define QUI_AIM_Y_CLR (float3_t) { 0.5f, 1.0f, 0.5f }
#define QUI_AIM_Z_CLR (float3_t) { 0.5f, 0.5f, 1.0f }
#define QUI_AIM_V_CLR (float3_t) { 1.0f, 1.0f, 1.0f }

#define QUI_AIM_M 12

#define QUI_AIM_SZ ((QUI_AIM_N + QUI_AIM_M) * 2 * sizeof(float3_t))

int qui_aim_mk() {
	/* visual objects */
	glCreateBuffers(1, &qui_aim_vbo);

	if (!qui_aim_vbo)
		goto ouch;

	glNamedBufferStorage(qui_aim_vbo, QUI_AIM_SZ, NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

	float3_t *v = glMapNamedBufferRange(qui_aim_vbo, 0, QUI_AIM_SZ, GL_MAP_WRITE_BIT);

	if (!v)
		goto ouch;

	/* dashed circle */
	for (int i = 0; i < QUI_AIM_N; ++i) {
		if (i & 2) {
			float phi = (float)i / (float)QUI_AIM_N * 2.0 * M_PI;

			*v++ = (float3_t) { QUI_AIM_R * cos(phi), QUI_AIM_R * sin(phi), 0.f };
			*v++ = QUI_AIM_V_CLR;

			phi = (float)(i + 1) / (float)QUI_AIM_N * 2.0 * M_PI;

			*v++ = (float3_t) { QUI_AIM_R * cos(phi), QUI_AIM_R * sin(phi), 0.f };
			*v++ = QUI_AIM_V_CLR;
		}
	}

	/* axes */
	*v++ = (float3_t) { QUI_AIM_O, 0.f, 0.f };
	*v++ = QUI_AIM_X_CLR;
	*v++ = (float3_t) { QUI_AIM_O + QUI_AIM_L, 0.f, 0.f };
	*v++ = QUI_AIM_X_CLR;
	*v++ = (float3_t) { -QUI_AIM_O, 0.f, 0.f };
	*v++ = QUI_AIM_X_CLR;
	*v++ = (float3_t) { -QUI_AIM_O - QUI_AIM_L, 0.f, 0.f };
	*v++ = QUI_AIM_X_CLR;

	*v++ = (float3_t) { 0.f, QUI_AIM_O, 0.f };
	*v++ = QUI_AIM_Y_CLR;
	*v++ = (float3_t) { 0.f, QUI_AIM_O + QUI_AIM_L, 0.f };
	*v++ = QUI_AIM_Y_CLR;
	*v++ = (float3_t) { 0.f, -QUI_AIM_O, 0.f };
	*v++ = QUI_AIM_Y_CLR;
	*v++ = (float3_t) { 0.f, -QUI_AIM_O - QUI_AIM_L, 0.f };
	*v++ = QUI_AIM_Y_CLR;

	*v++ = (float3_t) { 0.f, 0.f, QUI_AIM_O };
	*v++ = QUI_AIM_Z_CLR;
	*v++ = (float3_t) { 0.f, 0.f, QUI_AIM_O + QUI_AIM_L };
	*v++ = QUI_AIM_Z_CLR;
	*v++ = (float3_t) { 0.f, 0.f, -QUI_AIM_O };
	*v++ = QUI_AIM_Z_CLR;
	*v++ = (float3_t) { 0.f, 0.f, -QUI_AIM_O - QUI_AIM_L };
	*v++ = QUI_AIM_Z_CLR;

	glUnmapNamedBuffer(qui_aim_vbo);

	glGenVertexArrays(1, &qui_aim_vao);

	if (!qui_aim_vao)
		goto ouch;

	glBindVertexArray(qui_aim_vao);
	glBindBuffer(GL_ARRAY_BUFFER, qui_aim_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, 0, 24, 0);
	glEnableVertexAttribArray(0); 
	glVertexAttribPointer(1, 3, GL_FLOAT, 0, 24, (char*)12);
	glEnableVertexAttribArray(1); 
	glBindVertexArray(0);

	return 0;
ouch:
	qui_aim_rm();

	return -1;
}

int qui_aim_rm() {
	if (qui_aim_vbo)
		glDeleteBuffers(1, &qui_aim_vbo);

	if (qui_aim_vao)
		glDeleteVertexArrays(1, &qui_aim_vao);

	return 0;
}

int qui_aim(float3_t *p, float3_t *s, float3_t *t, float3_t *n) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t M = {
		s->x, s->y, s->z, 0.f,
		t->x, t->y, t->z, 0.f,
		n->x, n->y, n->z, 0.f,
		p->x, p->y, p->z, 1.f
	};
	float44_t VM = mul_float44(M, V);
	float fVM = frobenius_float33(float33_float44(VM)) / sqrt(3);
	float44_t B = {
		fVM, 0.f, 0.f, 0.f,
		0.f, fVM, 0.f, 0.f,
		0.f, 0.f, fVM, 0.f,
		VM.m30, VM.m31, VM.m32, 1.f
	};
	float44_t S = qui_mtrx_top(QUI_MTRX_S);
	float44_t PV = mul_float44(V, P);
	float44_t PVM = mul_float44(S, mul_float44(M, PV));
	float44_t PB = mul_float44(S,mul_float44(B, P));

	glUseProgram(qui_shdr_po);
	glBindVertexArray(qui_aim_vao);
	glLineWidth(qui_aim_lnwdth);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(0);
	glPolygonOffset(0, -16);
	glEnable(GL_POLYGON_OFFSET_LINE);

	glUniformMatrix4fv(qui_shdr_M, 1, 0, &PB.m[0][0]);
	glDrawArrays(GL_LINES, 0, QUI_AIM_N);

	glUniformMatrix4fv(qui_shdr_M, 1, 0, &PVM.m[0][0]);
	glDrawArrays(GL_LINES, QUI_AIM_N, QUI_AIM_M);

	glDisable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(0, 0);
	glDepthMask(1);
	glDisable(GL_DEPTH_TEST);
	glLineWidth(1);

	if (qui_in.prss & QUI_IN_RALT) {
		if (qui_in.prss & QUI_IN_LMB) {
			int k, typ = 0;
 
			float detPV = det_float44(PV);
                        float44_t iPV = invert_float44(PV, detPV);
                        float3_t o = float3_float4(cotransform_float44(iPV, (float4_t){ qui_in.p.x, qui_in.p.y, -1, 1 }));
                        float3_t r = m_float3(cotransform_float44(iPV, (float4_t){ 0, 0, 1, 0 }));

			k = qui_qr_find(QUI_QR_FND_BST, o, qui_aim_R, &typ);
			if (0 <= k) {
				*p = qui_qr[typ][k + qui_qr_s[typ] - 1];
				switch (typ) {
				case QUI_QR_TYP_PNTS:
					*s = (float3_t) { 1, 0, 0 };
					*t = (float3_t) { 0, 1, 0 };
					*n = (float3_t) { 0, 0, 1 };
					break;
				case QUI_QR_TYP_LNS:
					*s = normal_float3(sub_float3(qui_qr[typ][k], qui_qr[typ][k+1]));
					*t = normal_float3(cross_float3(*s, r));
					*n = normal_float3(cross_float3(*t, *s));
					break;
				case QUI_QR_TYP_TRIS:
					*n = normal_float3(cross_float3(sub_float3(qui_qr[typ][k], qui_qr[typ][k+1]), sub_float3(qui_qr[typ][k+1], qui_qr[typ][k+2])));
					*s = normal_float3(cross_float3(*n, r));
					*t = normal_float3(cross_float3(*n, *s));
					break;
				default:
					return 0;
				}
				return 1;
			}
		}

		if (qui_in.s) {
			*p = add_float3(*p, scale_float3(*n, qui_in.s / fVM * 0.03125));
		}
	}

	return 0;
}

#endif /* QUI_IMPL */
#endif /* QUI_AIM_H */
