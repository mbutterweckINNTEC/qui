#ifndef QUI_LN_H
#define QUI_LN_H
#define QUI_BZR_F 8

int qui_dots(int n, float2_t p[], float44_t M, float4_t c);

/* one long polyline */
int qui_plln(int n, float2_t p[], float44_t M, float4_t c);

/* separate 2-point lines */
int qui_lns(int n, float2_t p[], float44_t M, float4_t c);

/* bezier curve of order n */
int qui_bzr(int n, float2_t p[], float44_t M, float4_t c);

#ifdef QUI_IMPL

int qui_ln_(int n, float2_t p[], float44_t M, float4_t c, int mod) {
	float2_t *v;
	int s = n * sizeof(float2_t);
	int b;

	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t PVM = mul_float44(M, mul_float44(V, P));

	if (!p || n < 3)
		return -1;

	/* upload */
	v = qui_strm_map(s, &b);

	if (NULL == v) {
		return -1;
	}

	b /= sizeof(float2_t);

	memcpy(v, p, s);

	glUnmapNamedBuffer(qui_strm_vbo);

	/* draw setup */
	glUseProgram(qui_shdr_po);
	glUniformMatrix4fv(qui_shdr_M, 1, 0, &PVM.m[0][0]);
	glBindVertexArray(qui_strm_vao);
	glVertexAttrib4fv(1, (GLfloat*)&c);

	glDrawArrays(mod, b, n);

	return 0;
}

int qui_dots(int n, float2_t p[], float44_t M, float4_t c) {
	return qui_ln_(n, p, M, c, GL_POINTS);
}

int qui_plln(int n, float2_t p[], float44_t M, float4_t c) {
	return qui_ln_(n, p, M, c, GL_LINE_STRIP);
}

int qui_lns(int n, float2_t p[], float44_t M, float4_t c) {
	if (n & 1)
		return -1;

	return qui_ln_(n, p, M, c, GL_LINES);
}

static inline float nwtn(int n, int k) {
	float ns = 1.0;

	for (int i = 1; i <= k; ++i)
		ns *= (float)(n - i + 1) / (float)i;

	return ns;
}

static inline float2_t bzr(double t, int n, float2_t p[]) {
	float2_t b = { 0.f, 0.f };
	int m = n - 1;

	for (int i = 0; i < n; ++i)
		b = add_float2(b, scale_float2(p[i], nwtn(m, i) * pow(t, i) * pow(1.f - t, m - i)));

	return b;
}

int qui_bzr(int n, float2_t p[], float44_t M, float4_t c) {
	float2_t *v;
	int nt = n * QUI_BZR_F;
	int s = nt * sizeof(float2_t);
	int b;

	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float44_t V = qui_mtrx_top(QUI_MTRX_V);
	float44_t PVM = mul_float44(M, mul_float44(V, P));

	if (!p || n < 3)
		return -1;

	/* upload */
	v = qui_strm_map(s, &b);

	if (NULL == v)
		return -1;

	b /= sizeof(float2_t);

	for (int i = 0; i < nt; ++i) {
		v[i] = bzr((float)i / (float)(nt - 1), n, p);
//		printf("%f, ",(float)i / (float)(nt - 1));
		printf("[%f,%f],", v[i].x, v[i].y);
	}
	putchar('\n');

	glUnmapNamedBuffer(qui_strm_vbo);

	/* draw setup */
	glUseProgram(qui_shdr_po);
	glUniformMatrix4fv(qui_shdr_M, 1, 0, &PVM.m[0][0]);
	glBindVertexArray(qui_strm_vao);
	glVertexAttrib4fv(1, (GLfloat*)&c);

	glDrawArrays(GL_LINE_STRIP, b, nt);

	return 0;
}

#endif /* QUI_IMPL */
#endif /* QUI_LN_H */
