int qui_ngon(int n, float2_t p[], float44_t M, float4_t c) {
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

	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_ALWAYS, 0, 1);
	glColorMask(0, 0, 0, 0);

	/* clean mask */
	glStencilOp(GL_ZERO, GL_ZERO, GL_ZERO);
	glDrawArrays(GL_TRIANGLE_FAN, b, n);

	/* draw mask */
	glStencilOp(GL_INVERT, GL_INVERT, GL_INVERT);
	glDrawArrays(GL_TRIANGLE_FAN, b, n);

	/* draw when mask succeeds */
    glStencilFunc(GL_EQUAL, 1, 1);
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
    glColorMask(1, 1, 1, 1);

	glDrawArrays(GL_TRIANGLE_FAN, b, n);

	/* cleanup */
    glDisable(GL_STENCIL_TEST);

	return 0;    
}

int qui_ngon_strk(int n, float2_t p[], float44_t M, float4_t c) {
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

	glDrawArrays(GL_LINE_LOOP, b, n);

}