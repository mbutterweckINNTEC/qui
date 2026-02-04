int qui_ngon(struct qui_ctx *qc, int n, float2_t p[], float44_t M, float4_t c) {
	float2_t *v;
	int s = n * sizeof(float2_t);
	int b;

	if (!qc || !p || n < 3)
		return -1;


	/* upload */
	v = qui_strm_map(qc->strm_vbo, &qc->strm_n, s, &b);

	if (NULL == v) {
		printf("no mapping %d\n", glGetError());
		return -1;
	}

	b /= sizeof(float2_t);

	memcpy(v, p, s);

	glUnmapNamedBuffer(qc->strm_vbo);

	/* draw setup */
	glUseProgram(qc->po);
	glUniformMatrix4fv(qc->M, 1, 0, &M.m[0][0]);
	glBindVertexArray(qc->strm_vao);
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