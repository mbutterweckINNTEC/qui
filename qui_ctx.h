struct qui_fnt;

enum {
	QUI_FLGS_AA = 0x1	/* Setup for old-school antialiasing, affects colors. */
};

struct qui_ctx {
	/* glsl program */
	int po;

	/* uniform locations */
	int M;

	/* manipulator objects */
	int man_bo, man_vao;

	/* font */
	struct qui_fnt *fnt;

	/* general */
	int flgs;
};

int qui_man_ctx_mk(int *bo, int *vao);

int qui_ctx_mk(struct qui_ctx *qc) {
	if (!qc)
		return -1;

	qc->po = shdr_mk(qui_vsh, qui_fsh);

	if (!qc->po)
		return -1;

	qc->M = glGetUniformLocation(qc->po, "M");

	if (qc->M == -1)
		return -1;

	qui_man_ctx_mk(&qc->man_bo, &qc->man_vao);
	
	if (qc->man_vao == -1 || qc->man_bo)
		return -1;

	return 0;
}
