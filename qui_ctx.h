struct qui_fnt;
struct qui_in;

#define QUI_STRM_SZ 0x10000

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

	/* streaming buffer */
	int strm_vbo, strm_vao, strm_n;

	/* inputs */
	struct qui_in *in;

	/* general */
	int flgs;
};

int qui_man_ctx_mk(int *bo, int *vao);
int qui_strm_mk(int *strm_vbo, int *strm_vao);

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
	
	if (0 == qc->man_vao || 0 == qc->man_bo)
		return -1;

	qui_strm_mk(&qc->strm_vbo, &qc->strm_vao);

	if (0 == qc->strm_vao || 0 == qc->strm_vbo)
		return -1;

	qc->strm_n = 0;

	return 0;
}
