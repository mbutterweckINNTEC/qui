/* input facilities of qui */

enum {
	QUI_IN_LMB = 0x1,
	QUI_IN_X = 0x2,
	QUI_IN_Y = 0x4,
	QUI_IN_Z = 0x8,
	QUI_IN_0 = 0x10,
	QUI_IN_1 = 0x20,
	QUI_IN_2 = 0x40,
	QUI_IN_3 = 0x80,
	QUI_IN_4 = 0x100,
	QUI_IN_5 = 0x200,
	QUI_IN_6 = 0x400,
	QUI_IN_7 = 0x800,
	QUI_IN_8 = 0x1000,
	QUI_IN_9 = 0x2000,
	QUI_IN_DOT = 0x4000,
	QUI_IN_MINUS = 0x8000,
	QUI_IN_ESC = 0x10000,

	QUI_IN_ALL = ~0
};

struct qui_in {
	int prss;
	int rls;

	float2_t p;
	float2_t d;

	float s;
};

static int qui_in_prss(struct qui_in *qi, int bttn);
static int qui_in_rls(struct qui_in *qin, int bttn);
static int qui_in_mv(struct qui_in *qi, float2_t p);
static int qui_in_scrll(struct qui_in *qi, float scrll);
static int qui_in_nxt(struct qui_in *qi);

static int qui_in_prss(struct qui_in *qi, int bttn) {
	if (!qi)
		return -1;

	qi->prss |= bttn;
	qi->rls &=~ bttn;

	return 0;
}

static int qui_in_rls(struct qui_in *qi, int bttn) {
	if (!qi)
		return -1;

	qi->rls |= qi->prss & bttn;
	qi->prss &=~ bttn;

	return 0;
}

static int qui_in_mv(struct qui_in *qi, float2_t p) {
	if (!qi)
		return -1;

	qi->d = add_float2(qi->d, sub_float2(p, qi->p));
	qi->p = p;

	return 0;
}

static int qui_in_scrll(struct qui_in *qi, float scrll) {
	if (!qi)
		return -1;

	qi->s += scrll;

	return 0;
}

static int qui_in_nxt(struct qui_in *qi) {
	if (!qi)
		return -1;

	qi->rls = 0;
	qi->d = (float2_t){ 0.f, 0.f };
	qi->s = 0.f;

	return 0;
}
