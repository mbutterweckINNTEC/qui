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
	QUI_IN_BCK = 0x20000,
	QUI_IN_RET = 0x40000,
	QUI_IN_RMB = 0x80000,
	QUI_IN_LALT = 0x100000,
	QUI_IN_RALT = 0x200000,
	QUI_IN_LSHFT = 0x400000,
	QUI_IN_RSHFT = 0x800000,
	QUI_IN_LCTRL = 0x1000000,
	QUI_IN_RCTRL = 0x2000000,

	QUI_IN_NUM = 0x3ff0,

	QUI_IN_ALL = ~0
};

struct qui_in {
	int prss;
	int rls;

	float2_t p;
	float2_t d;

	float s;
};

extern struct qui_in qui_in;

int qui_in_prss(int bttn);
int qui_in_rls(int bttn);
int qui_in_mv(float2_t p);
int qui_in_scrll(float scrll);
int qui_in_nxt();

#ifdef QUI_IMPL

struct qui_in qui_in;

int qui_in_prss(int bttn) {
	qui_in.prss |= bttn;
	qui_in.rls &=~ bttn;

	return 0;
}

int qui_in_rls(int bttn) {
	qui_in.rls |= qui_in.prss & bttn;
	qui_in.prss &=~ bttn;

	return 0;
}

int qui_in_mv(float2_t p) {
	qui_in.d = add_float2(qui_in.d, sub_float2(p, qui_in.p));
	qui_in.p = p;

	return 0;
}

int qui_in_scrll(float scrll) {
	qui_in.s += scrll;

	return 0;
}

int qui_in_nxt() {
	qui_in.rls = 0;
	qui_in.d = (float2_t){ 0.f, 0.f };
	qui_in.s = 0.f;

	return 0;
}

#endif