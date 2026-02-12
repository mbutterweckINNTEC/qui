enum {
	QUI_MTRX_P,	/* projection matrix */
	QUI_MTRX_V,	/* view matrix */

	QUI_MTRX_N
};

int qui_mtrx_psh(int typ, float44_t M);
int qui_mtrx_pop(int typ);
float44_t qui_mtrx_top(int typ);

#ifdef QUI_IMPL

#define QUI_MTRX_S 8

float44_t qui_mtrx[QUI_MTRX_N][QUI_MTRX_S];
int qui_mtrx_s[QUI_MTRX_N];

int qui_mtrx_psh(int typ, float44_t M) {
	if (typ < QUI_MTRX_N) {
		if (qui_mtrx_s[typ] < QUI_MTRX_S) {
			qui_mtrx[typ][qui_mtrx_s[typ]++] = M;
			return 0;
		}
	}

	return -1;
}
int qui_mtrx_pop(int typ) {
	if (typ < QUI_MTRX_N) {
		if (qui_mtrx_s[typ]) {
			--qui_mtrx_s[typ];
			return 0;
		}
	}

	return -1;
}

float44_t qui_mtrx_top(int typ) {
	if (typ < QUI_MTRX_N) {
		if (qui_mtrx_s[typ]) {
			return qui_mtrx[typ][qui_mtrx_s[typ] - 1];
		}
	}

	return identity_sc;
}


#endif