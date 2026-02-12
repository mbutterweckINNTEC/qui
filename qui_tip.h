
extern float4_t qui_tip_clr;

int qui_tip(char *tip);
int qui_tip2(char *tip0, char *tip1);

#ifdef QUI_IMPL

float4_t qui_tip_clr = {0, 0, 0, 1};
float const qui_tip_scl = 0.03125;

int qui_tip(char *tip) {
	int ret = 0;

	float44_t T = {
		qui_tggl_scl, 0, 0, 0,
		0, qui_tggl_scl, 0, 0,
		0, 0, 1, 0,
		-1,-0.925, 0, 1
	};

	qui_mtrx_psh(QUI_MTRX_V, identity_sc);
	qui_txt(tip, T, qui_tip_clr);
	qui_mtrx_pop(QUI_MTRX_V);

	return ret;
}

int qui_tip2(char *tip0, char *tip1) {
	int ret = 0;

	float44_t T0 = {
		qui_tggl_scl, 0, 0, 0,
		0, qui_tggl_scl, 0, 0,
		0, 0, 1, 0,
		-1,-0.925, 0, 1
	};

	float44_t T1 = {
		qui_tggl_scl, 0, 0, 0,
		0, qui_tggl_scl, 0, 0,
		0, 0, 1, 0,
		-1,-0.95, 0, 1
	};
	qui_mtrx_psh(QUI_MTRX_V, identity_sc);
	qui_txt(tip0, T0, qui_tip_clr);
	qui_txt(tip1, T1, qui_tip_clr);
	qui_mtrx_pop(QUI_MTRX_V);

	return ret;
}
#endif /* QUI_IMPL */