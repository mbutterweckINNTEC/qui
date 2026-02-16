extern float4_t qui_tip_clr;

void qui_tip_bgn();

void qui_tip_end(char *tip);
void qui_tip_end2(char *tip0, char *tip1);

/* PROT */

enum {
	QUI_TIP_SGNL_FCS = 0x1	/* signal focus */
};

extern int qui_tip_sgnl;
int qui_tip_msk;

#ifdef QUI_IMPL

int qui_tip_sgnl;
int qui_tip_msk;
char *qui_tip_id;
time_t qui_tip_tm;

float4_t qui_tip_clr = {0, 0, 0, 1};
float const qui_tip_scl = 0.03125;

int qui_tip_tst(char *id) {
	if (qui_tip_sgnl & QUI_TIP_SGNL_FCS) {
		qui_tip_sgnl = 0;

		if (id != qui_tip_id) {
			qui_tip_tm = time(NULL);
			 qui_tip_id = id;
		}

		if (2 < time(NULL) - qui_tip_tm)
			return 1;
	} else {
		if (id == qui_tip_id)
			qui_tip_id = NULL;
	}

	return 0;
}

void qui_tip_end(char *tip) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float detP = det_float44(P);
	float4_t p = { -0.925,-0.925 };

	if (detP) {
		float44_t iP = invert_float44(P, detP);

		p = cotransform_float44(iP, p);
	}

	float44_t T = {
		qui_tip_scl, 0, 0, 0,
		0, qui_tip_scl, 0, 0,
		0, 0, 1, 0,
		p.x, p.y, 0, 1
	};

	if (qui_tip_tst(tip)) {
		qui_mtrx_psh(QUI_MTRX_V, identity_sc);
		qui_txt(tip, T, qui_tip_clr);
		qui_mtrx_pop(QUI_MTRX_V);
	}
}

void qui_tip_end2(char *tip0, char *tip1) {
	float44_t P = qui_mtrx_top(QUI_MTRX_P);
	float detP = det_float44(P);
	float4_t pu = { -0.925,-0.925 };
	float4_t pb = { -0.925,-0.95 };

	if (detP) {
		float44_t iP = invert_float44(P, detP);

		pu = cotransform_float44(iP, pu);
		pb = pu;
		pb.y -= 0.025;
	}

	float44_t T0 = {
		qui_tip_scl, 0, 0, 0,
		0, qui_tip_scl, 0, 0,
		0, 0, 1, 0,
		pu.x, pu.y, 0, 1
	};

	float44_t T1 = {
		qui_tip_scl, 0, 0, 0,
		0, qui_tip_scl, 0, 0,
		0, 0, 1, 0,
		pb.x, pb.y, 0, 1
	};

	if (qui_tip_tst(tip0)) {
		qui_mtrx_psh(QUI_MTRX_V, identity_sc);
		qui_txt(tip0, T0, qui_tip_clr);
		qui_txt(tip1, T1, qui_tip_clr);
		qui_mtrx_pop(QUI_MTRX_V);
	}
}

void qui_tip_bgn() {
	assert(0 == qui_tip_sgnl);

	qui_tip_sgnl = 0;
	qui_tip_msk = ~0;
}

#endif /* QUI_IMPL */