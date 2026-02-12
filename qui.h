enum {
	QUI_FLGS_AA = 0x1	/* Setup for old-school antialiasing, affects colors. */
};

extern int qui_flgs;

int qui_mk();
int qui_rm();

/* LIB */

#include "qui_mtrx.h"
#include "qui_def.h"
#include "qui_in.h"
#include "qui_shdr.h"
#include "qui_util.h"
#include "qui_strm.h"
#include "qui_fnt.h"
#include "qui_txt.h"
#include "qui_ngon.h"
#include "qui_val.h"
#include "qui_man.h"

#ifdef QUI_IMPL

int qui_flgs;

int qui_strm_mk();
int qui_strm_rm();

int qui_man_mk();
int qui_man_rm();

int qui_mk() {
	int r = 0;

	r |= qui_shdr_mk();
	r |= qui_man_mk();
	r |= qui_strm_mk();


	return r;
}

int qui_rm() {
	int r = 0;

	r |= qui_shdr_rm();
	r |= qui_man_rm();
	r |= qui_strm_rm();

	return r;
}

#endif /* QUI_IMPL */