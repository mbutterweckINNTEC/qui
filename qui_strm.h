#define QUI_STRM_SZ 0x10000

extern int qui_strm_vbo, qui_strm_vao, qui_strm_n;

void *qui_strm_map(int s, int *offst /* out */);

int qui_strm_mk();
int qui_strm_rm();

#ifdef QUI_IMPL

int qui_strm_vbo, qui_strm_vao, qui_strm_n;

int qui_strm_mk() {
	glGenVertexArrays(1, &qui_strm_vao);

	if (0 == qui_strm_vao)
		return -1;

	glCreateBuffers(1, &qui_strm_vbo);

	if (0 == qui_strm_vbo)
		return -1;

	glNamedBufferStorage(qui_strm_vbo, QUI_STRM_SZ, NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

	glBindVertexArray(qui_strm_vao);
	glBindBuffer(GL_ARRAY_BUFFER, qui_strm_vbo);
	glVertexAttribPointer(0, 2, GL_FLOAT, 0, 8, 0);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);

	return 0;
}

int qui_strm_rm() {
	if (qui_strm_vbo)
		glDeleteBuffers(1, &qui_strm_vbo);

	if (qui_strm_vao)
		glDeleteVertexArrays(1, &qui_strm_vao);

	return 0;
}


void *qui_strm_map(int s, int *offst /* out */) {
	void *v = NULL;

	if (QUI_STRM_SZ < qui_strm_n + s) {
		qui_strm_n = s;
		*offst = 0;
		v = glMapNamedBufferRange(qui_strm_vbo, 0, s, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT | GL_MAP_UNSYNCHRONIZED_BIT);
	} else {
		*offst = qui_strm_n;
		qui_strm_n += s;
		v = glMapNamedBufferRange(qui_strm_vbo, *offst, s, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_RANGE_BIT | GL_MAP_UNSYNCHRONIZED_BIT);
	}

	return v;
}

#endif /* QUI_IMPL */