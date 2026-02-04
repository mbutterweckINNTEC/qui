int qui_strm_mk(int *strm_vbo, int *strm_vao) {
	glGenVertexArrays(1, strm_vao);

	if (0 == *strm_vao)
		return -1;

	glCreateBuffers(1, strm_vbo);

	if (0 == *strm_vbo)
		return -1;

	glNamedBufferStorage(*strm_vbo, QUI_STRM_SZ, NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

	glBindVertexArray(*strm_vao);
	glBindBuffer(GL_ARRAY_BUFFER, *strm_vbo);
	glVertexAttribPointer(0, 2, GL_FLOAT, 0, 8, 0);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);

	return 0;
}

void *qui_strm_map(int strm_vbo, int *strm_n, int s, int *b) {
	void *v = NULL;

	if (QUI_STRM_SZ < *strm_n + s) {
		*strm_n = s;
		*b = 0;
		v = glMapNamedBufferRange(strm_vbo, 0, s, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT | GL_MAP_UNSYNCHRONIZED_BIT);
	} else {
		*b = *strm_n;
		*strm_n += s;
		v = glMapNamedBufferRange(strm_vbo, *b, s, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_RANGE_BIT | GL_MAP_UNSYNCHRONIZED_BIT);
	}

	return v;
}