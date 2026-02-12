CC = gcc
CFLGS = -ggdb -O0
LFLGS = -lglfw -lGLEW -lEGL -lGL -lGLU -lOpenGL -lm

qui_demo: qui_demo.c qui.h qui_bttn.h qui_def.h qui_fnt.h qui_in.h qui_man.h qui_mtrx.h qui_ngon.h qui_shdr.h qui_strm.h qui_tggl.h qui_txt.h qui_util.h qui_val.h
	$(CC) $(CFLGS) $< -o $@ $(LFLGS)

