CC = gcc
CFLGS = -ggdb -O0
LFLGS = -lglfw -lGLEW -lEGL -lGL -lGLU -lOpenGL -lm

qui_demo: qui_demo.c qui_ctx.h qui_def.h qui_fnt.h qui_in.h qui_man.h qui_shdr.h qui_txt.h qui_util.h
	$(CC) $(CFLGS) $< -o $@ $(LFLGS)

