all:
	gcc -o needleman_wunsch -Wall -lz nw_cmdline.c needleman_wunsch.c utility_lib.c string_buffer.c

clean:
	rm needleman_wunsch
