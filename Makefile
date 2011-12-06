ifdef DEBUG
	FLAGS=-DDEBUG=1 --debug
endif

all:
	gcc -o needleman_wunsch $(FLAGS) -O3 -Wall -lz nw_cmdline.c nw_load_scores.c needleman_wunsch.c utility_lib.c string_buffer.c

clean:
	if test -e needleman_wunsch; then rm needleman_wunsch; fi
	
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
