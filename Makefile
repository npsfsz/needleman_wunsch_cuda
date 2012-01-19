ifdef DEBUG
	FLAGS=-DDEBUG=1 --debug
else
	FLAGS=-O3
endif

UTILITY_LIB_PATH=../../libs/utility_lib
STRING_BUF_PATH=../../libs/string_buffer
SCORING_PATH=../../alignment_scoring

# Add compile time
FLAGS := $(FLAGS) -DCOMPILE_TIME='"$(shell date)"'

all:
#	gcc -o needleman_wunsch $(FLAGS) -Wall -lz nw_cmdline.c nw_load_scores.c needleman_wunsch.c utility_lib.c string_buffer.c
	gcc -o needleman_wunsch $(FLAGS) -Wall -lz -I $(UTILITY_LIB_PATH) -I $(STRING_BUF_PATH) -I $(SCORING_PATH) nw_cmdline.c nw_load_scores.c needleman_wunsch.c $(UTILITY_LIB_PATH)/utility_lib.c $(STRING_BUF_PATH)/string_buffer.c

clean:
	if test -e needleman_wunsch; then rm needleman_wunsch; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
