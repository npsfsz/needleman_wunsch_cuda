LIBS_PATH=../../libs

ifndef CC
	CC = gcc
endif

ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug -g
else
	CFLAGS := -O3
endif

UTILITY_LIB_PATH = $(LIBS_PATH)/utility_lib
STRING_BUF_PATH = $(LIBS_PATH)/string_buffer
BIT_ARRAY_PATH = $(LIBS_PATH)/bit_array
SEQ_FILE_PATH = $(LIBS_PATH)/seq_file
SCORING_PATH = $(LIBS_PATH)/alignment_scoring
SAMTOOLS_PATH = $(HOME)/bioinf/samtools-0.1.18

# Add data type for alignment scoring
CFLAGS := $(CFLAGS) -Wall -Wextra -DSCORE_TYPE='unsigned int' \
         -I$(SCORING_PATH) -I$(SEQ_FILE_PATH) -I$(UTILITY_LIB_PATH) \
         -I$(BIT_ARRAY_PATH) -I$(STRING_BUF_PATH) -I$(SAMTOOLS_PATH) -I.

LIB_INCS = -L$(SCORING_PATH) -L$(SEQ_FILE_PATH) -L$(UTILITY_LIB_PATH) \
           -L$(BIT_ARRAY_PATH) -L$(STRING_BUF_PATH) -L$(SAMTOOLS_PATH) -L.

LIB_LIST = -lseqfile -lstrbuf -lbitarr -lbam -lutil -lz

FILES = nw_cmdline.c needleman_wunsch.c $(SCORING_PATH)/*.c

all:
	$(CC) -o needleman_wunsch $(CFLAGS) $(LIB_INCS) $(FILES) $(LIB_LIST)

clean:
	rm -rf needleman_wunsch needleman_wunsch.dSYM needleman_wunsch.greg
