
TARGET = aofb
SRCS = aofb.c

ARCH ?= 
CROSS_COMPILE ?= #arm-xilinx-linux-gnueabi-

CC = $(CROSS_COMPILE)gcc
CFLAGS = -O3 -g

OBJS = $(SRCS:.c=.o)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $@

clean:
	rm -f $(OBJS) $(TARGET)
