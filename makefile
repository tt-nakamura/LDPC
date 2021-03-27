NTL = -lntl -lgmp -L/usr/local/lib
LDPC = LDPCCode.o indexx.o AddNoise.o
OTHERS = ConvolCode.o BCHCode.o RSCode.o GF2Xlib.o

example: example.o $(LDPC) $(OTHERS)
	g++ example.o $(LDPC) $(OTHERS) $(NTL)

fig1: fig1.o $(LDPC)
	g++ fig1.o $(LDPC) $(NTL)
fig2: fig2.o $(LDPC)
	g++ fig2.o $(LDPC) $(NTL)
fig3: fig3.o $(LDPC)
	g++ fig3.o $(LDPC) $(NTL)
fig4: fig4.o $(LDPC) $(OTHERS)
	g++ fig4.o $(LDPC) $(OTHERS) $(NTL)
table1: table1.o $(LDPC) $(OTHERS)
	g++ table1.o $(LDPC) $(OTHERS) $(NTL)