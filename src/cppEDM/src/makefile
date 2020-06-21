
## CXX = g++ # CXX defined in Makevars from R environment variables

## CXXFLAGS defined in Makevars from R environment variables
## CFLAGS += -std=c++11 -DCCM_THREADED -DMULTIVIEW_VALUES_OVERLOAD -O3 -fPIC
CFLAGS = $(CXXFLAGS) -DCCM_THREADED -DMULTIVIEW_VALUES_OVERLOAD

OBJ = Common.o AuxFunc.o DateTimeUtil.o Parameter.o Embed.o Interface.o\
	Neighbors.o Simplex.o Eval.o CCM.o Multiview.o SMap.o 

LIB = libEDM.a

all:	$(LIB)
	$(AR) -rcs $(LIB) $(OBJ)  # AR passed from Makevars
	cp $(LIB) ../lib/

clean:
	rm -f $(OBJ) $(LIB)

distclean:
	rm -f $(OBJ) $(LIB) ../lib/$(LIB) *~ *.bak *.csv

$(LIB): $(OBJ)

AuxFunc.o: AuxFunc.cc
	$(CXX) -c AuxFunc.cc $(CFLAGS)

DateTimeUtil.o: DateTimeUtil.cc
	$(CXX) -c DateTimeUtil.cc $(CFLAGS)

Common.o: Common.cc
	$(CXX) -c Common.cc $(CFLAGS)

Parameter.o: Parameter.cc
	$(CXX) -c Parameter.cc $(CFLAGS)

Embed.o: Embed.cc
	$(CXX) -c Embed.cc $(CFLAGS)

Interface.o: Interface.cc
	$(CXX) -c Interface.cc $(CFLAGS)

Neighbors.o: Neighbors.cc
	$(CXX) -c Neighbors.cc $(CFLAGS)

Simplex.o: Simplex.cc
	$(CXX) -c Simplex.cc $(CFLAGS)

Eval.o: Eval.cc
	$(CXX) -c Eval.cc $(CFLAGS)

CCM.o: CCM.cc
	$(CXX) -c CCM.cc $(CFLAGS)

Multiview.o: Multiview.cc
	$(CXX) -c Multiview.cc $(CFLAGS)

SMap.o: SMap.cc
	$(CXX) -c SMap.cc $(CFLAGS)

SRCS = `echo ${OBJ} | sed -e 's/.o /.cc /g'`
depend:
	@echo ${SRCS}
	makedepend -Y $(SRCS)
# DO NOT DELETE

Common.o: Common.h DataFrame.h
AuxFunc.o: AuxFunc.h Common.h DataFrame.h Neighbors.h Parameter.h Version.h
AuxFunc.o: Embed.h DateTime.h
DateTimeUtil.o: DateTime.h
Parameter.o: Parameter.h Common.h DataFrame.h Version.h
Embed.o: Embed.h Common.h DataFrame.h Parameter.h Version.h
Interface.o: Common.h DataFrame.h
Neighbors.o: Neighbors.h Common.h DataFrame.h Parameter.h Version.h
Simplex.o: Common.h DataFrame.h Parameter.h Version.h Neighbors.h Embed.h
Simplex.o: AuxFunc.h
Eval.o: Common.h DataFrame.h
CCM.o: Common.h DataFrame.h Embed.h Parameter.h Version.h AuxFunc.h
CCM.o: Neighbors.h
Multiview.o: Common.h DataFrame.h AuxFunc.h Neighbors.h Parameter.h Version.h
Multiview.o: Embed.h
