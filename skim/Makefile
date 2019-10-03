
override LINKFLAGS += -shared -L$(BELLE_TOP_DIR)/lib/so -L$(BELLE_TOP_DIR)/bin


override DEBUG =-O3
override COMPILEFLAGS += -I$(BELLE_TOP_DIR)/include -DHEP_SHORT_NAMES  -fpic -pipe
SRC = User_reco.cc
# If there are several .cc files to be compiled you must point them above 



%.o: %.cc %.h
	gcc -c $(COMPILEFLAGS) $(DEBUG) -MMD $(INC) -o $@ $<

User_reco.so: User_reco.o userinfo.o
	gcc $(LINKFLAGS) $(DEBUG)    $^ -o $@
clean:
	-rm *.o *.so *.d


-include $(DEP)
