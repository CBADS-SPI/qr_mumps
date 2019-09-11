all: setup lib examples

setup:
	(if test ! -d $(BUILD) ; then mkdir $(BUILD); fi)
	(if test ! -d $(BUILD)/src ; then mkdir $(BUILD)/src; fi)
	(if test ! -d $(BUILD)/lib ; then mkdir $(BUILD)/lib; fi)
	(if test ! -d $(BUILD)/include ; then mkdir $(BUILD)/include; fi)
	(if test ! -d $(BUILD)/examples ; then mkdir $(BUILD)/examples; fi)
	(if test ! -d $(BUILD)/testing ; then mkdir $(BUILD)/testing; fi)
	(cp makeincs/Make.inc.$(PLAT) $(BUILD)/Make.inc)
	(cp aux/Makefile.build $(BUILD)/Makefile)
	(cp aux/Makefile.rules.inc $(BUILD)/src)
	(cp aux/Makefile.src $(BUILD)/src/Makefile)
	(cp aux/Makefile.examples $(BUILD)/examples/Makefile)
	(cp aux/Makefile.testing $(BUILD)/testing/Makefile)


.PHONY: lib examples

lib: setup
	(cd $(BUILD); $(MAKE) lib)

examples: lib
	(cd $(BUILD); $(MAKE) examples)


