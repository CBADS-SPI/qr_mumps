all: setup lib examples

setup:
	(if test ! -d $(BUILD) ; then mkdir $(BUILD); fi)
	(if test ! -d $(BUILD)/src ; then mkdir $(BUILD)/src; fi)
	(if test ! -d $(BUILD)/lib ; then mkdir $(BUILD)/lib; fi)
	(if test ! -d $(BUILD)/include ; then mkdir $(BUILD)/include; fi)
	(if test ! -d $(BUILD)/examples ; then mkdir $(BUILD)/examples; fi)
	(if test ! -d $(BUILD)/testing ; then mkdir $(BUILD)/testing; fi)
	(cp makeincs/Make.inc.$(PLAT) $(BUILD)/Make.inc)
	(cp auxiliary/Makefile.build $(BUILD)/Makefile)
	(cp auxiliary/Makefile.rules.inc $(BUILD)/src)
	(cp auxiliary/Makefile.src $(BUILD)/src/Makefile)
	(cp auxiliary/Makefile.examples $(BUILD)/examples/Makefile)
	(cp auxiliary/Makefile.testing $(BUILD)/testing/Makefile)


.PHONY: lib examples

lib: setup
	(cd $(BUILD); $(MAKE) lib)

examples: lib
	(cd $(BUILD); $(MAKE) examples)


