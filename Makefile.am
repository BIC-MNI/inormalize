bin_PROGRAMS = inormalize
inormalize_SOURCES = InormalizeArgs.cc inormalize.cc

bin_SCRIPTS = nu_correct_norm normalize_mri lgmask headmask
EXTRA_DIST = nu_correct_norm.pl normalize_mri.pl lgmask.pl headmask.pl

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

CLEANFILES = nu_correct_norm normalize_mri lgmask headmask
