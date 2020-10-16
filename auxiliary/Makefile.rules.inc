STARPU=$(filter %have_starpu, $(FDEFS))

TEMPLATES     = $(topdir)/src/analysis/qrm_analyse.F90 \
                $(topdir)/src/analysis/qrm_compute_graph.F90 \
                $(topdir)/src/analysis/qrm_do_ordering.F90 \
                $(topdir)/src/analysis/qrm_do_colamd.F90 \
                $(topdir)/src/analysis/qrm_elim_tree.F90 \
                $(topdir)/src/analysis/qrm_rowcount.F90 \
                $(topdir)/src/analysis/qrm_rowperm.F90 \
                $(topdir)/src/analysis/qrm_do_metis.F90 \
                $(topdir)/src/analysis/qrm_ata_graph.F90 \
                $(topdir)/src/analysis/qrm_do_scotch.F90 \
                $(topdir)/src/analysis/qrm_symbolic.F90 \
                $(topdir)/src/analysis/qrm_compute_memory.F90 \
                $(topdir)/src/factorization/qrm_factorize.F90  \
	        $(topdir)/src/factorization/qrm_init_front.F90 \
	        $(topdir)/src/factorization/qrm_factorization_core.F90 \
	        $(topdir)/src/factorization/qrm_factorization_init.F90 \
	        $(topdir)/src/factorization/qrm_do_subtree.F90 \
	        $(topdir)/src/factorization/qrm_clean_front.F90 \
	        $(topdir)/src/factorization/qrm_assemble_block.F90 \
	        $(topdir)/src/factorization/qrm_activate_front.F90 \
	        $(topdir)/src/factorization/qrm_factorization_tasks.F90 \
                $(topdir)/src/methods/qrm_least_squares.F90 \
                $(topdir)/src/methods/qrm_min_norm.F90 \
	        $(topdir)/src/methods/qrm_residual_norm.F90 \
                $(topdir)/src/methods/qrm_residual_orth.F90 \
                $(topdir)/src/methods/qrm_get_r.F90 \
	        $(topdir)/src/modules/qrm_c_interface.F90 \
	        $(topdir)/src/modules/qrm_methods_mod.F90 \
                $(topdir)/src/solve/qrm_solve_kernels.F90 \
                $(topdir)/src/solve/qrm_solve.F90 \
	        $(topdir)/src/solve/qrm_apply.F90 \
                $(topdir)/src/solve/qrm_solve_tasks.F90 \
	        $(topdir)/src/solve/qrm_solve_sing_front.F90 \
                $(topdir)/src/utils/qrm_readmat.F90 \
                $(topdir)/src/utils/qrm_matmul.F90 \
                $(topdir)/src/utils/qrm_matnrm.F90 \
                $(topdir)/src/utils/qrm_vecnrm.F90 \
	        $(topdir)/src/utils/qrm_remap_pnt.F90 \
	        $(topdir)/src/F77/qrm_geqrt.F90 \
	        $(topdir)/src/F77/qrm_gemqrt.F90 \
	        $(topdir)/src/F77/qrm_tpqrt.F90 \
	        $(topdir)/src/F77/qrm_tpmqrt.F90


ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
TEMPLATES    += $(topdir)/src/StarPU/qrm_starpu_cpu_funcs.F90 \
	        $(topdir)/src/StarPU/qrm_starpu_codelets.c
endif

COMMON        = $(topdir)/src/analysis/qrm_amalg_tree.F90 \
	        $(topdir)/src/analysis/qrm_check_cperm.F90 \
	        $(topdir)/src/analysis/qrm_prune_tree.F90 \
	        $(topdir)/src/analysis/qrm_reorder_tree.F90 \
	        $(topdir)/src/analysis/qrm_postorder.F90 \
	        $(topdir)/src/analysis/qrm_print_tree.F90 \
	        $(topdir)/src/analysis/qrm_compress_data.F90 \
	        $(topdir)/src/modules/qrm_c_comm_interface.F90 \
	        $(topdir)/src/utils/qrm_prnt_array.F90 \
	        $(topdir)/src/C/qrm_get_wtime.c \
	        $(topdir)/src/C/qrm_metis_wrap.c \
	        $(topdir)/src/C/qrm_pthread_wrap.c \
	        $(topdir)/src/C/qrm_colamd_wrap.c \
	        $(topdir)/src/C/qrm_atomic.c


ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
COMMON       += $(topdir)/src/StarPU/starpu_f_wrappers.c
endif

MODTEMPLATES  = $(topdir)/src/modules/qrm_spmat_mod.F90 \
	        $(topdir)/src/modules/qrm_mod.F90 \
	        $(topdir)/src/modules/qrm_utils_mod.F90 \
	        $(topdir)/src/modules/qrm_analysis_mod.F90 \
	        $(topdir)/src/modules/qrm_factorization_mod.F90 \
	        $(topdir)/src/modules/qrm_fdata_mod.F90 \
	        $(topdir)/src/modules/qrm_sdata_mod.F90 \
	        $(topdir)/src/modules/qrm_solve_mod.F90

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
MODTEMPLATES += $(topdir)/src/StarPU/qrm_starpu_codelets_mod.F90
endif

MODCOMMON     = $(topdir)/src/modules/qrm_memhandling_mod.F90 \
	        $(topdir)/src/modules/qrm_const_mod.F90 \
	        $(topdir)/src/modules/qrm_error_mod.F90 \
	        $(topdir)/src/modules/qrm_sort_mod.F90 \
	        $(topdir)/src/modules/qrm_string_mod.F90 \
	        $(topdir)/src/modules/qrm_adata_mod.F90 \
	        $(topdir)/src/modules/qrm_pthread_mod.F90 \
	        $(topdir)/src/modules/qrm_dscr_mod.F90 \
	        $(topdir)/src/modules/qrm_common_mod.F90

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
MODCOMMON    += $(topdir)/src/StarPU/starpu_f_mod.F90
endif

CHEADERSINC   = $(topdir)/src/include/qrm_common.h \
	        $(topdir)/src/include/qrm_mumps.h

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
CHEADERSINC  += $(topdir)/src/include/qrm_starpu.h \
	        $(topdir)/src/include/qrm_starpu_prof.h
endif



TEMPLATES   += $(MODTEMPLATES)
COMMON      += $(MODCOMMON)

COMMONSRC    = $(notdir $(COMMON))
ARITHMSRC    = $(addprefix $(arith), $(notdir $(TEMPLATES)))
CHEADERS     = $(notdir $(CHEADERSINC))


ARITHMOBJS   = $(patsubst %.F90,%.o,$(patsubst %.f,%.o,$(patsubst %.c,%.o,$(ARITHMSRC))))
COMMONOBJS   = $(patsubst %.F90,%.o,$(patsubst %.f,%.o,$(patsubst %.c,%.o,$(COMMONSRC))))  qrm_mem_mod.o


CHEADERSSUBS=$(filter %$*, $(CHEADERSINC))
.SECONDEXPANSION:
$(CHEADERS): % : $$(CHEADERSSUBS)
	cp $< .



COMMONSUBS=$(filter %$*, $(COMMON))
.SECONDEXPANSION:
$(COMMONSRC): % : $$(COMMONSUBS)
	cp $< .



ARITHMSUBS=$(filter %$*, $(TEMPLATES))
.SECONDEXPANSION:
$(ARITHMSRC): $(arith)% : $$(ARITHMSUBS)
	perl $(topdir)/auxiliary/$(arith).pl $< > $@



# Special rule for generating the memory allocation/deallocation
# routines for allocatables and pointers and for all the arithmetics
# from a template
qrm_mem_mod.F90: $(topdir)/src/modules/qrm_mem_mod_template.F90
	perl $(topdir)/auxiliary/qrm_genmemmod.pl $< > $@






# Gruesome dependencies more or less automatically generated with gfortran -MM
$(addprefix $(arith), qrm_activate_front.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_activate_front.o):\
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_analyse.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod \
	qrm_string_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_analyse.o):\
	$(addprefix $(arith), qrm_starpu_codelets_mod.mod) \
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_analysis_mod.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_dscr_mod.mod 

$(addprefix $(arith), qrm_apply.o):\
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_dscr_mod.mod \
	qrm_string_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_apply.o):\
	$(addprefix $(arith), qrm_starpu_codelets_mod.mod) \
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_assemble_block.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common_mod.mod \
	qrm_error_mod.mod 

$(addprefix $(arith), qrm_ata_graph.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_c_interface.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_methods_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h 

$(addprefix $(arith), qrm_clean_front.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_mem_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_clean_front.o):\
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_compute_graph.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_compute_memory.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_memhandling_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_do_colamd.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_do_metis.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_do_ordering.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_const_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_do_scotch.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_do_subtree.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_elim_tree.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_factorization_core.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_factorization_init.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_factorization_init.o):\
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_factorization_mod.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_dscr_mod.mod 

$(addprefix $(arith), qrm_factorization_tasks.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common_mod.mod \
	qrm_dscr_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_factorization_tasks.o):\
	$(addprefix $(arith), qrm_starpu_codelets_mod.mod) 
endif

$(addprefix $(arith), qrm_factorize.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod \
	qrm_string_mod.mod 

$(addprefix $(arith), qrm_fdata_mod.o):\
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_memhandling_mod.mod \
	qrm_mem_mod.mod \
	qrm_pthread_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_fdata_mod.o):\
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_get_r.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_init_front.o):\
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod \
	qrm_sort_mod.mod 

$(addprefix $(arith), qrm_least_squares.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_methods_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod 

$(addprefix $(arith), qrm_matmul.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_string_mod.mod 

$(addprefix $(arith), qrm_matnrm.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_string_mod.mod 

$(addprefix $(arith), qrm_methods_mod.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) 

$(addprefix $(arith), qrm_min_norm.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_methods_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_mod.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_methods_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common_mod.mod \
	qrm_const_mod.mod \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod \
	qrm_string_mod.mod 

$(addprefix $(arith), qrm_readmat.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod 

$(addprefix $(arith), qrm_residual_norm.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod 

$(addprefix $(arith), qrm_residual_orth.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_rowcount.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_rowperm.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod 

$(addprefix $(arith), qrm_sdata_mod.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_sdata_mod.o):\
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_solve.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod \
	qrm_string_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_solve.o):\
	$(addprefix $(arith), qrm_starpu_codelets_mod.mod) \
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_solve_kernels.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod \
	qrm_string_mod.mod 

$(addprefix $(arith), qrm_solve_mod.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_dscr_mod.mod 

$(addprefix $(arith), qrm_solve_sing_front.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_string_mod.mod 

$(addprefix $(arith), qrm_solve_tasks.o):\
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_dscr_mod.mod \
	qrm_string_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_solve_tasks.o):\
	$(addprefix $(arith), qrm_starpu_codelets_mod.mod) 
endif

$(addprefix $(arith), qrm_spmat_mod.o):\
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod \
	qrm_sort_mod.mod \
	qrm_string_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_spmat_mod.o):\
	starpu_f_mod.mod 
endif

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_starpu_codelets_mod.o):\
	qrm_dscr_mod.mod 
endif

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_starpu_cpu_funcs.o):\
	$(addprefix $(arith), qrm_analysis_mod.mod) \
	$(addprefix $(arith), qrm_factorization_mod.mod) \
	$(addprefix $(arith), qrm_fdata_mod.mod) \
	$(addprefix $(arith), qrm_sdata_mod.mod) \
	$(addprefix $(arith), qrm_solve_mod.mod) \
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_dscr_mod.mod \
	$(addprefix $(arith), qrm_starpu_codelets_mod.mod) \
	starpu_f_mod.mod 
endif

$(addprefix $(arith), qrm_symbolic.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	$(addprefix $(arith), qrm_utils_mod.mod) \
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_sort_mod.mod 

$(addprefix $(arith), qrm_utils_mod.o):\
	$(addprefix $(arith), qrm_spmat_mod.mod) \
	qrm_adata_mod.mod \
	qrm_common_mod.mod 

$(addprefix $(arith), qrm_vecnrm.o):\
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_string_mod.mod 

qrm_adata_mod.o:\
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
qrm_adata_mod.o:\
	starpu_f_mod.mod 
endif

qrm_amalg_tree.o:\
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

qrm_c_comm_interface.o:\
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_dscr_mod.mod \
	qrm_error_mod.mod 

qrm_check_cperm.o:\
	qrm_common.h \
	qrm_error_mod.mod \
	qrm_mem_mod.mod 

qrm_common_mod.o:\
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_const_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod \
	qrm_string_mod.mod 

qrm_compress_data.o:\
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod 

qrm_dscr_mod.o:\
	qrm_common.h \
	qrm_const_mod.mod \
	qrm_error_mod.mod \
	qrm_memhandling_mod.mod \
	qrm_pthread_mod.mod \
	qrm_string_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
qrm_dscr_mod.o:\
	starpu_f_mod.mod 
endif

qrm_error_mod.o:\
	qrm_common.h \
	qrm_const_mod.mod 

qrm_memhandling_mod.o:\
	qrm_const_mod.mod \
	qrm_error_mod.mod \
	qrm_pthread_mod.mod 

qrm_mem_mod.o:\
	qrm_const_mod.mod \
	qrm_error_mod.mod \
	qrm_memhandling_mod.mod \
	qrm_pthread_mod.mod 

qrm_postorder.o:\
	qrm_common.h \
	qrm_mem_mod.mod 

qrm_print_tree.o:\
	qrm_adata_mod.mod \
	qrm_mem_mod.mod 

qrm_prune_tree.o:\
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_error_mod.mod \
	qrm_mem_mod.mod \
	qrm_sort_mod.mod 

qrm_reorder_tree.o:\
	qrm_adata_mod.mod \
	qrm_common.h \
	qrm_common_mod.mod \
	qrm_mem_mod.mod \
	qrm_sort_mod.mod 

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
$(addprefix $(arith), qrm_starpu_codelets.o):\
	qrm_starpu.h 
endif

ifeq ($(STARPU), $(DEFINE_PREPEND)-Dhave_starpu)
qrm_starpu_sched.o:\
	qrm_starpu.h 
endif
