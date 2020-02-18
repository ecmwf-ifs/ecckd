noinst_LIBRARIES = libecckd.a
libecckd_a_SOURCES = \
	adept_scalar.h \
	constants.h \
	heating_rate.h \
	planck_function.h \
	planck_function.cpp \
	radiative_transfer_lw.h \
	radiative_transfer_lw.cpp \
	calc_cost_function_lw.h \
	calc_cost_function_lw.cpp \
	read_spectrum.h \
	read_spectrum.cpp \
	read_merged_spectrum.h \
	read_merged_spectrum.cpp \
	write_order.h \
	write_order.cpp \
	ckd_model.h \
	ckd_model.cpp \
	average_optical_depth.h \
	average_optical_depth.cpp \
	lbl_fluxes.h \
	lbl_fluxes.cpp \
	solve_lbfgs_lw.h \
	solve_lbfgs_lw.cpp \
	single_gas_data.h \
	single_gas_data.cpp

bin_PROGRAMS = reorder_spectrum_lw find_g_points_lw merge_spectra create_look_up_table \
	lw_ckd lw_spectra optimize_lut

reorder_spectrum_lw_SOURCES = reorder_spectrum_lw.cpp

find_g_points_lw_SOURCES = find_g_points_lw.cpp

merge_spectra_SOURCES = merge_spectra.cpp

create_look_up_table_SOURCES = create_look_up_table.cpp

lw_ckd_SOURCES = lw_ckd.cpp

lw_spectra_SOURCES = lw_spectra.cpp

optimize_lut_SOURCES = optimize_lut.cpp

LDADD = ./libecckd.a ../tools/libtools.a ../lbfgs/liblbfgs.a
