noinst_LIBRARIES = libecckd.a
libecckd_a_SOURCES = \
	adept_scalar.h \
	constants.h \
	heating_rate.h \
	planck_function.h \
	planck_function.cpp \
	radiative_transfer_lw.h \
	radiative_transfer_lw.cpp \
	radiative_transfer_sw.h \
	radiative_transfer_sw.cpp \
	calc_cost_function_lw.h \
	calc_cost_function_lw.cpp \
	calc_cost_function_sw.h \
	calc_cost_function_sw.cpp \
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
	solve_lbfgs.h \
	solve_lbfgs.cpp \
	single_gas_data.h \
	single_gas_data.cpp \
	equipartition.h \
	equipartition.cpp \
	write_standard_attributes.h \
	rayleigh_scattering.h \
	Timer.h \
	cumsum.h

bin_PROGRAMS = reorder_spectrum find_g_points merge_spectra \
	create_look_up_table run_ckd lw_spectra optimize_lut \
	test_equipartition scale_lut reorder_cloud_spectrum

reorder_spectrum_SOURCES = reorder_spectrum.cpp

reorder_cloud_spectrum_SOURCES = reorder_cloud_spectrum.cpp

find_g_points_SOURCES = find_g_points.cpp

merge_spectra_SOURCES = merge_spectra.cpp

create_look_up_table_SOURCES = create_look_up_table.cpp

run_ckd_SOURCES = run_ckd.cpp

lw_spectra_SOURCES = lw_spectra.cpp

optimize_lut_SOURCES = optimize_lut.cpp

test_equipartition_SOURCES = test_equipartition.cpp

scale_lut_SOURCES = scale_lut.cpp

LDADD = ./libecckd.a ../tools/libtools.a ../lbfgs/liblbfgs.a
