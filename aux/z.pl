# This file is used to generate the double-precision, complex version
# of all the typed files/routines

while (<>) {
    s/\b_xnrm2\b/dznrm2/g;
    s/\b_qrm_transp\b/'c'/g;
    s/\b_qrm_data_c\b/double _Complex/g;
    s/\b_qrm_data_fc\b/complex(c_double_complex)/g;
    s/\b_qrm_real_c\b/double/g;
    s/\b_qrm_real_fc\b/real(c_double)/g;
    s/\b_qrm_one_c\b/(double _Complex)1.0 + (double _Complex)(0.0*_Complex_I)/g;
    s/\b_qrm_zero_c\b/(double _Complex)0.0 + (double _Complex)(0.0*_Complex_I)/g;
    s/\b_qrm_sprec\b/real(kind(1.e0))/g;
    s/\b_qrm_dprec\b/real(kind(1.d0))/g;
    s/\b_qrm_cprec\b/complex(kind(1.e0))/g;
    s/\b_qrm_zprec\b/complex(kind(1.d0))/g;
    s/\b_qrm_sizeof_i\b/4/g;
    s/\b_qrm_sizeof_i8\b/8/g;
    s/\b_qrm_sizeof_s\b/4/g;
    s/\b_qrm_sizeof_d\b/8/g;
    s/\b_qrm_sizeof_c\b/8/g;
    s/\b_qrm_sizeof_z\b/16/g;
    s/\b_qrm_zero\b/cmplx(0.d0,0.d0, kind(1.d0))/g;
    s/\b_qrm_one\b/cmplx(1.d0,0.d0, kind(1.d0))/g;
    s/\b_qrm_mone\b/cmplx(-1.d0,0.d0, kind(1.d0))/g;
    s/\b_qrm_two\b/cmplx(2.d0,0.d0, kind(1.d0))/g;
    s/\b_qrm_rzero\b/0.d0/g;
    s/\b_qrm_rone\b/1.d0/g;
    s/\b_qrm_rtwo\b/2.d0/g;
    s/\b_qrm_sizeof_data\b/16/g;
    s/\b_qrm_real\b/real(kind(1.d0))/g;
    s/\b_qrm_data\b/complex(kind(1.d0))/g;
    s/\b_qrm/zqrm/g;
    s/_conjg/conjg/g;
    s/\b_rx/dz/g;
    s/\b_x/z/g;
    print;
}
