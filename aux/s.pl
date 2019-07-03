# This file is used to generate the single-precision, real version
# of all the typed files/routines

while (<>) {
    s/\b_xnrm2\b/snrm2/g;
    s/\b_qrm_transp\b/'t'/g;
    s/\b_qrm_data_c\b/float/g;
    s/\b_qrm_data_fc\b/real(c_float)/g;
    s/\b_qrm_real_c\b/float/g;
    s/\b_qrm_real_fc\b/real(c_float)/g;
    s/\b_qrm_one_c\b/(float)1.0/g;
    s/\b_qrm_zero_c\b/(float)0.0/g;
    s/\b_qrm_sprec\b/real(kind(1.e0))/g;
    s/\b_qrm_dprec\b/real(kind(1.d0))/g;
    s/\b_qrm_cprec\b/complex(kind((1.e0,1.e0)))/g;
    s/\b_qrm_zprec\b/complex(kind((1.d0,1.d0)))/g;
    s/\b_qrm_sizeof_i\b/4/g;
    s/\b_qrm_sizeof_i8\b/8/g;
    s/\b_qrm_sizeof_s\b/4/g;
    s/\b_qrm_sizeof_d\b/8/g;
    s/\b_qrm_sizeof_c\b/8/g;
    s/\b_qrm_sizeof_z\b/16/g;
    s/\b_qrm_zero\b/0.e0/g;
    s/\b_qrm_mone\b/-1.e0/g;
    s/\b_qrm_one\b/1.e0/g;
    s/\b_qrm_two\b/2.e0/g;
    s/\b_qrm_rzero\b/0.e0/g;
    s/\b_qrm_rone\b/1.e0/g;
    s/\b_qrm_rtwo\b/2.e0/g;
    s/\b_qrm_sizeof_data\b/4/g;
    s/\b_qrm_real\b/real(kind(1.e0))/g;
    s/\b_qrm_data\b/real(kind(1.e0))/g;
    s/\b_qrm/sqrm/g;
    s/_conjg//g;
    s/\b_rx/s/g;
    s/\b_x/s/g;
    print;
}
