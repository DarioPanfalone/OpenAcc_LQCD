#ifndef PACKER_H
#define PACKER_H

void smartpack_gauge(float out[2*12*no_links] , const Conf *in);
void smartunpack_gauge(Conf *out, const float in[2*12*no_links]);

void smartpack_fermion(float out[6*sizeh], const Fermion *in);
void smartpack_fermion_d(float out[6*sizeh*2], const Fermion *in);
void smartunpack_fermion_d(Fermion *out, const float in[6*sizeh*2]);

void make_shift_table(int table[8*size]);

#endif
